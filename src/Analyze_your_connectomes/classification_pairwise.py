import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.svm import SVC
from sklearn.ensemble import RandomForestClassifier
from sklearn.neural_network import MLPClassifier  # NEW: MLP
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.cross_decomposition import PLSRegression
from sklearn.inspection import permutation_importance
from sklearn.metrics import (
    balanced_accuracy_score, f1_score, cohen_kappa_score,
    confusion_matrix, accuracy_score
)
from sklearn.base import BaseEstimator, TransformerMixin
from sklearn.preprocessing import OneHotEncoder
from sklearn.cross_decomposition import PLSRegression
from sklearn.feature_selection import SelectKBest, f_classif
from sklearn.model_selection import (
    StratifiedShuffleSplit, StratifiedKFold, GridSearchCV,
    train_test_split, PredefinedSplit  # NEW: simple CV helpers
)
from sklearn.utils import parallel_backend
from tqdm import tqdm
from imblearn.under_sampling import CondensedNearestNeighbour  # (usato se resampling=True)
from itertools import combinations
import json
import re
from sklearn.decomposition import PCA, NMF, FastICA, IncrementalPCA
from sklearn.decomposition import TruncatedSVD
from sklearn.cross_decomposition import CCA
from sklearn.neighbors import NeighborhoodComponentsAnalysis
from sklearn.preprocessing import OneHotEncoder
from sklearn.base import BaseEstimator, TransformerMixin
from imblearn.under_sampling import CondensedNearestNeighbour
import warnings
from sklearn.exceptions import ConvergenceWarning
from sklearn.linear_model import LogisticRegression
from sklearn.svm import LinearSVC
from sklearn.discriminant_analysis import LinearDiscriminantAnalysis
import os

# Read environment variables passed by the job launcher
threshold_strategy = os.environ.get('THRESHOLD_STRATEGY', 'proportional_edges')
classifier_name = os.environ.get('CLASSIFIER_NAME', 'LIN_SVM')

if threshold_strategy == "consistency":
    threshold = float(os.environ.get('THRESHOLD_RATIO', 0.2))
else:  # proportionaledges
    threshold = float(os.environ.get('DENSITY', 0.2))

EXCLUDE_PATIENTS = True  # Imposta False per disabilitare il filtro
HARMONIZE_FEATURES = False
# === CONFIG ===
machine="grid"  # "local" | "grid"
if machine=="local":
    base_path = "Y:/project/shared/data/sample_dataset/analysis/harmonization_final/metrics"
    excel_path = "Y:/project/shared/data/sample_dataset/clinical_data/dataset_summary_only_ALL.xlsx"
    exclude_excel_path = "Y:/project/shared/data/sample_dataset/analysis/harmonization_final/data/subjects_harmonization.xlsx"
elif machine=="cluster":
    base_path  = "/home/usero/empenn_group_storage/private/usero/TEST/DATASET/sample_dataset/Multilayer/metrics"
    excel_path = "/workspace/project/data/sample_dataset/multilayer/dataset_summary_only_ALL.xlsx"
    exclude_excel_path = "/workspace/project/data/sample_dataset/multilayer/subjects_harmonization.xlsx"
os.makedirs(base_path, exist_ok=True)
#threshold_strategy = "consistency"      # "consistency" "proportional_edges"
#threshold = 0.2

harmonizations = ["combat",None] #"covbat", 
modes = [ "structural+functional","multilayer", "structural","functional"]
kernel = "rbf"  # oppure "linear"       
n_splits = 10   # folds per la cv esterna (solo per "nested")
use_feature_selection = False
feature_generator = "pls"   # "none" | "pls"
#generators= ["pls", "pca", "svd", "ica", "nca"]
#feature_generator = "nca"    # "none" | "pls" | "pca" | "svd" | "ica"  | "nca"
#pls_n_components_grid = [10, 25, 50, 100]  # usato se feature_generator=="pls"
resampling = False
resampling_str = "_resampling" if resampling else "_noresampling"

# === NEW CV OPTIONS ===
cv_strategy = "nested" # "nested"| "repeated_tvt"
repeats = 10                   # quante ripetizioni train/val/test
val_size = 0.20                 # quota validation nella parte non di test
test_size = 0.10                # quota test sul totale
simple_val_size = 0.2   # quota validation nel "simple"
simple_test_size = 0.2  # quota test nel "simple"
base_random_state = 42         # seed di base; per ogni repeat uso base_random_state + repeat

# === CLASSIFIER ===
#classifier_name = "LOGREG"  # scegli tra "SVM", "RF", "MLP", "LOGREG", "LIN_SVM", "LDA"
if HARMONIZE_FEATURES:
    saving_path = os.path.join(base_path, "CLASSIFICATION_PAIRWISE_FEAT_HARM", classifier_name)
else:
    saving_path = os.path.join(base_path, "CLASSIFICATION_PAIRWISE", classifier_name)
os.makedirs(saving_path, exist_ok=True)
if threshold_strategy == "consistency":
    tag = f"thr{str(threshold)}"   # identico a prima
else:
    tag = f"dens{str(threshold)}"

def _harm_to_str(h):
    return h if (h is not None and h != "") else "noharm"
## Utilities
def norm_label(s: str) -> str:
    # normalize for comparison: string, strip, single spaces
    s = str(s)
    s = s.strip()
    s = re.sub(r"\s+", " ", s)
    return s
def parse_only_pair(s: str):
    # accepts "A_vs_B" or "A vs B"
    parts = re.split(r"_vs_| vs ", s.strip())
    if len(parts) != 2:
        raise ValueError(f"Invalid ONLY_PAIR format: {s!r} (use 'A_vs_B')")
    a, b = map(norm_label, parts)
    # returns the **canonical** pair in alphabetical order
    return tuple(sorted([a, b]))


# ================================
# Utilities
# ================================
def _comp_importance(components_pxk, comp_weights=None):
    """
    components_pxk: (p x k) component weights/loadings on the original features.
    comp_weights:   (k,) pesi dei componenti (varianza spiegata, var(scores), corr canoniche, …).
    Returns (p,) normalized importance.
    """
    if components_pxk is None:
        return None
    p, k = components_pxk.shape
    if comp_weights is None:
        comp_weights = np.ones(k) / max(k, 1)
    comp_weights = np.asarray(comp_weights).reshape(-1)
    comp_weights = comp_weights / (comp_weights.sum() + 1e-12)
    W2 = components_pxk**2                 # (p x k)
    imp = W2 @ comp_weights.reshape(-1, 1) # (p x 1)
    imp = imp.ravel()
    return imp / (imp.sum() + 1e-12)

def _scores_variance(estimator, X, y=None):
    """
    Variance of component scores (as a proxy for "component weight")
    when explained_variance_ratio_ is not available.
    """
    try:
        Z = estimator.transform(X)
    except Exception:
        return None
    if Z is None or Z.ndim != 2:
        return None
    v = np.var(Z, axis=0)
    if np.allclose(v.sum(), 0):
        return None
    return v / v.sum()
class PLSDAFeatureGenerator(BaseEstimator, TransformerMixin):
    """
    sklearn transformer: generates PLS-DA components from X with categorical y.
    - fit: one-hot encode y, then PLSRegression on (X, Y_oh)
    - transform: returns PLS scores (X -> T)
    - get_feature_names_out: ['PLS_comp1', ..., 'PLS_compA']
    """
    def __init__(self, n_components=10, scale=False, copy=True):
        self.n_components = n_components
        self.scale = scale
        self.copy = copy

    def fit(self, X, y):
        # one-hot encode y (strings -> numbers)
        self._ohe_ = OneHotEncoder(sparse_output=False, handle_unknown="ignore")
        Y = self._ohe_.fit_transform(np.asarray(y).reshape(-1, 1))
        # PLS on X, Y_oh
        self.pls_ = PLSRegression(n_components=self.n_components, scale=self.scale, copy=self.copy)
        self.pls_.fit(X, Y)
        return self

    def transform(self, X):
        # Returns the observation scores (components)
        T = self.pls_.transform(X)  # shape: (n_samples, n_components)
        return T

    def get_feature_names_out(self, input_features=None):
        return np.array([f"PLS_comp{i+1}" for i in range(self.n_components)])

    # utility for VIP/diagnostics
    def onehot(self, y):
        return self._ohe_.transform(np.asarray(y).reshape(-1,1))
    def classes_(self):
        return self._ohe_.categories_[0]
class PCAFeatureGenerator(BaseEstimator, TransformerMixin):
    def __init__(self, n_components=25, svd_solver="auto", whiten=False, random_state=42):
        self.n_components = n_components
        self.svd_solver = svd_solver
        self.whiten = whiten
        self.random_state = random_state
    def fit(self, X, y=None):
        self.pca_ = PCA(n_components=self.n_components, svd_solver=self.svd_solver,
                        whiten=self.whiten, random_state=self.random_state)
        self.pca_.fit(X)
        # components_: (k x p) -> transpose -> (p x k)
        self.components_ = self.pca_.components_.T
        self.comp_weights_ = getattr(self.pca_, "explained_variance_ratio_", None)
        return self
    def transform(self, X):
        return self.pca_.transform(X)
    def get_feature_names_out(self, input_features=None):
        k = self.components_.shape[1]
        return np.array([f"PCA_comp{i+1}" for i in range(k)])

class TruncatedSVDFeatureGenerator(BaseEstimator, TransformerMixin):
    """Recommended for sparse matrices."""
    def __init__(self, n_components=50, random_state=42):
        self.n_components = n_components
        self.random_state = random_state
    def fit(self, X, y=None):
        self.svd_ = TruncatedSVD(n_components=self.n_components, random_state=self.random_state)
        self.svd_.fit(X)
        self.components_ = self.svd_.components_.T
        self.comp_weights_ = getattr(self.svd_, "explained_variance_ratio_", None)
        return self
    def transform(self, X):
        return self.svd_.transform(X)
    def get_feature_names_out(self, input_features=None):
        return np.array([f"SVD_comp{i+1}" for i in range(self.n_components)])

class NMFFeatureGenerator(BaseEstimator, TransformerMixin):
    def __init__(self, n_components=25, init="nndsvda", max_iter=400, random_state=42):
        self.n_components = n_components
        self.init = init
        self.max_iter = max_iter
        self.random_state = random_state
    def fit(self, X, y=None):
        self.nmf_ = NMF(n_components=self.n_components, init=self.init,
                        max_iter=self.max_iter, random_state=self.random_state)
        W = self.nmf_.fit_transform(X)     # (n x k)
        H = self.nmf_.components_          # (k x p)
        self.components_ = H.T             # (p x k)
        # pesi dai punteggi:
        self.comp_weights_ = _scores_variance(self.nmf_, X)
        return self
    def transform(self, X):
        return self.nmf_.transform(X)
    def get_feature_names_out(self, input_features=None):
        return np.array([f"NMF_comp{i+1}" for i in range(self.n_components)])

class ICAFeatureGenerator(BaseEstimator, TransformerMixin):
    def __init__(self, n_components=25, max_iter=1000, tol=1e-3, random_state=42):
        self.n_components = n_components
        self.max_iter = max_iter
        self.tol = tol
        self.random_state = random_state
    def fit(self, X, y=None):
        self.ica_ = FastICA(n_components=self.n_components, max_iter=self.max_iter,
                            tol=self.tol, random_state=self.random_state)
        Z = self.ica_.fit_transform(X)         # (n x k)
        # In ICA, mixing_ is the inverse of components_; use mixing_ (p x k)
        self.components_ = getattr(self.ica_, "mixing_", None)
        if self.components_ is None:
            # fallback: use components_.T if mixing_ is not available
            self.components_ = self.ica_.components_.T
        self.comp_weights_ = _scores_variance(self.ica_, X)
        return self
    def transform(self, X):
        return self.ica_.transform(X)
    def get_feature_names_out(self, input_features=None):
        return np.array([f"ICA_comp{i+1}" for i in range(self.n_components)])

class CCAFeatureGenerator(BaseEstimator, TransformerMixin):
    def __init__(self, n_components=2, scale=False, max_iter=5000):
        self.n_components = n_components
        self.scale = scale
        self.max_iter = max_iter
    def fit(self, X, y):
        self._ohe_ = OneHotEncoder(sparse_output=False, handle_unknown="ignore")
        Y = self._ohe_.fit_transform(np.asarray(y).reshape(-1, 1))
        self.cca_ = CCA(n_components=self.n_components, scale=self.scale, max_iter=self.max_iter)
        self.cca_.fit(X, Y)
        self.components_ = self.cca_.x_weights_         # (p x k)
        U, V = self.cca_.transform(X, Y)
        cors = []
        for i in range(U.shape[1]):
            c = np.corrcoef(U[:, i], V[:, i])[0, 1]
            cors.append(max(0.0, float(c)))
        cors = np.array(cors)
        self.comp_weights_ = cors / (cors.sum() + 1e-12)
        return self
    def transform(self, X):
        # usa Y dummy solo per shape; sklearn ignora Y in transform
        Yd = np.zeros((X.shape[0], len(self._ohe_.categories_[0])))
        U, _ = self.cca_.transform(X, Yd)
        return U
    def get_feature_names_out(self, input_features=None):
        return np.array([f"CCA_comp{i+1}" for i in range(self.n_components)])

class NCAFeatureGenerator(BaseEstimator, TransformerMixin):
    def __init__(self, n_components=10, max_iter=1000, random_state=42):
        self.n_components = n_components
        self.max_iter = max_iter
        self.random_state = random_state
    def fit(self, X, y):
        self.nca_ = NeighborhoodComponentsAnalysis(n_components=self.n_components,
                                                   max_iter=self.max_iter,
                                                   random_state=self.random_state)
        self.nca_.fit(X, y)
        # sklearn exposes components_ (k x p)
        self.components_ = self.nca_.components_.T   # (p x k)
        self.comp_weights_ = _scores_variance(self.nca_, X, y)
        return self
    def transform(self, X):
        return self.nca_.transform(X)
    def get_feature_names_out(self, input_features=None):
        return np.array([f"NCA_comp{i+1}" for i in range(self.n_components)])
def _compute_vip(pls_or_wrapper, X, Y):
    """
    VIP per PLS-DA multi-classe (PLS2), con normalizzazione per ||w_a||^2.
    Ritorna: vip (p,)
    """
    import numpy as np
    pls = getattr(pls_or_wrapper, "pls_", pls_or_wrapper)

    T = pls.x_scores_      # (n, A)
    W = pls.x_weights_     # (p, A)
    Q = pls.y_loadings_    # (c, A)

    p, A = W.shape
    if A == 0:
        return np.zeros(p)

    # Quota di varianza di Y spiegata dal comp. a (proxy SSY)
    SSY_t = np.sum(T**2, axis=0)        # (A,)
    SSY_q = np.sum(Q**2, axis=0)        # (A,)
    SSY   = SSY_t * SSY_q               # (A,)

    # Normalizzazione per la lunghezza del vettore pesi del componente
    w_norm2 = np.sum(W**2, axis=0)      # (A,)
    w_norm2[w_norm2 == 0] = 1.0

    num = ( (W**2) / w_norm2[None, :] ) @ SSY.reshape(-1, 1)   # (p,1)
    den = SSY.sum() if SSY.sum() != 0 else 1.0

    vip = np.sqrt( (p * num.ravel()) / den )
    return vip
def make_k_grid(n_features_available, base=[50, 100, 200, 500, 1000, 1500, 2000]):
    k_list = sorted({k for k in base if k <= n_features_available})
    if n_features_available > 200:
        k_list += [int(0.05*n_features_available), int(0.1*n_features_available)]
    k_list = sorted({k for k in k_list if 2 <= k <= n_features_available})
    if not k_list:
        k_list = [min(50, n_features_available)]
    return k_list

def cv_suffix(cv_strategy, n_splits, repeats):
    if cv_strategy == "nested":
        return f"_cv{n_splits}_nested"
    elif cv_strategy == "simple":
        return "_cv1_simple"
    elif cv_strategy == "repeated_tvt":
        return f"_cv{repeats}_repeatedtvt"
    else:
        return "_cvUNK"

def cv_num_runs(cv_strategy, n_splits, repeats):
    if cv_strategy == "nested":
        return n_splits
    elif cv_strategy == "simple":
        return 1
    elif cv_strategy == "repeated_tvt":
        return repeats
    else:
        return 0
def _comp_importance(components_pxk, comp_weights=None):
    if components_pxk is None:
        return None
    p, k = components_pxk.shape
    if k == 0:
        return None
    if comp_weights is None:
        comp_weights = np.ones(k) / k
    comp_weights = np.asarray(comp_weights).reshape(-1)
    s = comp_weights.sum()
    if s <= 0:
        comp_weights = np.ones(k) / k
    else:
        comp_weights = comp_weights / s
    W2 = components_pxk ** 2            # (p x k)
    imp = W2 @ comp_weights[:, None]    # (p x 1)
    imp = imp.ravel()
    s = imp.sum()
    if s > 0:
        imp = imp / s
    return imp

def _scores_variance(estimator, X, y=None):
    try:
        Z = estimator.transform(X)
    except Exception:
        return None
    if Z is None or Z.ndim != 2 or Z.shape[1] == 0:
        return None
    v = np.var(Z, axis=0)
    if not np.isfinite(v).any():
        return None
    s = v.sum()
    if s <= 0:
        return None
    return v / s

def top_k_from_importances(imp, n_components_hint=None, p=None):
    """Chooses a reasonable K for reporting selected_features."""
    if imp is None or imp.size == 0:
        return []
    if n_components_hint is None:
        # fallback: 25 or 5% of features, capped at 100
        if p is None: p = imp.size
        k = max(1, min(100, max(25, int(0.05 * p))))
    else:
        k = int(max(1, n_components_hint))
    top_idx = np.argsort(imp)[::-1][:k]
    return top_idx.tolist()

def _current_fold_index():
    # try to use the variables if they exist; fallback to 0
    return (fold_idx if 'fold_idx' in globals() or 'fold_idx' in locals() else
            (repeat if 'repeat' in globals() or 'repeat' in locals() else 0))
def log_feature_importances_for_generators(best_estimator, X_train, y_train,
                                           mode, harm_str, label_a, label_b):
    global feature_importance_logs, selected_features_logs

    if "feature_importance_logs" not in globals():
        feature_importance_logs = []
    if "selected_features_logs" not in globals():
        selected_features_logs = []

    pair_name = f"{label_a}_vs_{label_b}"
    fold_id = _current_fold_index()
    features = X_train.columns

    # --- 1) PLS with VIP ---
    pls_step = best_estimator.named_steps.get("plsdafeaturegenerator")
    if pls_step is not None:
        X_for_vip = X_train  # evitiamo locals().get; X_gs spesso è X_train standardizzato internamente
        y_for_vip = y_train

        # one-hot encoding consistent with the wrapper
        try:
            Y_oh = pls_step.onehot(y_for_vip)
            vip = _compute_vip(pls_step, X_for_vip.values, Y_oh)
            vip = np.asarray(vip).ravel()
            if vip.size == len(features):
                vip_df_local = pd.DataFrame({
                    "mode": mode, "harmonization": harm_str, "pair": pair_name,
                    "fold": fold_id,
                    "metric": "VIP", "feature": features, "importance_mean": vip, "importance_std": np.nan
                })
                feature_importance_logs.append(vip_df_local)

                # top-K coerente con n_components
                k_selected = int(getattr(pls_step, "n_components", 10))
                top_idx = top_k_from_importances(vip, n_components_hint=k_selected, p=len(features))
                # header + righe
                selected_features_logs.append({
                    "mode": mode, "harmonization": harm_str, "pair": pair_name,
                    "fold": fold_id, "k_selected": len(top_idx), "feature": None
                })
                for feat in features[top_idx]:
                    selected_features_logs.append({
                        "mode": mode, "harmonization": harm_str, "pair": pair_name,
                        "fold": fold_id, "k_selected": len(top_idx), "feature": feat
                    })
        except Exception as e:
            print(f"[WARN] PLS VIP not computed: {e}")

    # --- 2) Other interpretable generators ---
    step_specs = [
        ("pcafeaturegenerator",          "PCA"),
        ("truncatedsvdfeaturegenerator", "SVD"),
        ("nmffeaturegenerator",          "NMF"),
        ("icafeaturegenerator",          "ICA"),
        ("ccafeaturegenerator",          "CCA"),
        ("ncafeaturegenerator",          "NCA"),
    ]
    for step_name, tag in step_specs:
        step = best_estimator.named_steps.get(step_name)
        if step is None:
            continue

        comps = getattr(step, "components_", None)     # (p x k) atteso
        if comps is None:
            # ICA sometimes exposes mixing_ (p x k)
            comps = getattr(step, "mixing_", None)
        if comps is None:
            continue

        # component weights: prefer explicit attribute, otherwise use var(scores)
        compw = getattr(step, "comp_weights_", None)
        if compw is None:
            compw = _scores_variance(step, X_train, y_train)

        imp = _comp_importance(comps, comp_weights=compw)
        if imp is None or imp.size != len(features):
            # fallback to absolute weight magnitude only
            try:
                imp = np.abs(comps).sum(axis=1)
                s = imp.sum()
                if s > 0: imp = imp / s
            except Exception:
                continue

        # log importanze
        fi_df = pd.DataFrame({
            "mode": mode, "harmonization": harm_str, "pair": pair_name,
            "fold": fold_id,
            "metric": f"{tag}_IMP",
            "feature": features,
            "importance_mean": imp,
            "importance_std": np.nan
        })
        feature_importance_logs.append(fi_df)

        # reporting-oriented selection: k ~ n_components
        n_comp_hint = getattr(step, "n_components", None)
        top_idx = top_k_from_importances(imp, n_components_hint=n_comp_hint, p=len(features))

        selected_features_logs.append({
            "mode": mode, "harmonization": harm_str, "pair": pair_name,
            "fold": fold_id, "k_selected": len(top_idx), "feature": None
        })
        for feat in features[top_idx]:
            selected_features_logs.append({
                "mode": mode, "harmonization": harm_str, "pair": pair_name,
                "fold": fold_id, "k_selected": len(top_idx), "feature": feat
            })
CV_TAG = cv_suffix(cv_strategy, n_splits, repeats if 'repeats' in globals() else 1)
NUM_RUNS = cv_num_runs(cv_strategy, n_splits, repeats if 'repeats' in globals() else 1)

# ================================
# Pipelines and hyperparameter grids
# ================================
def build_pipe_and_param_grid(classifier_name, kernel, X_train, use_feature_selection=True):
    nF = X_train.shape[1]
    k_grid = make_k_grid(nF)

    # === Step opzionali: feature generator + selettore ===
    # Generator: o "passthrough" o PLSRegression (componenti supervisionate)
    if feature_generator == "pls":
        generator = PLSDAFeatureGenerator(n_components=10, scale=False)
    elif feature_generator == "pca":
        generator = PCAFeatureGenerator(n_components=25)
    elif feature_generator == "svd":
        generator = TruncatedSVDFeatureGenerator(n_components=50)
    elif feature_generator == "nmf":
        generator = NMFFeatureGenerator(n_components=25)
    elif feature_generator == "ica":
        generator = ICAFeatureGenerator(n_components=25)
    elif feature_generator == "cca":
        generator = CCAFeatureGenerator(n_components=3, scale=False)
    elif feature_generator == "nca":
        generator = NCAFeatureGenerator(n_components=10)
    elif feature_generator == "none":
        generator = "passthrough"
    else:
        raise ValueError("feature_generator must be in {'none','pls','pca','svd','nmf','ica','cca','nca'}")

    # Selettore su input (o sulle componenti, se vuoi): qui lo manteniamo SOLO quando feature_generator=='none'
    if use_feature_selection and feature_generator == "none":
        selector = SelectKBest(f_classif)
        selector_grid = {"selectkbest__k": k_grid}
    else:
        selector = "passthrough"
        selector_grid = {}
    name = classifier_name.upper()
    # === Classifier pipes ===
    if name == "SVM":
        pipe = make_pipeline(
            StandardScaler(),
            generator,
            selector,
            SVC(kernel=kernel, class_weight='balanced', probability=True, random_state=42)
        )
        if kernel == "rbf":
            clf_grid = {
                "svc__C": np.logspace(-3, 3, 7),
                "svc__gamma": np.logspace(-4, 0, 5),
            }
        else:
            raise ValueError("Kernel SVM non gestito in questa griglia.")

    elif name == "RF":
        pipe = make_pipeline(
            StandardScaler(),   # innocuo
            generator,
            selector,
            RandomForestClassifier(
                n_estimators=300, class_weight='balanced',
                random_state=42, n_jobs=1  # evita parallelismo annidato
            )
        )
        clf_grid = {
            "randomforestclassifier__max_depth": [None, 10, 20, 40],
            "randomforestclassifier__min_samples_leaf": [1, 2, 4],
            "randomforestclassifier__max_features": ["sqrt", 0.1, 0.25, 0.5],
        }

    elif name == "MLP":
        pipe = make_pipeline(
            StandardScaler(),
            generator,
            selector,
            MLPClassifier(
                max_iter=1000,
                random_state=42,
                early_stopping=True,
                learning_rate="constant",
                solver="lbfgs"
            )
        )
        clf_grid = {
            "mlpclassifier__hidden_layer_sizes": [
                (50),(25), (25,10), (25, 10, 5), (50,25), (50,25,10), (50,25,10,5)
            ],
            "mlpclassifier__activation": ["relu", "tanh"],
            "mlpclassifier__alpha": [1e-4, 1e-3, 1e-2, 1e-1],
            "mlpclassifier__max_iter": [100, 500, 1000],
        }
    elif name in ("LOGREG", "LR"):
        # For pairwise 2-class, liblinear works well
        pipe = make_pipeline(
            StandardScaler(),
            generator,
            selector,
            LogisticRegression(
                max_iter=1000,
                class_weight="balanced",
                solver="liblinear"  # supporta L1 e L2
            )
        )
        clf_grid = {
            "logisticregression__C":      np.logspace(-3, 3, 7),
            "logisticregression__penalty": ["l1", "l2"]
        }

    # === NEW: Separate Linear SVM (LinearSVC) ===
    elif name in ("LIN_SVM", "LSVM", "LINEAR_SVM"):
        # LinearSVC does not have probability=True, but only predict is used here
        pipe = make_pipeline(
            StandardScaler(),
            generator,
            selector,
            LinearSVC(
                class_weight="balanced",
                dual="auto",          
                max_iter=10000,
                random_state=42
            )
        )
        clf_grid = {
            "linearsvc__C": np.logspace(-3, 3, 7)
        }

    # === NEW: Linear Discriminant Analysis (LDA) ===
    elif name == "LDA":
        # LDA does not support class_weight; class imbalance can be handled with resampling
        pipe = make_pipeline(
            StandardScaler(),   # spesso aiuta se le feature hanno scale diverse
            generator,
            selector,
            LinearDiscriminantAnalysis()
        )
        # Minimal reasonable grid: solver + shrinkage
        # shrinkage è usato solo se solver in {'lsqr','eigen'}
        clf_grid = {
            "lineardiscriminantanalysis__solver":   ["svd", "lsqr"],
            "lineardiscriminantanalysis__shrinkage": [None, "auto"]  # used only if solver != 'svd'
        }

    else:
        raise ValueError(
            "classifier_name must be one of {'SVM','RF','MLP','LOGREG','LIN_SVM','LDA'}"
        )


    # === Full parameter grid ===
    param_grid = {**selector_grid, **clf_grid}
    if feature_generator == "pls":
        param_grid.update({"plsdafeaturegenerator__n_components": [5,10, 25]})
    elif feature_generator == "pca":
        param_grid.update({"pcafeaturegenerator__n_components": [10, 25, 50]})
    elif feature_generator == "svd":
        param_grid.update({"truncatedsvdfeaturegenerator__n_components": [10,25,50]})
    elif feature_generator == "nmf":
        param_grid.update({"nmffeaturegenerator__n_components": [10, 25, 50]})
    elif feature_generator == "ica":
        param_grid.update({"icafeaturegenerator__n_components": [10, 25, 50]})
    elif feature_generator == "cca":
        param_grid.update({"ccafeaturegenerator__n_components": [5,10,25]})
    elif feature_generator == "nca":
        param_grid.update({"ncafeaturegenerator__n_components": [5, 10, 25]})

    return pipe, param_grid
def harmonize_features(
    df: pd.DataFrame,
    train_ids,
    *,
    id_col: str = "id",
    label_col: str = "label",
    batch_col: str = "batch",          # <- MUST exist in df (e.g. site/scanner)
    covar_cols=["Age","Sex"],                  # e.g. ["age", "sex"] if you want to preserve them
    smooth_terms=None,                # e.g. ["age"] for ComBat-GAM (if using GAM)
    eb: bool = True,
    ref_batch=None,
    return_model: bool = False,
):
    from neuroHarmonize import harmonizationLearn, harmonizationApply  # API NeuroHarmonize [web:2][web:6]

    if covar_cols is None:
        covar_cols = []

    df_h = df.copy()

    if id_col not in df_h.columns:
        raise ValueError(f"Missing column '{id_col}' in df.")
    if batch_col not in df_h.columns:
        raise ValueError(
            f"Missing batch column '{batch_col}'. "
            f"A site/scanner variable is required for ComBat/NeuroHarmonize."
        )

    # Features = everything except id/label/batch/covariates
    exclude = {id_col, label_col, batch_col, *covar_cols}
    feature_cols = [c for c in df_h.columns if c not in exclude]

    # Covariate matrix (batch + covariates) aligned with df
    covars = df_h[[batch_col] + list(covar_cols)].copy()
    covars = covars.rename(columns={batch_col: "SITE"})
    # Encode categorical variables -> numeric codes (NeuroHarmonize requires numeric covariates)
    for c in covars.columns:
        if not pd.api.types.is_numeric_dtype(covars[c]):
            covars[c] = pd.Categorical(covars[c]).codes

    # Check NaN values explicitly before harmonization
    if df_h[feature_cols].isna().any().any():
        raise ValueError("NaN values found in features: handle them before harmonization.")
    if covars.isna().any().any():
        raise ValueError("NaN values found in covariates: handle them before harmonization.")

    train_ids = set(train_ids)
    train_mask = df_h[id_col].isin(train_ids)


    if train_mask.sum() == 0:
        raise ValueError("train_ids does not match any subject in df.")

    X_train = df_h.loc[train_mask, feature_cols].to_numpy(dtype=float)
    cov_train = covars.loc[train_mask].reset_index(drop=True)
    var_train = np.var(X_train, axis=0)
    ok = np.isfinite(var_train) & (var_train > 0)

    feat_ok = [c for c, keep in zip(feature_cols, ok) if keep]
    feat_bad = [c for c, keep in zip(feature_cols, ok) if not keep]

    # harmonize only feat_ok
    Xtr_ok = df_h.loc[train_mask, feat_ok].to_numpy(float)
    Xall_ok = df_h.loc[:, feat_ok].to_numpy(float)

    model, _ = harmonizationLearn(Xtr_ok, cov_train)
    Xall_ok_adj = harmonizationApply(Xall_ok, covars, model)

    df_h.loc[:, feat_ok] = Xall_ok_adj
    
    df_h = df_h.drop(columns=feat_bad)
    if return_model:
        return df_h, model
    else:
        return df_h

# === Load labels ===
info = pd.read_excel(excel_path, sheet_name="ALL")[["Subject", "Diagnosis_amyloid","Batch_variable_simplified","Age","Sex"]]
info = info.rename(columns={"Subject": "id", "Diagnosis_amyloid": "label", "Batch_variable_simplified": "batch"})
if EXCLUDE_PATIENTS:
    exclude_df = pd.read_excel(exclude_excel_path, sheet_name="Sheet1")
    excluded_ids = set(exclude_df['Subject'].dropna().unique())
    print(f"\n{'='*60}")
    print(f"FILTER ACTIVE: {len(excluded_ids)} subjects will be excluded")
    print(f"Examples of excluded IDs: {list(excluded_ids)[:5]}")
    print(f"{'='*60}\n")
else:
    excluded_ids = set()
    print("\nFILTER DISABLED: tutti i pazienti saranno inclusi\n")

results = []
all_confusion_matrices = {}
selected_features_logs = []
feature_importance_logs = []
param_logs = []
for mode in modes:
    for harm in harmonizations:
        if threshold_strategy == "consistency":
            thr_str = f"thr{str(threshold)}"   # identico a prima
        else:
            thr_str = f"dens{str(threshold)}"
        print(f"Processing mode={mode}, harmonization={harm}, harmonization_features={HARMONIZE_FEATURES} thr={thr_str}, classifier={classifier_name}, use_feature_selection={use_feature_selection}, feature_generator={feature_generator}, cv={cv_strategy}")
        harm_str = harm if harm else "noharm"
        if threshold_strategy == "consistency":
            thr_str = f"thr{str(threshold)}"   # identico a prima
        else:
            thr_str = f"dens{str(threshold)}"

        # === Caricamento dati ===
        if mode == "structural+functional":
            if HARMONIZE_FEATURES:
                struct_nodal_file = f"nodal_metrics_structural_noharm_{thr_str}.csv"
                #struct_global_file = f"global_metrics_structural_noharm_{thr_str}.csv"
                func_nodal_file   = f"nodal_metrics_functional_noharm_{thr_str}.csv"
                #func_global_file  = f"global_metrics_functional_noharm_{thr_str}.csv"
            else:
                struct_nodal_file = f"nodal_metrics_structural_{harm_str}_{thr_str}.csv"
                #struct_global_file = f"global_metrics_structural_{harm_str}_{thr_str}.csv"
                func_nodal_file   = f"nodal_metrics_functional_{harm_str}_{thr_str}.csv"
                #func_global_file  = f"global_metrics_functional_{harm_str}_{thr_str}.csv"
            struct_nodal_path = os.path.join(base_path, struct_nodal_file)
            #struct_global_path= os.path.join(base_path, struct_global_file)
            func_nodal_path   = os.path.join(base_path, func_nodal_file)
            #func_global_path  = os.path.join(base_path, func_global_file)
            #if not all(map(os.path.exists, [struct_nodal_path, struct_global_path, func_nodal_path, func_global_path])):
            if not all(map(os.path.exists, [struct_nodal_path, func_nodal_path])):
                print(f"Skipping {mode} - {harm_str}: one or more files missing.")
                continue
            struct_nodal_df = pd.read_csv(struct_nodal_path)
            #struct_global_df= pd.read_csv(struct_global_path)
            func_nodal_df   = pd.read_csv(func_nodal_path)
            #func_global_df  = pd.read_csv(func_global_path)
            #struct_data = pd.merge(struct_nodal_df, struct_global_df, on="id")
            #func_data   = pd.merge(func_nodal_df, func_global_df, on="id")
            struct_data = struct_nodal_df
            func_data   = func_nodal_df
            struct_features = struct_data.drop(columns=["id"]).add_prefix("struct_")
            func_features   = func_data.drop(columns=["id"]).add_prefix("func_")
            combined_data = (struct_data[["id"]].join(struct_features).merge(func_data[["id"]].join(func_features), on="id", how="inner"))
            data = pd.merge(combined_data, info, on="id").dropna()
        else:
            if HARMONIZE_FEATURES:
                nodal_file  = f"nodal_metrics_{mode}_noharm_{thr_str}.csv"
                #global_file = f"global_metrics_{mode}_noharm_{thr_str}.csv"
            else:
                nodal_file  = f"nodal_metrics_{mode}_{harm_str}_{thr_str}.csv"
                #global_file = f"global_metrics_{mode}_{harm_str}_{thr_str}.csv"
            nodal_path  = os.path.join(base_path, nodal_file)
            #global_path = os.path.join(base_path, global_file)
            #   if not (os.path.exists(nodal_path) and os.path.exists(global_path)):
            if not os.path.exists(nodal_path):
                print(f"Skipping {mode} - {harm_str}: file missing.")
                continue
            nodal_df  = pd.read_csv(nodal_path)
            #global_df = pd.read_csv(global_path)
            #data = pd.merge(nodal_df, global_df, on="id")
            data = nodal_df.copy()
            data = pd.merge(data, info, on="id").dropna()
        if HARMONIZE_FEATURES and harm_str=="combat":
            data= harmonize_features(data, excluded_ids)

        
        if EXCLUDE_PATIENTS and len(excluded_ids) > 0:
            n_before = len(data)
            data = data[~data['id'].isin(excluded_ids)]
            n_after = len(data)
            print(f"  Subject filter: {n_before} -> {n_after} soggetti ({n_before - n_after} esclusi)")
        data = data.sort_values('id').reset_index(drop=True)
        X = data.drop(columns=["id", "label", "batch", "Age", "Sex"])
        y = data["label"]
        pairs = list(combinations(sorted(y.unique()), 2))

        for (label_a, label_b) in pairs:
            print(f"\n=== Pairwise {label_a} vs {label_b} | mode={mode}, harmonization={harm_str} ===")
            pair_key = f"{mode}_{harm_str}_{label_a}_vs_{label_b}"
            mask = y.isin([label_a, label_b])
            X_pair = X[mask].copy()
            y_pair = y[mask].copy()

            # initialize accumulator
            all_confusion_matrices[pair_key] = None

            if cv_strategy == "nested":
                # === OUTER CV ===
                sss = StratifiedKFold(n_splits=n_splits, shuffle=True, random_state=42)
                inner_cv = StratifiedShuffleSplit(n_splits=10, test_size=0.2, random_state=42)

                for fold_idx, (train_idx, test_idx) in enumerate(tqdm(sss.split(X_pair, y_pair), total=n_splits)):
                    X_train, X_test = X_pair.iloc[train_idx], X_pair.iloc[test_idx]
                    y_train, y_test = y_pair.iloc[train_idx], y_pair.iloc[test_idx]

                    # Resampling only on the training set
                    if resampling:
                        cc = CondensedNearestNeighbour(random_state=42)
                        X_train, y_train = cc.fit_resample(X_train, y_train)

                    pipe, param_grid = build_pipe_and_param_grid(classifier_name, kernel, X_train, use_feature_selection)
                    gs = GridSearchCV(
                        estimator=pipe,
                        param_grid=param_grid,
                        cv=inner_cv,
                        scoring="balanced_accuracy",
                        n_jobs=-1,   # evita problemi di pickling su alcuni ambienti
                        refit=True,
                        verbose=0
                    )
                    #with parallel_backend("threading"):
                    gs.fit(X_train, y_train)
                    y_pred = gs.predict(X_test)
                    # --- LOG selected parameters ---
                    best_params = gs.best_params_.copy()
                    best_val_score = getattr(gs, "best_score_", np.nan)

                    # compact console output
                    print(f"[{pair_key}] fold={fold_idx if 'fold_idx' in locals() else (repeat if 'repeat' in locals() else 0)} "
                        f"| best_val={best_val_score:.3f} | params={best_params}")
                    
                    # log row: meta + json + flattened parameters (Pandas will insert NaN where missing)
                    meta = {
                        "mode": mode,
                        "harmonization": harm_str,
                        "pair": f"{label_a}_vs_{label_b}",
                        "fold": fold_idx if 'fold_idx' in locals() else (repeat if 'repeat' in locals() else 0),
                        "cv_strategy": cv_strategy,
                        "classifier": classifier_name,
                        "best_val_score": best_val_score,
                        "best_params_json": json.dumps(best_params)
                    }
                    flatten = {k: v for k, v in best_params.items()}  # es. selectkbest__k, svc__C, mlpclassifier__alpha, ...
                    param_logs.append({**meta, **flatten})
                    best = gs.best_estimator_
                    log_feature_importances_for_generators(best, X_train, y_train, mode, harm_str, label_a, label_b)

                    # metrics
                    bal_acc_test = balanced_accuracy_score(y_test, y_pred)
                    f1_test      = f1_score(y_test, y_pred, average="weighted")
                    kappa_test   = cohen_kappa_score(y_test, y_pred)
                    accuracy_test= accuracy_score(y_test, y_pred)

                    cm_fold = confusion_matrix(y_test, y_pred, labels=[label_a, label_b])
                    all_confusion_matrices[pair_key] = cm_fold if all_confusion_matrices[pair_key] is None else (all_confusion_matrices[pair_key] + cm_fold)

                    results.append({
                       "mode": mode, "harmonization": harm_str, "threshold": threshold,
                        "fold": fold_idx, "type": "test",
                        "cv_strategy": cv_strategy,
                        "accuracy": accuracy_test, "balanced_accuracy": bal_acc_test,
                        "f1_score": f1_test, "cohen_kappa": kappa_test,
                        "classifier": classifier_name, "pair": f"{label_a}_vs_{label_b}"
                    }) 

            elif cv_strategy == "repeated_tvt":
                # === REPEATED TRAIN/VAL/TEST  ===
                pbar = tqdm(range(repeats), total=repeats,
                desc=f"{pair_key} repeats", leave=False)
                for repeat in range(repeats):
                    seed = base_random_state + repeat

                    # 1) Split TEST once (stratified)
                    X_tmp, X_test, y_tmp, y_test = train_test_split(
                        X_pair, y_pair,
                        test_size=test_size, stratify=y_pair, random_state=seed
                    )
                    # 2) Split VALIDATION from the remaining data (stratified)
                    val_ratio_of_tmp = val_size / (1.0 - test_size)
                    X_train, X_val, y_train, y_val = train_test_split(
                        X_tmp, y_tmp,
                        test_size=val_ratio_of_tmp, stratify=y_tmp, random_state=seed
                    )

                    # 3) (optional) resampling only on the training set
                    if resampling:
                        cc = CondensedNearestNeighbour(random_state=seed)
                        X_train, y_train = cc.fit_resample(X_train, y_train)

                    # 4) GridSearch using ONLY validation via PredefinedSplit
                    #    Costruiamo X_gs = [X_train; X_val] con maschera test_fold: -1=train, 0=val
                    X_gs = pd.concat([X_train, X_val], axis=0)
                    y_gs = pd.concat([y_train, y_val], axis=0).reset_index(drop=True)
                    test_fold = np.array([-1]*len(X_train) + [0]*len(X_val))
                    ps = PredefinedSplit(test_fold)

                    pipe, param_grid = build_pipe_and_param_grid(classifier_name, kernel, X_train, use_feature_selection)
                    gs = GridSearchCV(
                        estimator=pipe,
                        param_grid=param_grid,
                        cv=ps,                        # valuta sulla validation
                        scoring="balanced_accuracy",
                        n_jobs=8,                     # evita parallelismo annidato
                        refit=True,                   # refit su train+val
                        verbose=0
                    )
                    #with parallel_backend("threading"):
                    gs.fit(X_gs, y_gs)

                    # 5) Final test on the test split
                    y_pred = gs.predict(X_test)
                    # --- LOG selected parameters ---
                    best_params = gs.best_params_.copy()
                    best_val_score = getattr(gs, "best_score_", np.nan)

                    # compact console output
                    print(f"[{pair_key}] fold={fold_idx if 'fold_idx' in locals() else (repeat if 'repeat' in locals() else 0)} "
                        f"| best_val={best_val_score:.3f} | params={best_params}")
                    
                    # log row: meta + json + flattened parameters (Pandas will insert NaN where missing)
                    meta = {
                        "mode": mode,
                        "harmonization": harm_str,
                        "pair": f"{label_a}_vs_{label_b}",
                        "fold": fold_idx if 'fold_idx' in locals() else (repeat if 'repeat' in locals() else 0),
                        "cv_strategy": cv_strategy,
                        "classifier": classifier_name,
                        "best_val_score": best_val_score,
                        "best_params_json": json.dumps(best_params)
                    }
                    flatten = {k: v for k, v in best_params.items()}  # es. selectkbest__k, svc__C, mlpclassifier__alpha, ...
                    param_logs.append({**meta, **flatten})
                    # Log selected features (fold = repeat)
                    best = gs.best_estimator_
                    log_feature_importances_for_generators(best, X_train, y_train, mode, harm_str, label_a, label_b)
                    # metrics (fold=repeat)
                    bal_acc_test = balanced_accuracy_score(y_test, y_pred)
                    f1_test      = f1_score(y_test, y_pred, average="weighted")
                    kappa_test   = cohen_kappa_score(y_test, y_pred)
                    accuracy_test= accuracy_score(y_test, y_pred)

                    cm_fold = confusion_matrix(y_test, y_pred, labels=[label_a, label_b])
                    all_confusion_matrices[pair_key] = cm_fold if all_confusion_matrices[pair_key] is None else (all_confusion_matrices[pair_key] + cm_fold)

                    results.append({
                        "mode": mode, "harmonization": harm_str, "threshold": threshold,
                        "fold": repeat, "type": "test",
                        "cv_strategy": cv_strategy,
                        "accuracy": accuracy_test, "balanced_accuracy": bal_acc_test,
                        "f1_score": f1_test, "cohen_kappa": kappa_test,
                        "classifier": classifier_name, "pair": f"{label_a}_vs_{label_b}"
                    })
                    # NEW: update progress bar with a quick summary
                    try:
                        pbar.set_postfix({
                            "rep": repeat+1,
                            "BA_test": f"{balanced_accuracy_score(y_test, y_pred):.3f}",
                            "BA_val":  f"{gs.best_score_:.3f}"  # score CV sulla validation
                        })
                    except Exception:
                        pass
                    pbar.update(1)
                pbar.close()
            else:
                raise ValueError("cv_strategy must be 'nested' or 'simple'.")

# === SAVING OUTPUTS ===
results_df = pd.DataFrame(results)
base_fs_str = "_fs" if use_feature_selection else "_nofs"

# 2) add the Feature Generator (static)
if feature_generator == "none":
    fg_part = ""
else:
    fg_part = f"_FG-{feature_generator}"

# 3) initial fs_str used BEFORE fitting (without FG hyperparameters)
fs_str = f"{base_fs_str}{fg_part}"
if threshold_strategy == "consistency":
    thr_str = f"thr{str(threshold)}"   # identico a prima
else:
    thr_str = f"dens{str(threshold)}"

# metrics (PAIRWISE)
results_csv = os.path.join(
    saving_path,
    f"{classifier_name}_PAIRWISE_metrics_{thr_str}{fs_str}{resampling_str}{CV_TAG}_pairwise.csv"
)
results_df.to_csv(results_csv, index=False)

# Confusion matrix plots
def plot_confusion_matrix(cm, labels, title, save_path):
    plt.figure(figsize=(8, 6))
    with np.errstate(divide='ignore', invalid='ignore'):
        cm_normalized = cm.astype('float') / cm.sum(axis=1, keepdims=True) * 100
        cm_normalized = np.nan_to_num(cm_normalized, nan=0.0, posinf=0.0, neginf=0.0)
    ax = sns.heatmap(cm, annot=True, fmt='d', cmap='Blues',
                    xticklabels=labels, yticklabels=labels,
                    cbar_kws={'label': 'Count'})
    for i in range(len(labels)):
        for j in range(len(labels)):
            ax.text(j + 0.5, i + 0.7, f'({cm_normalized[i, j]:.1f}%)',
                    ha='center', va='center', fontsize=9, color='red')
    plt.title(title)
    plt.ylabel('True Label')
    plt.xlabel('Predicted Label')
    plt.tight_layout()
    plt.savefig(save_path, dpi=300, bbox_inches='tight')
    plt.close()

print("Generating pairwise confusion matrix plots...")
for key, cm in all_confusion_matrices.items():
    if cm is None: 
        continue
    parts = key.split('_')
    mode_name, harm_name = parts[0], parts[1]
    pair_name = '_'.join(parts[2:])
    labels_pair = pair_name.split('_vs_')
    title = (f"Confusion Matrix - {mode_name.title()} - {harm_name.title()} - "
            f"{labels_pair[0]} vs {labels_pair[1]}\n"
            f"(Sum of {NUM_RUNS} runs) - "
            f"{classifier_name} - {thr_str} - {fs_str} - {resampling_str} - {cv_strategy}")
    save_path = os.path.join(
        saving_path,
        f"{classifier_name}_PAIRWISE_confusion_matrix_plot_{key}_{thr_str}{fs_str}{resampling_str}{CV_TAG}.png"
    )
    plot_confusion_matrix(cm, labels_pair, title, save_path)

# Selected features (LONG)
features_df = pd.DataFrame(selected_features_logs)
features_csv_path = os.path.join(
    saving_path,
    f"{classifier_name}_PAIRWISE_selected_features_{thr_str}{fs_str}{resampling_str}{CV_TAG}_pairwise.csv"
)
features_df.to_csv(features_csv_path, index=False)

# Selected features summary
summary_df = (
    features_df[features_df["feature"].isna()]
    .groupby(["mode","harmonization","pair"])["k_selected"]
    .agg(["mean","std","min","max"])
    .reset_index()
)
summary_csv_path = os.path.join(
    saving_path,
    f"{classifier_name}_PAIRWISE_selected_features_SUMMARY_{thr_str}{fs_str}{resampling_str}{CV_TAG}_summary.csv"
)
summary_df.to_csv(summary_csv_path, index=False)
# === SAVING SELECTED PARAMETERS ===
params_df = pd.DataFrame(param_logs)

# Keep consistency with your filenames; if you already use a CV_TAG, reuse it here:
# Example without CV_TAG; add it if already defined:

params_csv_path = os.path.join(
    saving_path,
    f"{classifier_name}_PAIRWISE_best_params_{thr_str}{fs_str}{resampling_str}_{cv_strategy}.csv"
)
params_df.to_csv(params_csv_path, index=False)

if "feature_importance_logs" in globals() and len(feature_importance_logs) > 0:
    fi_df = pd.concat(feature_importance_logs, ignore_index=True)
    fi_csv_path = os.path.join(
        saving_path,
        f"{classifier_name}_PAIRWISE_feature_importance_{thr_str}{fs_str}{resampling_str}{CV_TAG}.csv"
    )
    fi_df.to_csv(fi_csv_path, index=False)
    print(f"Feature importance: {fi_csv_path}")
print("Pairwise classification completed successfully!")
print(f"Results: {results_csv}")
print(f"Features: {features_csv_path}")
print(f"Summary : {summary_csv_path}")
print(f"Best-params : {params_csv_path}")
