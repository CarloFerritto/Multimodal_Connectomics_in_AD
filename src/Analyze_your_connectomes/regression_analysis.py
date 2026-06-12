import os
import warnings
import numpy as np
import pandas as pd

from tqdm import tqdm

from sklearn.model_selection import KFold, StratifiedKFold
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.cross_decomposition import PLSRegression
from sklearn.linear_model import ElasticNetCV
from sklearn.metrics import mean_squared_error, mean_absolute_error, r2_score

warnings.filterwarnings("ignore", category=UserWarning)

# =========================
# CONFIG
# =========================

EXCLUDE_PATIENTS = True

modes = [ "structural+functional","multilayer","structural", "functional"] #

outer_splits = 10
inner_splits = 5
random_state = 42

base_components_grid = list(range(1, 26))

# ElasticNetCV grid (same as the current version)
ENET_ALPHAS = np.array([1e-4, 1e-3, 1e-2, 1e-1, 1.0, 10, 100], dtype=float)
ENET_L1RATIOS = np.array([0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1], dtype=float)

# Threshold from environment variable
thr_str = os.environ.get("THR_STR", "dens0.25").strip()
if not thr_str:
    thr_str = "dens0.25"

# Machine
machine = "cluster" # "local" | "cluster"
if machine == "local":
    base_path = "Y:/project/shared/data/sample_dataset/analysis/harmonization_final/metrics"
    excel_path = "Y:/project/shared/data/sample_dataset/clinical_data/dataset_summary_only_ALL.xlsx"
    exclude_excel_path = "Y:/project/shared/data/sample_dataset/analysis/harmonization_final/data/subjects_harmonization.xlsx"
elif machine == "cluster":
    base_path = "/workspace/project/data/sample_dataset/multilayer/metrics"
    excel_path = "/workspace/project/data/sample_dataset/multilayer/dataset_summary_only_ALL.xlsx"
    exclude_excel_path = "/workspace/project/data/sample_dataset/multilayer/subjects_harmonization.xlsx"
else:
    raise ValueError("MACHINE must be 'local' or 'cluster'")

saving_path = os.path.join(base_path, "REG_PLS_ENET_NESTEDSTRATIFIEDCV", thr_str)
os.makedirs(saving_path, exist_ok=True)

# =========================
# Utilities
# =========================

def _compute_vip(pls_or_wrapper, X, Y):
    pls = getattr(pls_or_wrapper, "pls_", pls_or_wrapper)
    T = pls.x_scores_
    W = pls.x_weights_
    Q = pls.y_loadings_
    p, A = W.shape
    if A == 0:
        return np.zeros(p)

    SSY_t = np.sum(T**2, axis=0)
    SSY_q = np.sum(Q**2, axis=0)
    SSY = SSY_t * SSY_q

    w_norm2 = np.sum(W**2, axis=0)
    w_norm2[w_norm2 == 0] = 1.0

    num = (((W**2) / w_norm2[None, :]) @ SSY.reshape(-1, 1))
    den = SSY.sum() if SSY.sum() != 0 else 1.0
    vip = np.sqrt((p * num.ravel()) / den)
    return vip

def safe_components_grid(n_train, n_features):
    max_a = max(1, min(n_train - 1, n_features))
    grid = [a for a in base_components_grid if 1 <= a <= max_a]
    if len(grid) == 0:
        grid = [max_a]
    return grid

def corr_nan_safe(a, b):
    a = np.asarray(a).ravel()
    b = np.asarray(b).ravel()
    if np.std(a) < 1e-12 or np.std(b) < 1e-12:
        return np.nan
    return float(np.corrcoef(a, b)[0, 1])

def tune_pls_components_mse(X_train_df, y_train, inner_cv, comp_grid):
    """
    Select n_components by minimizing mean inner-CV MSE.
    """
    best_a = None
    best_mse = np.inf

    for a in comp_grid:
        mses = []
        for tr2, va2 in inner_cv.split(X_train_df):
            X_tr2 = X_train_df.iloc[tr2].to_numpy(float)
            y_tr2 = y_train[tr2]
            X_va2 = X_train_df.iloc[va2].to_numpy(float)
            y_va2 = y_train[va2]

            pipe = Pipeline([
                ("scaler", StandardScaler()),
                ("pls", PLSRegression(n_components=a, scale=False))
            ])
            pipe.fit(X_tr2, y_tr2)
            y_hat = pipe.predict(X_va2).ravel()
            mses.append(mean_squared_error(y_va2, y_hat))

        mean_mse = float(np.mean(mses))
        if mean_mse < best_mse:
            best_mse = mean_mse
            best_a = a

    return best_a, best_mse

def tune_pls_components_mse_stratified(X_train_df, y_train, y_train_bins, inner_cv, comp_grid):
    best_a = None
    best_mse = np.inf

    for a in comp_grid:
        mses = []
        for tr2, va2 in inner_cv.split(X_train_df, y_train_bins):
            X_tr2 = X_train_df.iloc[tr2].to_numpy(float)
            y_tr2 = y_train[tr2]
            X_va2 = X_train_df.iloc[va2].to_numpy(float)
            y_va2 = y_train[va2]

            pipe = Pipeline([
                ("scaler", StandardScaler()),
                ("pls", PLSRegression(n_components=a, scale=False))
            ])
            pipe.fit(X_tr2, y_tr2)
            y_hat = pipe.predict(X_va2).ravel()
            mses.append(mean_squared_error(y_va2, y_hat))

        mean_mse = float(np.mean(mses))
        if mean_mse < best_mse:
            best_mse = mean_mse
            best_a = a

    return best_a, best_mse


def load_mode_data(mode: str, thr_str: str):
    if mode == "structural+functional":
        struct_nodal_file = f"nodal_metrics_structural_combat_{thr_str}.csv"
        #struct_global_file = f"global_metrics_structural_combat_{thr_str}.csv"
        func_nodal_file = f"nodal_metrics_functional_combat_{thr_str}.csv"
        #func_global_file = f"global_metrics_functional_combat_{thr_str}.csv"

        struct_nodal_path = os.path.join(base_path, struct_nodal_file)
        #struct_global_path = os.path.join(base_path, struct_global_file)
        func_nodal_path = os.path.join(base_path, func_nodal_file)
        #func_global_path = os.path.join(base_path, func_global_file)

        if not all(map(os.path.exists, [struct_nodal_path, func_nodal_path])):
            return None

        struct_nodal_df = pd.read_csv(struct_nodal_path)
        #struct_global_df = pd.read_csv(struct_global_path)
        func_nodal_df = pd.read_csv(func_nodal_path)
        #func_global_df = pd.read_csv(func_global_path)

        #struct_data = pd.merge(struct_nodal_df, struct_global_df, on="id")
        #func_data = pd.merge(func_nodal_df, func_global_df, on="id")
    
        struct_features = struct_nodal_df.drop(columns=["id"]).add_prefix("struct_")
        func_features = func_nodal_df.drop(columns=["id"]).add_prefix("func_")

        combined = struct_nodal_df[["id"]].join(struct_features).merge(
            func_nodal_df[["id"]].join(func_features),
            on="id", how="inner"
        )
        return combined
    else:
        nodal_file = f"nodal_metrics_{mode}_combat_{thr_str}.csv"
        #global_file = f"global_metrics_{mode}_combat_{thr_str}.csv"

        nodal_path = os.path.join(base_path, nodal_file)
        #global_path = os.path.join(base_path, global_file)

        if not os.path.exists(nodal_path):
            return None

        nodal_df = pd.read_csv(nodal_path)
        #global_df = pd.read_csv(global_path)
        data = nodal_df.copy()
        return data

# =========================
# Load clinical + excluded
# =========================

info_full = pd.read_excel(excel_path, sheet_name="ALL")

required_cols = ["Subject", "Batch_variable_simplified"]
missing_req = [c for c in required_cols if c not in info_full.columns]
if missing_req:
    raise ValueError(f"Missing columns in excel file: {missing_req}")

info_full = info_full.rename(columns={
    "Subject": "id",
    "Batch_variable_simplified": "batch",
})
if "Diagnosis_amyloid" not in info_full.columns:
    raise ValueError("Missing Diagnosis_amyloid column in the excel file.")
if EXCLUDE_PATIENTS:
    exclude_df = pd.read_excel(exclude_excel_path, sheet_name="Sheet1")
    excluded_ids = set(exclude_df["Subject"].dropna().unique())
    print(f"FILTER ACTIVE: {len(excluded_ids)} excluded")
else:
    excluded_ids = set()
    print("FILTER DISABLED")

# Scores from environment variable
score_env = os.environ.get("SCORES", "MOCA,PET_CENTILOID,CDRSB,ADAS13,MMSE").strip() # MOCA, MMSE, ADAS13, CDRSB, PET_CENTILOID
score_cols = [s.strip() for s in score_env.split(",") if s.strip()]
if len(score_cols) == 0:
    raise ValueError("SCORES is empty: set env SCORES='ADAS13,CDRSB,...'")

print("thr_str:", thr_str)
print("Analyzed scores:", score_cols)

# =========================
# RUN nested CV
# =========================

metrics_rows = []
vip_rows = []
params_rows = []
enet_coef_rows = []
pred_rows = []
for mode in modes:
    print(f"\n==== MODE: {mode} ====")

    mode_data = load_mode_data(mode, thr_str=thr_str)
    if mode_data is None:
        print(f"[SKIP] Missing files for {mode} (thr={thr_str})")
        continue

    metric_cols = [c for c in mode_data.columns if c != "id"]

    for score_col in score_cols:
        if score_col not in info_full.columns:
            print(f"[SKIP] Score {score_col} not present in the excel file")
            continue

        info_min = info_full[["id", "batch", "Diagnosis_amyloid", score_col,]].copy()
        data = pd.merge(mode_data, info_min, on="id", how="inner")
        data = data[data["Diagnosis_amyloid"] != "MCI-"].reset_index(drop=True)
        if EXCLUDE_PATIENTS and len(excluded_ids) > 0:
            n_before = len(data)
            data = data[~data["id"].isin(excluded_ids)].reset_index(drop=True)
            
        data = data.dropna(subset=["batch", score_col,"Diagnosis_amyloid"]).reset_index(drop=True)
        if len(data) < outer_splits + 5:
            print(f"[SKIP] Too few subjects for score={score_col} mode={mode}: n={len(data)}")
            continue
        data = data.sort_values("id").reset_index(drop=True)
        # X = metrics only, numeric coercion
        X_all = data[metric_cols].copy().apply(pd.to_numeric, errors="coerce")
        X_all[metric_cols] = X_all[metric_cols].fillna(0)
        bad_cols = X_all.columns[X_all.isna().all(axis=0)].tolist()
        if bad_cols:
            print(f"[WARN] Dropping non-convertible columns: {bad_cols[:20]}{'...' if len(bad_cols)>20 else ''}")
            X_all = X_all.drop(columns=bad_cols)

        y_all = pd.to_numeric(data[score_col], errors="coerce")
        
        ok_rows = y_all.notna() & X_all.notna().all(axis=1)

        X_all = X_all.loc[ok_rows].reset_index(drop=True)
        y_all = y_all.loc[ok_rows].reset_index(drop=True)
        ids_all = data.loc[ok_rows, "id"].reset_index(drop=True)
        print(f"Subject filter ({score_col}): {n_before} -> {len(y_all)}")
        n_bins = 10  # tipico 5–10; con 280 e 5-fold, 10 di solito ok
        y_bins = pd.qcut(y_all, q=n_bins, labels=False, duplicates="drop")
        y_bins = y_bins.to_numpy()
        if len(X_all) < outer_splits + 5:
            print(f"[SKIP] Too few subjects after dropna for score={score_col} mode={mode}: n={len(X_all)}")
            continue

        outer_cv = StratifiedKFold(n_splits=outer_splits, shuffle=True, random_state=random_state)

        for fold_idx, (train_idx, test_idx) in enumerate(
        tqdm(outer_cv.split(X_all, y_bins), total=outer_splits, desc=f"{mode} | {score_col}", leave=False)
        ):
            X_train_df = X_all.iloc[train_idx].copy()
            X_test_df  = X_all.iloc[test_idx].copy()
            y_train = y_all.iloc[train_idx].to_numpy()
            y_test  = y_all.iloc[test_idx].to_numpy()

            #inner_cv = KFold(n_splits=inner_splits, shuffle=True, random_state=random_state + fold_idx + 1)
            n_bins_inner = 10  # puoi anche usare n_bins
            y_train_bins = pd.qcut(y_train, q=n_bins_inner, labels=False, duplicates="drop")

            # stratified inner CV
            inner_cv = StratifiedKFold(
                n_splits=inner_splits,
                shuffle=True,
                random_state=random_state + fold_idx + 1
            )
            # ===== PLS: tuning on MSE =====
            comp_grid = safe_components_grid(n_train=len(train_idx), n_features=X_train_df.shape[1])
            #best_k, best_inner_mse = tune_pls_components_mse(X_train_df, y_train, inner_cv, comp_grid)
            best_k, best_inner_mse = tune_pls_components_mse_stratified(X_train_df, y_train, y_train_bins, inner_cv, comp_grid)
            pls_model = Pipeline([
                ("scaler", StandardScaler()),
                ("pls", PLSRegression(n_components=best_k, scale=False))
            ])
            pls_model.fit(X_train_df.to_numpy(float), y_train)
            y_pred_pls = pls_model.predict(X_test_df.to_numpy(float)).ravel()

            mse_pls = float(mean_squared_error(y_test, y_pred_pls))
            rmse_pls = float(np.sqrt(mse_pls))
            mae_pls = float(mean_absolute_error(y_test, y_pred_pls))
            r2_pls = float(r2_score(y_test, y_pred_pls)) if np.std(y_test) > 1e-12 else np.nan
            corr_pls = corr_nan_safe(y_test, y_pred_pls)
            for sid, yt, yp in zip(ids_all.iloc[test_idx].tolist(), y_test, y_pred_pls):
                pred_rows.append({
                    "thr_str": thr_str,
                    "model": "PLS",
                    "mode": mode,
                    "score": score_col,
                    "fold": fold_idx,
                    "id": sid,
                    "y_true": float(yt),
                    "y_pred": float(yp)
                })
            metrics_rows.append({
                "thr_str": thr_str,
                "model": "PLS",
                "mode": mode,
                "score": score_col,
                "fold": fold_idx,
                "n_train": int(len(train_idx)),
                "n_test": int(len(test_idx)),
                "n_features": int(X_train_df.shape[1]),
                "inner_best_mse": float(best_inner_mse),
                "best_n_components": int(best_k),
                "mse_test": mse_pls,
                "rmse_test": rmse_pls,
                "mae_test": mae_pls,
                "r2_test": r2_pls,
                "corr_test": corr_pls
            })

            params_rows.append({
                "thr_str": thr_str,
                "model": "PLS",
                "mode": mode,
                "score": score_col,
                "fold": fold_idx,
                "best_n_components": int(best_k),
                "inner_best_mse": float(best_inner_mse),
            })

            # PLS VIP
            scaler_pls = pls_model.named_steps["scaler"]
            pls_step = pls_model.named_steps["pls"]
            X_train_scaled = scaler_pls.transform(X_train_df.to_numpy(float))
            vip = _compute_vip(pls_step, X_train_scaled, y_train.reshape(-1, 1)).ravel()

            feat_names = list(X_train_df.columns)
            if len(vip) == len(feat_names):
                for f, v in zip(feat_names, vip):
                    vip_rows.append({
                        "thr_str": thr_str,
                        "model": "PLS",
                        "mode": mode,
                        "score": score_col,
                        "fold": fold_idx,
                        "param": f"n_components={int(best_k)}",
                        "feature": f,
                        "importance": float(v),
                        "importance_type": "VIP"
                    })

            # ===== ElasticNetCV: tuning on MSE ====

            inner_splits_list = list(inner_cv.split(X_train_df, y_train_bins))

            enet_model = Pipeline([
                ("scaler", StandardScaler()),
                ("enetcv", ElasticNetCV(
                    l1_ratio=ENET_L1RATIOS,
                    alphas=ENET_ALPHAS,
                    cv=inner_splits_list,   # <- qui la stratificazione è rispettata
                    max_iter=20000,
                    n_jobs=-1,
                ))
            ])
            '''
            enet_model = Pipeline([
                ("scaler", StandardScaler()),
                ("enetcv", ElasticNetCV(
                    l1_ratio=ENET_L1RATIOS,
                    alphas=ENET_ALPHAS,
                    cv=inner_cv,
                    max_iter=20000,
                    n_jobs=-1,
                ))
            ])
            '''
            enet_model.fit(X_train_df.to_numpy(float), y_train)

            enetcv = enet_model.named_steps["enetcv"]
            y_pred_enet = enet_model.predict(X_test_df.to_numpy(float)).ravel()

            mse_en = float(mean_squared_error(y_test, y_pred_enet))
            rmse_en = float(np.sqrt(mse_en))
            mae_en = float(mean_absolute_error(y_test, y_pred_enet))
            r2_en = float(r2_score(y_test, y_pred_enet)) if np.std(y_test) > 1e-12 else np.nan
            corr_en = corr_nan_safe(y_test, y_pred_enet)
            for sid, yt, yp in zip(ids_all.iloc[test_idx].tolist(), y_test, y_pred_enet):
                pred_rows.append({
                    "thr_str": thr_str,
                    "model": "ElasticNetCV",
                    "mode": mode,
                    "score": score_col,
                    "fold": fold_idx,
                    "id": sid,
                    "y_true": float(yt),
                    "y_pred": float(yp)
                })
            metrics_rows.append({
                "thr_str": thr_str,
                "model": "ElasticNetCV",
                "mode": mode,
                "score": score_col,
                "fold": fold_idx,
                "n_train": int(len(train_idx)),
                "n_test": int(len(test_idx)),
                "n_features": int(X_train_df.shape[1]),
                "inner_best_mse": float(np.min(np.mean(enetcv.mse_path_, axis=2))),  # mean over folds, then min over grid
                "best_alpha": float(enetcv.alpha_),
                "best_l1_ratio": float(enetcv.l1_ratio_),
                "mse_test": mse_en,
                "rmse_test": rmse_en,
                "mae_test": mae_en,
                "r2_test": r2_en,
                "corr_test": corr_en
            })

            params_rows.append({
                "thr_str": thr_str,
                "model": "ElasticNetCV",
                "mode": mode,
                "score": score_col,
                "fold": fold_idx,
                "best_alpha": float(enetcv.alpha_),
                "best_l1_ratio": float(enetcv.l1_ratio_),
            })

            # ElasticNetCV coefficients (on scaled X)
            coef = np.asarray(enetcv.coef_).ravel()
            for f, b in zip(feat_names, coef):
                enet_coef_rows.append({
                    "thr_str": thr_str,
                    "model": "ElasticNetCV",
                    "mode": mode,
                    "score": score_col,
                    "fold": fold_idx,
                    "alpha": float(enetcv.alpha_),
                    "l1_ratio": float(enetcv.l1_ratio_),
                    "feature": f,
                    "coef": float(b),
                    "abs_coef": float(abs(b))
                })

# =========================
# SAVE
# =========================
metrics_df = pd.DataFrame(metrics_rows)
vip_df = pd.DataFrame(vip_rows)
params_df = pd.DataFrame(params_rows)
enet_coef_df = pd.DataFrame(enet_coef_rows)
metrics_df.to_csv(os.path.join(saving_path, f"metrics_{thr_str}.csv"), index=False)
vip_df.to_csv(os.path.join(saving_path, f"pls_vip_{thr_str}.csv"), index=False)
params_df.to_csv(os.path.join(saving_path, f"bestparams_{thr_str}.csv"), index=False)
enet_coef_df.to_csv(os.path.join(saving_path, f"enet_coefs_{thr_str}.csv"), index=False)
# =========================
# SUMMARY (mean ± std)
# =========================

# numeric metrics to summarize (ignored if a column is missing)
metric_cols_summary = [c for c in ["mse_test", "rmse_test", "mae_test", "r2_test", "corr_test"] if c in metrics_df.columns]

summary_df = (
    metrics_df
    .groupby(["thr_str", "model", "mode", "score"], as_index=False)[metric_cols_summary]
    .agg(["mean", "std"])
)

# Flatten MultiIndex columns: mse_test_mean, mse_test_std, ...
summary_df.columns = [
    f"{col}_{stat}" if stat else col
    for col, stat in summary_df.columns.to_flat_index()
]

summary_csv = os.path.join(saving_path, f"summary_{thr_str}.csv")
summary_df.to_csv(summary_csv, index=False)
pred_df = pd.DataFrame(pred_rows)
pred_df.to_csv(os.path.join(saving_path, f"predictions_{thr_str}.csv"), index=False)
print("DONE")
print("Saved in:", saving_path)
