import os
import pandas as pd
import numpy as np
import networkx as nx
from nilearn import plotting
from scipy.stats import mode
from sklearn.preprocessing import MinMaxScaler
from concurrent.futures import ProcessPoolExecutor
from tqdm import tqdm

# === CONFIG ===
machine = "local"  # "local" or "grid"
modality = os.getenv("MODALITY", "multilayer")  # "functional", "structural", "multilayer"
if modality == "all":
    modes = ["multilayer","functional","structural"]  # "functional", "structural", "multilayer"
else:
    modes = [modality]

# Read from environment variables, with reasonable defaults
threshold_strategy = os.getenv("THRESHOLD_STRATEGY", "consistency")  # "consistency" or "proportional_edges"
threshold_ratio = float(os.getenv("THRESHOLD_RATIO", "0.2"))         # used if strategy == "consistency"
density = float(os.getenv("DENSITY", "0.5"))                         # used if strategy == "proportional_edges"

harm_env = os.getenv("HARMONIZATION", "")  # "noharm", "combat", "covbat"
harmonization = None if harm_env == "" or harm_env.lower() == "noharm" else harm_env  # None / "combat" / "covbat"

NUM_REGIONS = 416

FUNC_FILENAME = "func_new.csv" if harmonization is None else f"func_{harmonization}.csv"
STRU_FILENAME = "struct_streamlines_log_new.csv" if harmonization is None else f"struct_streamlines_log_{harmonization}.csv"
STRU_ORIG_FILENAME = "struct_streamlines_new.csv"
if machine == "local":
    BASE_PATH = "C:/project/data/harmonized_multilayer/data"
    EXCEL_PATH ="Y:/project/shared/data/sample_dataset/clinical_data/dataset_summary.xlsx"
    REF_PATH = os.path.join(BASE_PATH, "subjects_harmonization.xlsx")
    FUNC_PATH = os.path.join(BASE_PATH,FUNC_FILENAME)
    STRU_PATH = os.path.join(BASE_PATH,STRU_FILENAME)
    STRU_ORIG_PATH = os.path.join(BASE_PATH,STRU_ORIG_FILENAME)
    SAVE_PATH = os.path.join(BASE_PATH,"metrics")
elif machine=="grid":
    BASE_PATH = "/workspace/project/data/sample_dataset/multilayer"
    EXCEL_PATH = "/workspace/project/data/sample_dataset/multilayer/dataset_summary.xlsx"
    REF_PATH = os.path.join(BASE_PATH, "subjects_harmonization.xlsx")
    FUNC_PATH = os.path.join(BASE_PATH,FUNC_FILENAME)
    STRU_PATH = os.path.join(BASE_PATH,STRU_FILENAME)
    STRU_ORIG_PATH = os.path.join(BASE_PATH,STRU_ORIG_FILENAME)
    SAVE_PATH = os.path.join(BASE_PATH,"metrics")
os.makedirs(SAVE_PATH, exist_ok=True)

# === UTILS ===
def load_matrix_excel(path, mode):
    df = pd.read_csv(path)
    ids = df.columns
    matrices = {}
    for col in ids:
        mat = np.zeros((NUM_REGIONS, NUM_REGIONS))
        iu = np.triu_indices(NUM_REGIONS, 1)
        mat[iu] = df[col].values
        mat += mat.T
        if mode == "functional":
            mat[mat < 0] = 0  # Zero negative values
        scaler = MinMaxScaler()
        mat_flat = mat[np.triu_indices(NUM_REGIONS, k=1)]
        mat_flat_scaled = scaler.fit_transform(mat_flat.reshape(-1, 1)).flatten()
        mat_scaled = np.zeros_like(mat)
        mat_scaled[iu] = mat_flat_scaled
        mat_scaled += mat_scaled.T
        np.fill_diagonal(mat_scaled, 0)  # Ensure diagonal is zero after scaling
        matrices[col] = mat_scaled
    return matrices
def consistency_mask_CN(matrices, subject_info, ref_subjects, threshold_ratio):
    """
    MODIFICATION: Use only CN- subjects that are present in both subject_info and ref_subjects
    """
    # Filter only CN- subjects present in subjects_harmonization.xlsx
    CN_ids_all = subject_info[subject_info["Diagnosis_amyloid"] == "CN-"]["Subject"].values
    CN_ids = [pid for pid in CN_ids_all if pid in ref_subjects]  # NEW FILTER
    
    # print(f"Total CN-: {len(CN_ids_all)}, CN- in reference: {len(CN_ids)}")
    
    CN_mats = [matrices[pid] for pid in CN_ids if pid in matrices]
    CN_stack = np.stack(CN_mats)
    mean_mat = np.nanmean(CN_stack, axis=0)
    std_mat = np.nanstd(CN_stack, axis=0)
    
    with np.errstate(divide='ignore', invalid='ignore'):
        consistency = std_mat / mean_mat
        consistency[np.isnan(consistency)] = np.inf
    
    flat = consistency[np.triu_indices_from(consistency, 1)]
    thresh = np.percentile(flat, 100 * threshold_ratio)
    mask = consistency <= thresh
    
    print(f"Sparsity after threshold: {1 - np.sum(mask) / mask.size:.3f}")
    return mask

def proportional_threshold_edges(matrix, density: float):
    """
    Optimized version using argpartition instead of argsort.
    Complexity: O(n) instead of O(n log n)
    """
    assert 0 < density <= 1, "density must be in (0,1]"
    
    n = matrix.shape[0]
    iu = np.triu_indices(n, 1)
    w = matrix[iu]
    
    m_total = w.size
    k = int(np.round(density * m_total))
    k = max(0, min(k, m_total))
    
    if k == 0:
        out = np.zeros_like(matrix)
        np.fill_diagonal(out, 0)
        return out
    
    # Use argpartition instead of argsort - much faster!
    # Find the indices of the k largest elements
    kth = m_total - k  # k-esimo elemento (da destra)
    keep_idx = np.argpartition(w, kth)[kth:]  # top-k elementi
    
    out = np.zeros_like(matrix)
    out[iu[0][keep_idx], iu[1][keep_idx]] = w[keep_idx]
    out = out + out.T
    np.fill_diagonal(out, 0)
    
    return out
def apply_mask(matrix, mask):
    return matrix * mask

def plot_mask(mask, title="Threshold Mask"):
    plotting.plot_matrix(mask.astype(float), figure=(8, 6), labels=None, vmax=1.0, title=title)

def compute_M_max(M):
    num_nodes_SG = M.shape[0]
    n_nodes = num_nodes_SG // 2
    M_max = np.zeros((n_nodes, n_nodes), dtype=float)
    for source in range(n_nodes):
        for target in range(n_nodes):
            tmp = np.array([
                M[source, target],
                M[source, target + n_nodes],
                M[source + n_nodes, target],
                M[source + n_nodes, target + n_nodes]
            ])
            tmp = tmp[np.nonzero(tmp)]
            M_max[source, target] = np.max(tmp) if tmp.size > 0 else 0
    np.fill_diagonal(M_max, 0)
    return M_max

def participation_coefficient_multilayer(func_mat, stru_mat):
    strenght_func = np.sum(func_mat, axis=1)
    strenght_stru = np.sum(stru_mat, axis=1)
    total_strenght= strenght_stru+ strenght_func
    coeff = np.zeros_like(total_strenght, dtype=float)
    for i in range(len(total_strenght)):
        if total_strenght[i] > 0:
            p_func = strenght_func[i] / total_strenght[i]
            p_stru = strenght_stru[i] / total_strenght[i]
            coeff[i] = 2*(1 - (p_func ** 2 + p_stru ** 2))
    return dict(enumerate(coeff))

def multilayer_clustering_coefficient(func_mat, stru_mat):
    N = func_mat.shape[0]
    M = 2
    clustering = {}
    for i in range(N):
        numerator = 0
        denominator = 0
        for alpha in [0, 1]:
            A_alpha = func_mat if alpha == 0 else stru_mat
            k_alpha = np.sum(A_alpha[i] > 0)
            denominator += k_alpha * (k_alpha - 1)
            for beta in [0, 1]:
                if beta == alpha:
                    continue
                A_beta = func_mat if beta == 0 else stru_mat
                for j in range(N):
                    for m in range(N):
                        if i != m and i != j and j != m:
                            val = (A_alpha[i, j] * A_beta[j, m] * A_alpha[m, i]) ** (1 / 3)
                            numerator += val
        denom = (M - 1) * denominator if denominator > 0 else 1
        clustering[i] = numerator / denom
    return clustering
def compute_multilayer_clustering_for_node(args):
    i, func_mat, stru_mat = args
    N = func_mat.shape[0]
    M = 2
    numerator = 0
    denominator = 0
    for alpha in [0, 1]:
        A_alpha = func_mat if alpha == 0 else stru_mat
        k_alpha = np.sum(A_alpha[i] > 0)
        denominator += k_alpha * (k_alpha - 1)
        for beta in [0, 1]:
            if beta == alpha:
                continue
            A_beta = func_mat if beta == 0 else stru_mat
            for j in range(N):
                for m in range(N):
                    if i != m and i != j and j != m:
                        val = (A_alpha[i, j] * A_beta[j, m] * A_alpha[m, i]) ** (1 / 3)
                        numerator += val
    denom = (M - 1) * denominator if denominator > 0 else 1
    return (i, numerator / denom)

def compute_graph_metrics(matrix):
    G = nx.from_numpy_array(matrix)
    inv_weights = {(u, v): 1 / d["weight"] if d["weight"] > 0 else 1e6 for u, v, d in G.edges(data=True)}
    nx.set_edge_attributes(G, inv_weights, name="inv_weight")
    try:
        avg_path = nx.average_shortest_path_length(G, weight="inv_weight") if nx.is_connected(G) else np.nan
    except:
        avg_path = np.nan
    metrics = {
        "strength": dict(G.degree(weight="weight")),
        "betweenness": nx.betweenness_centrality(G, weight="inv_weight", normalized=True),
        "closeness": nx.closeness_centrality(G, distance="inv_weight"),
        "clustering": nx.clustering(G, weight="weight"),
        "local_efficiency": local_efficiency_weighted_parallel(matrix,max_workers=64),
        "assortativity": nx.degree_pearson_correlation_coefficient(G),
        "avg_clustering": nx.average_clustering(G, weight="weight"),
        "avg_shortest_path": avg_path,
        "global_efficiency": global_efficiency(G),
        "transitivity": nx.transitivity(G)
    }
    return metrics

def compute_multilayer_metrics(func_mat, stru_mat):
    combined_strength = dict(enumerate(np.sum(func_mat + stru_mat, axis=1)))
    part_coeff = participation_coefficient_multilayer(func_mat, stru_mat)
    clustering = multilayer_clustering_coefficient_parallel(func_mat, stru_mat,max_workers=64)
    combined = np.block([
        [func_mat, np.identity(NUM_REGIONS)],
        [np.identity(NUM_REGIONS), stru_mat]
    ])
    Mmax = compute_M_max(combined)
    G = nx.from_numpy_array(Mmax)
    inv_weights = {(u, v): 1 / d["weight"] if d["weight"] > 0 else 1e6 for u, v, d in G.edges(data=True)}
    nx.set_edge_attributes(G, inv_weights, name="inv_weight")
    try:
        avg_path = nx.average_shortest_path_length(G, weight="inv_weight") if nx.is_connected(G) else np.nan
    except:
        avg_path = np.nan

    metrics = {
        "strength": combined_strength,
        "participation": part_coeff,
        "clustering": clustering,
        "betweenness": nx.betweenness_centrality(G, weight="inv_weight", normalized=True),
        "closeness": nx.closeness_centrality(G, distance="inv_weight"),
        "local_efficiency": local_efficiency_weighted_parallel(Mmax,max_workers=64),
        "assortativity": nx.degree_pearson_correlation_coefficient(G),
        "avg_clustering": np.mean(list(clustering.values())),
        "avg_shortest_path": avg_path,
        "global_efficiency": global_efficiency(G),
        "transitivity": nx.transitivity(G)
    }
    return metrics
def multilayer_clustering_coefficient_parallel(func_mat, stru_mat, max_workers=None):
    N = func_mat.shape[0]
    args = [(i, func_mat, stru_mat) for i in range(N)]
    clustering = {}
    with ProcessPoolExecutor(max_workers=max_workers) as executor:
        for i, val in executor.map(compute_multilayer_clustering_for_node, args):
            clustering[i] = val
    return clustering


def local_efficiency_weighted_parallel(matrix, max_workers=None):
    nodes = range(matrix.shape[0])
    args = [(n, matrix) for n in nodes]

    local_eff = {}
    with ProcessPoolExecutor(max_workers=max_workers) as executor:
        for node, eff in executor.map(compute_local_efficiency_for_node, args):
            local_eff[node] = eff

    return local_eff

def compute_local_efficiency_for_node(args):
    node, matrix = args
    G = nx.from_numpy_array(matrix)
    
    # Invert weights
    inv_weights = {(u, v): 1 / d["weight"] if d["weight"] > 0 else 1e6
                   for u, v, d in G.edges(data=True)}
    nx.set_edge_attributes(G, inv_weights, name="inv_weight")

    neighbors = list(G.neighbors(node))
    if len(neighbors) < 2:
        return (node, 0.0)

    subgraph = G.subgraph(neighbors).copy()
    path_lengths = dict(nx.all_pairs_dijkstra_path_length(subgraph, weight="inv_weight"))

    total = 0.0
    count = 0

    for u in neighbors:
        for v in neighbors:
            if u != v:
                try:
                    d = path_lengths[u][v]
                    total += 1 / d
                    count += 1
                except KeyError:
                    continue

    return (node, total / count if count > 0 else 0.0)

def global_efficiency(G):
    nodes = list(G.nodes)
    path_lengths = dict(nx.all_pairs_dijkstra_path_length(G, weight="inv_weight"))
    total = 0
    count = 0
    for u in nodes:
        for v in nodes:
            if u != v:
                try:
                    d = path_lengths[u][v]
                    total += 1 / d
                    count += 1
                except:
                    continue
    return total / count if count > 0 else 0.0

def build_multilayer_matrix(func_matrices, stru_matrices, func_mask, stru_mask):
    multilayer_matrices = {}
    for pid in func_matrices:
        if pid in stru_matrices:
            func = apply_mask(func_matrices[pid], func_mask)
            stru = apply_mask(stru_matrices[pid], stru_mask)
            combined = np.block([
                [func, np.identity(NUM_REGIONS)],
                [np.identity(NUM_REGIONS), stru]
            ])
            multilayer_matrices[pid] = combined
    return multilayer_matrices

def restore_zeros(matrices, original_matrices):
    for k in matrices:
        matrices[k][original_matrices[k] == 0] = 0
    return matrices
# === MAIN ===

def main():
    ref_df = pd.read_excel(REF_PATH)
    ref_subjects = set(ref_df['Subject'].values)
    for mode in modes:
        info = pd.read_excel(EXCEL_PATH, sheet_name="ALL")
        if mode == "functional":
            matrices = load_matrix_excel(FUNC_PATH, mode="functional")
        elif mode == "structural":
            matrices = load_matrix_excel(STRU_PATH, mode="structural")
            matrices_orig= load_matrix_excel(STRU_ORIG_PATH, mode="structural")
        elif mode == "multilayer":
            func_matrices = load_matrix_excel(FUNC_PATH, mode="functional")
            stru_matrices = load_matrix_excel(STRU_PATH, mode="structural")
            stru_matrices_orig= load_matrix_excel(STRU_ORIG_PATH, mode="structural")
            matrices = {}  # filled after mask creation
        else:
            raise ValueError("Invalid mode")

        if mode == "multilayer":
            if threshold_strategy == "consistency":
                func_mask = consistency_mask_CN(func_matrices, info, ref_subjects, threshold_ratio)
                stru_mask = consistency_mask_CN(stru_matrices, info, ref_subjects, threshold_ratio)
                plot_mask(func_mask, title=f"func_{harmonization or 'noharm'}_mask_thr_{str(threshold_ratio)}")
                plot_mask(stru_mask, title=f"struct_{harmonization or 'noharm'}_mask_thr_{str(threshold_ratio)}")
            # in any case, first restore structural zeros
            stru_matrices = restore_zeros(stru_matrices, stru_matrices_orig)

        elif mode == "structural":
            if threshold_strategy == "consistency":
                mask = consistency_mask_CN(matrices, info, ref_subjects, threshold_ratio)
                plot_mask(mask, title=f"{mode}_{harmonization or 'noharm'}_mask_thr_{str(threshold_ratio)}")
            matrices = restore_zeros(matrices, matrices_orig)

        elif mode == "functional":
            if threshold_strategy == "consistency":
                mask = consistency_mask_CN(matrices, info, ref_subjects, threshold_ratio)
                plot_mask(mask, title=f"{mode}_{harmonization or 'noharm'}_mask_thr_{str(threshold_ratio)}")

        all_metrics = []
        print(f"metrics extraction for {mode} modality and {(harmonization or 'NO')} harmonization, "
              f"strategy={threshold_strategy}, thr={threshold_ratio}, dens={density}")
        if mode != "multilayer":
            for pid in tqdm(matrices):
                try:
                    mat = matrices[pid]
                    if threshold_strategy == "consistency":
                        mat_eff = apply_mask(mat, mask)

                    elif threshold_strategy == "proportional_edges":
                        mat_eff = proportional_threshold_edges(mat, density)
                    else:
                        raise ValueError("threshold_strategy must be 'consistency' or 'proportional_edges'")
                    metrics = compute_graph_metrics(mat_eff)
                    metrics["id"] = pid
                    all_metrics.append(metrics)
                except Exception as e:
                    print(f"Error in {pid}: {e}")
        else:
            for pid in tqdm(func_matrices):
                try:
                    fmat = func_matrices[pid]
                    smat = stru_matrices[pid]
                    if threshold_strategy == "consistency":
                        f_eff = apply_mask(fmat, func_mask)
                        s_eff = apply_mask(smat, stru_mask)
                    elif threshold_strategy == "proportional_edges":
                        f_eff = proportional_threshold_edges(fmat, density)
                        s_eff = proportional_threshold_edges(smat, density)
                    else:
                        raise ValueError("threshold_strategy must be 'consistency' or 'proportional_edges'")
                    metrics = compute_multilayer_metrics(f_eff, s_eff)
                    metrics["id"] = pid
                    all_metrics.append(metrics)
                except Exception as e:
                    print(f"Error in {pid}: {e}")


        if mode == "multilayer":
            df_nodes = pd.DataFrame([{
                "id": m["id"],
                **{f"{k}_{n}": v for k in ["strength", "betweenness", "closeness", "clustering", "local_efficiency", "participation"] if k in m for n, v in m[k].items()}
            } for m in all_metrics])
        else:
            df_nodes = pd.DataFrame([{  
                "id": m["id"],
                **{f"{k}_{n}": v for k in ["strength", "betweenness", "closeness", "clustering", "local_efficiency"] if k in m for n, v in m[k].items()}
            } for m in all_metrics])

        df_global = pd.DataFrame([{
            "id": m["id"],
            "assortativity": m["assortativity"],
            "avg_clustering": m["avg_clustering"],
            "global_efficiency": m["global_efficiency"],
            "avg_shortest_path": m["avg_shortest_path"],
            "transitivity": m["transitivity"]
        } for m in all_metrics])

        harm_str = harmonization if harmonization else "noharm"
        if threshold_strategy == "consistency":
            tag = f"thr{str(threshold_ratio)}"   # identico a prima
        else:
            tag = f"dens{str(density)}"          # nuova strategia

        df_nodes.to_pickle(os.path.join(SAVE_PATH, f"nodal_metrics_{mode}_{harm_str}_{tag}.pkl"))
        df_global.to_pickle(os.path.join(SAVE_PATH, f"global_metrics_{mode}_{harm_str}_{tag}.pkl"))

        df_nodes.to_csv(os.path.join(SAVE_PATH, f"nodal_metrics_{mode}_{harm_str}_{tag}.csv"), index=False)
        df_global.to_csv(os.path.join(SAVE_PATH, f"global_metrics_{mode}_{harm_str}_{tag}.csv"), index=False)

if __name__ == "__main__":
    main()
