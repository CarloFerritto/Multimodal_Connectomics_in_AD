import numpy as np
import pandas as pd
from scipy.stats import entropy
from sklearn.model_selection import StratifiedShuffleSplit

import numpy as np
import pandas as pd
from scipy.stats import entropy

def select_calibration_set_sa_global(
    df,
    batch_col="Batch_variable_simplified",
    diagnosis_col="Diagnosis",
    age_col="Age",
    sex_col="Sex",
    cn_label="CN-",
    min_subjects_per_batch=2,
    safe_threshold=5,
    max_batch_contribution=0.6,
    n_age_bins=10,
    target_total=None,            # if None: use only min/max per batch + optional size penalty
    lambda_batch=0.2,             # batch penalty weight (soft)
    lambda_size=0.0,              # if target_total=None, set >0 to prefer smaller sets
    iterations=5000,
    T_start=1.0,
    T_min=1e-3,
    alpha=0.995,
    random_seed=42,
    verbose_every=500,
):

    rng = np.random.default_rng(random_seed)

    # --- 1) Filter CN- and reset index (to work with internal 0..N-1 indices) ---
    df_cn = df[df[diagnosis_col] == cn_label].copy()
    df_cn = df_cn.reset_index(drop=False)  # keep original index in 'index'
    orig_index = df_cn["index"].values

    if df_cn.empty:
        raise ValueError("No CN- subjects found.")

    # --- 2) Prepare batch -> list of local indices ---
    batches = {}
    excluded_batches = []
    for b, g in df_cn.groupby(batch_col):
        idx = g.index.to_numpy()
        if len(idx) < min_subjects_per_batch:
            excluded_batches.append((b, len(idx)))
            continue
        batches[b] = idx

    if not batches:
        raise ValueError("No eligible batch (all < min_subjects_per_batch).")

    # Min/max constraints per batch
    min_b = {}
    max_b = {}
    for b, idx in batches.items():
        n_b = len(idx)
        min_b[b] = min(safe_threshold, n_b)
        max_b[b] = max(min_b[b], int(np.floor(max_batch_contribution * n_b)))

    # --- 3) Feature arrays for objective ---
    ages = df_cn[age_col].astype(float).to_numpy()
    sexes_raw = df_cn[sex_col].astype(str).to_numpy()

    # Map sex to 0/1 (adjust if different encodings are used)
    # accepts 'M','F' and also '0','1' etc. if you want to extend
    sex_map = {"M": 0, "F": 1}
    if not set(np.unique(sexes_raw)).issubset(set(sex_map.keys())):
        raise ValueError(f"Sex must be in {list(sex_map.keys())}, found: {np.unique(sexes_raw)}")
    sexes = np.vectorize(sex_map.get)(sexes_raw)

    # Global age bins (uniform range)
    age_min, age_max = float(np.min(ages)), float(np.max(ages))
    if age_min == age_max:
        # edge case: all ages equal
        age_edges = np.array([age_min, age_max + 1e-6])
        n_age_bins_eff = 1
    else:
        age_edges = np.linspace(age_min, age_max, n_age_bins + 1)
        n_age_bins_eff = n_age_bins

    # Feasibility matrix for joint distribution: cells (age_bin, sex) existing in pool
    age_bin_all = np.clip(np.digitize(ages, age_edges, right=False) - 1, 0, n_age_bins_eff - 1)
    feasible = np.zeros((n_age_bins_eff, 2), dtype=int)
    for ab, sx in zip(age_bin_all, sexes):
        feasible[ab, sx] += 1
    feasible_mask = feasible > 0
    n_feasible_cells = int(feasible_mask.sum())

    # Uniform target over feasible cells
    target_joint = np.zeros_like(feasible, dtype=float)
    target_joint[feasible_mask] = 1.0 / n_feasible_cells

    # --- 4) Helper to count selected subjects per batch ---
    # batch_of[i] = batch label for subject i (local index)
    batch_of = df_cn[batch_col].to_numpy()

    def counts_by_batch(selected_set):
        c = {b: 0 for b in batches}
        for i in selected_set:
            b = batch_of[i]
            if b in c:
                c[b] += 1
        return c

    # --- 5) Objective function ---
    def objective(selected_set):
        sel = np.fromiter(selected_set, dtype=int)
        if sel.size == 0:
            return np.inf

        # Joint histogram (age_bin x sex) over selected subjects
        ab = age_bin_all[sel]
        sx = sexes[sel]
        joint = np.zeros((n_age_bins_eff, 2), dtype=float)
        for a, s in zip(ab, sx):
            joint[a, s] += 1.0
        joint = joint / joint.sum()

        # KL only on feasible cells (avoids chasing impossible combinations)
        p = joint[feasible_mask] + 1e-12
        q = target_joint[feasible_mask] + 1e-12
        kl_joint = entropy(p, q)

        # Batch penalty (soft): KL(p_batch || uniform)
        cb = counts_by_batch(selected_set)
        pb = np.array([cb[b] for b in batches], dtype=float)
        pb = pb / pb.sum()
        ub = np.ones_like(pb) / len(pb)
        kl_batch = entropy(pb + 1e-12, ub + 1e-12)

        # Size penalty (optional)
        if target_total is not None:
            # if target_total is fixed, this is usually not needed; leaving as 0
            size_pen = 0.0
        else:
            # pushes toward smaller set (if lambda_size>0)
            size_pen = lambda_size * (len(selected_set) / len(df_cn))

        return kl_joint + lambda_batch * kl_batch + size_pen

    # --- 6) Initialization: build a set satisfying min per batch ---
    selected = set()
    for b, idx in batches.items():
        k = min_b[b]
        chosen = rng.choice(idx, size=k, replace=False)
        selected.update(chosen.tolist())

    # If target_total is fixed, add random subjects up to target (without violating max_b)
    if target_total is not None:
        if target_total < len(selected):
            raise ValueError(f"target_total={target_total} < minimum required={len(selected)}")
        # candidate pool for add
        all_candidates = set(np.concatenate(list(batches.values())).tolist())
        while len(selected) < target_total:
            # pick a candidate that does not violate max for its batch
            remaining = np.array(list(all_candidates - selected), dtype=int)
            if remaining.size == 0:
                break
            i = int(rng.choice(remaining))
            b = batch_of[i]
            # check max
            cb = counts_by_batch(selected)
            if cb[b] < max_b[b]:
                selected.add(i)
            else:
                # if batch is full, try another
                continue

    current_score = objective(selected)
    best = set(selected)
    best_score = current_score

    # --- 7) SA moves: swap / add / drop ---
    def propose_move(sel_set):
        cb = counts_by_batch(sel_set)

        # possible move types based on constraints
        can_add = any(cb[b] < max_b[b] for b in batches)
        can_drop = any(cb[b] > min_b[b] for b in batches)
        # swap always possible if at least one candidate is not selected
        all_candidates = set(np.concatenate(list(batches.values())).tolist())
        can_swap = len(all_candidates - sel_set) > 0 and len(sel_set) > 0

        moves = []
        probs = []

        if can_swap:
            moves.append("swap")
            probs.append(0.6)
        if target_total is None and can_add:
            moves.append("add")
            probs.append(0.2)
        if target_total is None and can_drop:
            moves.append("drop")
            probs.append(0.2)

        probs = np.array(probs, dtype=float)
        probs = probs / probs.sum()
        move = rng.choice(moves, p=probs)

        new_set = set(sel_set)

        all_candidates_arr = np.concatenate(list(batches.values()))
        all_candidates_set = set(all_candidates_arr.tolist())

        if move == "swap":
            # remove one and add one (respecting min/max)
            out = int(rng.choice(list(new_set)))
            b_out = batch_of[out]

            # if removing out would violate min, retry with a valid out
            valid_out = [i for i in new_set if cb[batch_of[i]] > min_b[batch_of[i]]]
            if not valid_out:
                return sel_set  # no removable subject, no move
            out = int(rng.choice(valid_out))
            b_out = batch_of[out]

            # pick a non-selected subject that does not exceed max for its batch
            not_sel = list(all_candidates_set - new_set)
            rng.shuffle(not_sel)
            chosen_in = None
            for inn in not_sel:
                b_in = batch_of[inn]
                # after removing out, cb changes:
                cb_tmp = cb.copy()
                cb_tmp[b_out] -= 1
                if cb_tmp[b_in] < max_b[b_in]:
                    chosen_in = inn
                    break
            if chosen_in is None:
                return sel_set

            new_set.remove(out)
            new_set.add(int(chosen_in))
            return new_set

        if move == "add":
            # add a subject from a batch not at max
            candidates = [i for i in (all_candidates_set - new_set) if cb[batch_of[i]] < max_b[batch_of[i]]]
            if not candidates:
                return sel_set
            new_set.add(int(rng.choice(candidates)))
            return new_set

        if move == "drop":
            # remove a subject from a batch above minimum
            removable = [i for i in new_set if cb[batch_of[i]] > min_b[batch_of[i]]]
            if not removable:
                return sel_set
            new_set.remove(int(rng.choice(removable)))
            return new_set

        return sel_set

    # --- 8) SA loop ---
    T = T_start
    for it in range(1, iterations + 1):
        candidate = propose_move(selected)
        cand_score = objective(candidate)

        # accept?
        if cand_score < current_score or rng.random() < np.exp(-(cand_score - current_score) / max(T, 1e-12)):
            selected = candidate
            current_score = cand_score
            if current_score < best_score:
                best = set(selected)
                best_score = current_score

        T = max(T_min, T * alpha)

        if verbose_every and it % verbose_every == 0:
            print(f"Iter {it}/{iterations}  T={T:.4g}  score={current_score:.6f}  best={best_score:.6f}  n={len(selected)}")

    # --- 9) Output: use ORIGINAL indices from input df ---
    best_local = np.array(list(best), dtype=int)
    best_orig = orig_index[best_local]  # indices nel df originale

    calibration_df = df.loc[best_orig].copy()
    analysis_df = df.drop(best_orig).copy()

    info = {
        "best_score": float(best_score),
        "n_selected": int(len(best)),
        "excluded_batches": excluded_batches,
        "min_per_batch": min_b,
        "max_per_batch": max_b,
        "selected_counts": dict(pd.Series(batch_of[best_local]).value_counts()),
    }
    return calibration_df, analysis_df, info
data=pd.read_excel("dataset_summary_ALL.xlsx", sheet_name="ALL")

calib,analysis, info =select_calibration_set_sa_global(
    data,
    batch_col="Batch_variable_simplified",
    diagnosis_col="Diagnosis_amyloid",
    age_col="Age",
    sex_col="Sex",
    cn_label="CN-",
    min_subjects_per_batch=2,
    safe_threshold=5,
    max_batch_contribution=0.6,
    n_age_bins=20,
    target_total=100,            # if None: use only min/max per batch + optional size penalty
    lambda_batch=0.2,             # batch penalty weight (soft)
    lambda_size=0.0,              # if target_total=None, set >0 to prefer smaller sets
    iterations=5000,
    T_start=1.0,
    T_min=1e-3,
    alpha=0.995,
    random_seed=42,
    verbose_every=500,
)
# --- SAVING OUTPUT ---
output_calibration = "subjects_harmonization.xlsx"
output_analysis = "subjects_classification.xlsx"

calib.to_excel(output_calibration, index=False)
analysis.to_excel(output_analysis, index=False)

print("\n--- RESULTS ---")
print("Sex Distribution:")
print(calib['Sex'].value_counts())
print("\nAge Distribution (quantiles):")
print(calib['Age'].describe())
print("\nBatch representation:")
print(calib['Batch_variable_simplified'].value_counts())


