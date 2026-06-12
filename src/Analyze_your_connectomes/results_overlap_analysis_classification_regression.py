#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import re
import glob
import json
import warnings
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl

from scipy.stats import spearmanr, pearsonr
from matplotlib.patches import Patch
from matplotlib.colors import ListedColormap

warnings.filterwarnings("ignore")

try:
    import nibabel as nib
    from nilearn import image
    from nilearn.plotting import plot_glass_brain
    HAVE_NILEARN = True
except Exception:
    HAVE_NILEARN = False


# =============================================================================
# CONFIG
# =============================================================================

CONFIG = {
    "classification_pattern": r"/workspace/project/shared_data",
    "regression_pattern": r"/workspace/project/shared_data",

    "atlas_nii": r"/workspace/project/shared_data",
    "node_label_mapping": r"/workspace/project/shared_data",
    "label_base": "zero",

    "main_pairs": ["CN-_vs_Dementia+", "CN-_vs_MCI+", "CN+_vs_CN-"],
    "main_scores": ["ADAS13", "CDRSB"],

    "correlation_method": "spearman",
    "threshold_overlap": 0.5,

    "output_dir": "./OVERLAP_ANALYSIS_MAINPAPER",
    "dpi": 300,

    "glass_display_mode": "lzry",
    "glass_plot_abs": False,

    "color_struct": "#ef4444",
    "color_func": "#3b82f6",
    "color_both": "#10b981",
    "color_none": "#a855f7",
    "color_single_struct": "#ef4444",
    "color_single_func": "#3b82f6",
}


# =============================================================================
# NORMALIZATION
# =============================================================================

def normalize_pair_name(x):
    if pd.isna(x):
        return x
    s = str(x).strip()
    mapping = {
        "CNminus_vs_Dementiaplus": "CN-_vs_Dementia+",
        "CNminus_vs_MCIplus": "CN-_vs_MCI+",
        "CNplus_vs_CNminus": "CN+_vs_CN-",
        "CNplus_vs_Dementiaplus": "CN+_vs_Dementia+",
        "CNplus_vs_MCIplus": "CN+_vs_MCI+",
        "Dementiaplus_vs_MCIplus": "Dementia+_vs_MCI+",
    }
    return mapping.get(s, s)


def normalize_score_name(x):
    if pd.isna(x):
        return x
    s = str(x).strip()
    mapping = {
        "CDRSB": "CDRSB",
        "ADAS13": "ADAS13",
    }
    return mapping.get(s, s)


# =============================================================================
# UTILS
# =============================================================================

def ensure_dir(path):
    os.makedirs(path, exist_ok=True)
    return path


def slugify(text):
    s = str(text).strip()
    s = s.replace("+", "plus").replace("-", "minus")
    s = re.sub(r"\s+", "_", s)
    s = re.sub(r"[^A-Za-z0-9_]+", "_", s)
    s = re.sub(r"_+", "_", s).strip("_")
    return s


def detect_pair_from_classification_path(path):
    name = os.path.basename(path)
    m = re.match(r"ROI_importance_all_modalities_(.*)_for_paper\.csv$", name)
    if not m:
        return None
    return normalize_pair_name(m.group(1))


def detect_score_from_regression_path(path):
    name = os.path.basename(path)
    m = re.match(r"ROI_importance_[^_]+_(.*)\.csv$", name)
    if not m:
        return None
    return normalize_score_name(m.group(1))


def load_mapping(mapping_csv, label_base="zero", max_node=415):
    if mapping_csv is None or not os.path.exists(mapping_csv):
        start = 1 if label_base == "one" else 0
        node_to_label = {i: i for i in range(start, max_node + 1)}
        node_to_name = {i: None for i in range(start, max_node + 1)}
        return node_to_label, node_to_name

    m = pd.read_csv(mapping_csv)
    if not {"node", "atlas_label"} <= set(m.columns):
        raise ValueError("node_label_mapping must contain columns: node, atlas_label")

    node_to_label = dict(zip(m["node"].astype(int), m["atlas_label"].astype(int)))
    if "label" in m.columns:
        node_to_name = dict(zip(m["node"].astype(int), m["label"].astype(str)))
    else:
        node_to_name = {int(n): None for n in m["node"].astype(int)}
    return node_to_label, node_to_name


def build_roi_img(atlas_img, node_to_label, values_per_node):
    data = atlas_img.get_fdata().copy()
    out = np.zeros_like(data, dtype=float)
    for node, val in values_per_node.items():
        lab = node_to_label.get(int(node), None)
        if lab is None:
            continue
        out[data == lab] = val
    return image.new_img_like(atlas_img, out)


def read_importance_csv(path):
    df = pd.read_csv(path)

    rename_map = {}
    for c in df.columns:
        cl = c.lower().strip()
        if cl == "nodelabel":
            rename_map[c] = "label"
        elif cl == "node":
            rename_map[c] = "node"
        elif cl == "label":
            rename_map[c] = "label"
        elif cl == "importance_structural":
            rename_map[c] = "importance_structural"
        elif cl == "importance_functional":
            rename_map[c] = "importance_functional"
        elif cl == "importance_sf_structural":
            rename_map[c] = "importance_sf_structural"
        elif cl == "importance_sf_functional":
            rename_map[c] = "importance_sf_functional"
        elif cl == "importance_multilayer":
            rename_map[c] = "importance_multilayer"

    df = df.rename(columns=rename_map)

    required = [
        "node",
        "importance_structural",
        "importance_functional",
        "importance_sf_structural",
        "importance_sf_functional",
        "importance_multilayer",
    ]
    missing = [c for c in required if c not in df.columns]
    if missing:
        raise ValueError(f"{os.path.basename(path)} missing columns: {missing}")

    keep = ["node"] + [c for c in ["label"] if c in df.columns] + required[1:]
    df = df[keep].copy()
    df["node"] = df["node"].astype(int)

    for c in required[1:]:
        df[c] = pd.to_numeric(df[c], errors="coerce")

    return df.sort_values("node").reset_index(drop=True)


def load_classification_tables(pattern):
    files = sorted(glob.glob(pattern))
    if not files:
        raise ValueError(f"No classification files found: {pattern}")

    tables = {}
    for f in files:
        raw_df = pd.read_csv(f, nrows=5)
        df = read_importance_csv(f)

        pair = None
        if "pair" in raw_df.columns and len(raw_df) > 0:
            vals = raw_df["pair"].dropna().astype(str)
            if len(vals):
                pair = normalize_pair_name(vals.iloc[0])

        if pair is None:
            pair = detect_pair_from_classification_path(f)

        if pair is None:
            print(f"WARNING: could not detect pair for {os.path.basename(f)}")
            continue

        tables[pair] = df
        print(f"[CLASS] {os.path.basename(f)} -> {pair} ({len(df)} rows)")

    return tables


def load_regression_tables(pattern):
    files = sorted(glob.glob(pattern))
    if not files:
        raise ValueError(f"No regression files found: {pattern}")

    tables = {}
    for f in files:
        raw_df = pd.read_csv(f, nrows=5)
        df = read_importance_csv(f)

        score = None
        if "score" in raw_df.columns and len(raw_df) > 0:
            vals = raw_df["score"].dropna().astype(str)
            if len(vals):
                score = normalize_score_name(vals.iloc[0])

        if score is None:
            score = detect_score_from_regression_path(f)

        if score is None:
            print(f"WARNING: could not detect score for {os.path.basename(f)}")
            continue

        tables[score] = df
        print(f"[REGR]  {os.path.basename(f)} -> {score} ({len(df)} rows)")

    return tables


def get_modality_profile(df, modality):
    if modality == "structural":
        return df[["node", "importance_structural"]].rename(
            columns={"importance_structural": "score"}
        ).copy()

    if modality == "functional":
        return df[["node", "importance_functional"]].rename(
            columns={"importance_functional": "score"}
        ).copy()

    if modality == "multilayer":
        return df[["node", "importance_multilayer"]].rename(
            columns={"importance_multilayer": "score"}
        ).copy()

    if modality == "structural+functional":
        return df[["node", "importance_sf_structural", "importance_sf_functional"]].copy()

    raise ValueError(f"Unknown modality: {modality}")


# =============================================================================
# CORRELATIONS
# =============================================================================

def compute_correlations(class_tables, reg_tables, cfg):
    rows = []

    for pair in cfg["main_pairs"]:
        if pair not in class_tables:
            print(f"[MISS CLASS] {pair}")
            continue
        cdf_all = class_tables[pair]

        for score in cfg["main_scores"]:
            if score not in reg_tables:
                print(f"[MISS REGR] {score}")
                continue
            rdf_all = reg_tables[score]

            for mod in ["structural", "functional", "structural+functional", "multilayer"]:
                cdf = get_modality_profile(cdf_all, mod)
                rdf = get_modality_profile(rdf_all, mod)

                if mod == "structural+functional":
                    merged = pd.merge(cdf, rdf, on="node", suffixes=("_class", "_regr"))
                    x = np.concatenate([
                        merged["importance_sf_structural_class"].values,
                        merged["importance_sf_functional_class"].values
                    ])
                    y = np.concatenate([
                        merged["importance_sf_structural_regr"].values,
                        merged["importance_sf_functional_regr"].values
                    ])
                else:
                    merged = pd.merge(cdf, rdf, on="node", suffixes=("_class", "_regr"))
                    x = merged["score_class"].values
                    y = merged["score_regr"].values

                valid = np.isfinite(x) & np.isfinite(y)
                x = x[valid]
                y = y[valid]

                if len(x) < 3:
                    corr, pval = np.nan, np.nan
                else:
                    if cfg["correlation_method"] == "pearson":
                        corr, pval = pearsonr(x, y)
                    else:
                        corr, pval = spearmanr(x, y)

                rows.append({
                    "pair": pair,
                    "score": score,
                    "modality": mod,
                    "correlation": corr,
                    "pvalue": pval,
                    "n_values": int(len(x)),
                    "n_rois": int(len(merged)),
                })

                if pd.notna(corr):
                    print(f"[CORR] {pair} | {score} | {mod} -> n={len(merged)} corr={corr:.4f}")
                else:
                    print(f"[CORR] {pair} | {score} | {mod} -> n={len(merged)} corr=NA")

    return pd.DataFrame(rows)


# =============================================================================
# OVERLAP LOGIC
# =============================================================================

def get_sets_for_four_panel_plot(class_df, regr_df, thr):
    out = {}

    merged_s = pd.merge(
        class_df[["node", "importance_structural"]],
        regr_df[["node", "importance_structural"]],
        on="node", suffixes=("_class", "_regr")
    )
    out["structural"] = set(
        merged_s.loc[
            (merged_s["importance_structural_class"] >= thr) &
            (merged_s["importance_structural_regr"] >= thr),
            "node"
        ].astype(int).tolist()
    )

    merged_f = pd.merge(
        class_df[["node", "importance_functional"]],
        regr_df[["node", "importance_functional"]],
        on="node", suffixes=("_class", "_regr")
    )
    out["functional"] = set(
        merged_f.loc[
            (merged_f["importance_functional_class"] >= thr) &
            (merged_f["importance_functional_regr"] >= thr),
            "node"
        ].astype(int).tolist()
    )

    merged_sf = pd.merge(
        class_df[["node", "importance_sf_structural", "importance_sf_functional"]],
        regr_df[["node", "importance_sf_structural", "importance_sf_functional"]],
        on="node", suffixes=("_class", "_regr")
    )

    sf_struct_sel = set(
        merged_sf.loc[
            (merged_sf["importance_sf_structural_class"] >= thr) &
            (merged_sf["importance_sf_structural_regr"] >= thr),
            "node"
        ].astype(int).tolist()
    )

    sf_func_sel = set(
        merged_sf.loc[
            (merged_sf["importance_sf_functional_class"] >= thr) &
            (merged_sf["importance_sf_functional_regr"] >= thr),
            "node"
        ].astype(int).tolist()
    )

    sf_both = sf_struct_sel & sf_func_sel
    sf_struct_only = sf_struct_sel - sf_both
    sf_func_only = sf_func_sel - sf_both

    out["sf_struct_only"] = sf_struct_only
    out["sf_func_only"] = sf_func_only
    out["sf_both"] = sf_both

    merged_ml = pd.merge(
        class_df[["node", "importance_multilayer"]],
        regr_df[["node", "importance_multilayer"]],
        on="node", suffixes=("_class", "_regr")
    )
    ml_sel = set(
        merged_ml.loc[
            (merged_ml["importance_multilayer_class"] >= thr) &
            (merged_ml["importance_multilayer_regr"] >= thr),
            "node"
        ].astype(int).tolist()
    )

    ml_fs = ml_sel & sf_struct_sel & sf_func_sel
    ml_f = (ml_sel & sf_func_sel) - sf_struct_sel
    ml_s = (ml_sel & sf_struct_sel) - sf_func_sel
    ml_none = ml_sel - (ml_fs | ml_f | ml_s)

    out["ml_f"] = ml_f
    out["ml_s"] = ml_s
    out["ml_fs"] = ml_fs
    out["ml_none"] = ml_none

    return out


# =============================================================================
# BRAIN PLOT
# =============================================================================

def _plot_four_panel_brain(atlas_img, node_to_label, sets_dict, title, out_png, cfg):
    saved = []

    if not HAVE_NILEARN or atlas_img is None or node_to_label is None:
        return saved

    structural_nodes = sets_dict["structural"]
    functional_nodes = sets_dict["functional"]

    sf_struct = sets_dict["sf_struct_only"]
    sf_func = sets_dict["sf_func_only"]
    sf_both = sets_dict["sf_both"]

    ml_f = sets_dict["ml_f"]
    ml_s = sets_dict["ml_s"]
    ml_fs = sets_dict["ml_fs"]
    ml_none = sets_dict["ml_none"]

    img_struct = build_roi_img(atlas_img, node_to_label, {n: 1.0 for n in structural_nodes}) if structural_nodes else None
    img_func = build_roi_img(atlas_img, node_to_label, {n: 1.0 for n in functional_nodes}) if functional_nodes else None

    sf_values = {}
    for n in sf_func:
        sf_values[n] = 1.0
    for n in sf_struct:
        sf_values[n] = 2.0
    for n in sf_both:
        sf_values[n] = 3.0
    img_sf = build_roi_img(atlas_img, node_to_label, sf_values) if sf_values else None

    ml_values = {}
    for n in ml_f:
        ml_values[n] = 1.0
    for n in ml_s:
        ml_values[n] = 2.0
    for n in ml_fs:
        ml_values[n] = 3.0
    for n in ml_none:
        ml_values[n] = 4.0
    img_ml = build_roi_img(atlas_img, node_to_label, ml_values) if ml_values else None

    fig = plt.figure(figsize=(10, 10))
    gs = fig.add_gridspec(5, 1, height_ratios=[1, 1, 1, 1, 0.18], hspace=0.04)

    ax_ml = fig.add_subplot(gs[3, 0])
    ax_sf = fig.add_subplot(gs[2, 0])
    ax_sc = fig.add_subplot(gs[0, 0])
    ax_fc = fig.add_subplot(gs[1, 0])
    ax_leg = fig.add_subplot(gs[4, 0])
    ax_leg.axis("off")

    if img_ml is not None:
        cmap_ml = ListedColormap([
            cfg["color_func"],
            cfg["color_struct"],
            cfg["color_both"],
            cfg["color_none"],
        ])
        plot_glass_brain(
            img_ml,
            display_mode=cfg["glass_display_mode"],
            colorbar=False,
            plot_abs=cfg["glass_plot_abs"],
            vmax=4.5,
            cmap=cmap_ml,
            axes=ax_ml
            
        )
    else:
        ax_ml.set_title("Multilayer (none)")

    if img_sf is not None:
        cmap_sf = ListedColormap([
            cfg["color_func"],
            cfg["color_struct"],
            cfg["color_both"],
        ])
        plot_glass_brain(
            img_sf,
            display_mode=cfg["glass_display_mode"],
            colorbar=False,
            plot_abs=cfg["glass_plot_abs"],
            vmax=3.5,
            cmap=cmap_sf,
            axes=ax_sf
            
        )
    else:
        ax_sf.set_title("Structural+Functional (none)")

    if img_struct is not None:
        cmap_struct = ListedColormap([
            (1, 1, 1, 0),
            mpl.colors.to_rgba(cfg["color_single_struct"])
        ])
        plot_glass_brain(
            img_struct,
            display_mode=cfg["glass_display_mode"],
            colorbar=False,
            plot_abs=cfg["glass_plot_abs"],
            vmax=1.0,
            cmap=cmap_struct,
            axes=ax_sc
            
        )
    else:
        ax_sc.set_title("Structural (none)")

    if img_func is not None:
        cmap_func = ListedColormap([
            (1, 1, 1, 0),
            mpl.colors.to_rgba(cfg["color_single_func"])
        ])
        plot_glass_brain(
            img_func,
            display_mode=cfg["glass_display_mode"],
            colorbar=False,
            plot_abs=cfg["glass_plot_abs"],
            vmax=1.0,
            cmap=cmap_func,
            axes=ax_fc
            
        )
    else:
        ax_fc.set_title("Functional (none)")

    handles = [
        Patch(facecolor=cfg["color_func"], label="FC-only"),
        Patch(facecolor=cfg["color_struct"], label="SC-only"),
        Patch(facecolor=cfg["color_both"], label="SC+FC"),
        Patch(facecolor=cfg["color_none"], label="ML-only"),
    ]
    ax_leg.legend(handles=handles, loc="center", frameon=False, ncol=4)

    fig.suptitle(title, fontsize=13, y=0.995)
    fig.savefig(out_png, dpi=cfg["dpi"], bbox_inches="tight", pad_inches=0.02)
    plt.close(fig)
    saved.append(out_png)

    return saved


# =============================================================================
# HEATMAP
# =============================================================================

def plot_correlation_heatmaps(corr_df, outdir, cfg):
    if corr_df is None or corr_df.empty:
        print("[HEATMAP] corr_df empty")
        return []

    saved = []

    pair_display = {
        "CN-_vs_Dementia+": "CN Aβ− vs AD Aβ+",
        "CN-_vs_MCI+": "CN Aβ− vs MCI Aβ+",
        "CN+_vs_CN-": "CN Aβ+ vs CN Aβ−",
    }

    modality_order = [ "structural", "functional","structural+functional","multilayer" ]
    modality_display = {
        "multilayer": "Multilayer",
        "structural+functional": "SC+FC",
        "structural": "SC",
        "functional": "FC",
    }

    score_order = ["ADAS13", "CDRSB"]

    df = corr_df.copy()
    df = df[df["pair"].isin(cfg["main_pairs"]) & df["score"].isin(cfg["main_scores"])].copy()

    print(f"[HEATMAP] input rows after filter: {len(df)}")

    if df.empty:
        return saved

    pivot_rows = []
    for pair in cfg["main_pairs"]:
        sub_pair = df[df["pair"] == pair].copy()
        row = {"pair": pair_display.get(pair, pair)}
        for score in score_order:
            for mod in modality_order:
                vv = sub_pair.loc[
                    (sub_pair["score"] == score) &
                    (sub_pair["modality"] == mod),
                    "correlation"
                ]
                pv = sub_pair.loc[
                    (sub_pair["score"] == score) &
                    (sub_pair["modality"] == mod),
                    "pvalue"
                ]
                row[(score, modality_display[mod])] = float(vv.iloc[0]) if len(vv) else np.nan
                row[(f"{score}_p", modality_display[mod])] = float(pv.iloc[0]) if len(pv) else np.nan
        pivot_rows.append(row)

    heat_df = pd.DataFrame(pivot_rows).set_index("pair")

    col_tuples = [(s, modality_display[m]) for s in score_order for m in modality_order]
    val_df = pd.DataFrame(index=heat_df.index)
    p_df = pd.DataFrame(index=heat_df.index)

    for c in col_tuples:
        val_df[c] = heat_df[c]
        p_df[c] = heat_df[(f"{c[0]}_p", c[1])]

    fig, ax = plt.subplots(figsize=(10, 4.8))
    vals = val_df.to_numpy(dtype=float)
    im = ax.imshow(vals, cmap="RdBu_r", vmin=-0.4, vmax=0.4, aspect="auto")

    ax.set_yticks(np.arange(len(val_df.index)))
    ax.set_yticklabels(val_df.index, fontsize=10)
    ax.set_xticks(np.arange(len(val_df.columns)))
    ax.set_xticklabels([c[1] for c in val_df.columns], fontsize=10)

    n_mod = len(modality_order)
    ax.axvline(n_mod - 0.5, color="black", linewidth=1.2)

    ax.set_xticks(np.arange(-.5, vals.shape[1], 1), minor=True)
    ax.set_yticks(np.arange(-.5, vals.shape[0], 1), minor=True)
    ax.grid(which="minor", color="white", linestyle="-", linewidth=1.5)
    ax.tick_params(which="minor", bottom=False, left=False)

    top_y = -1.15
    ax.text((0 + n_mod - 1) / 2, top_y, "ADAS13", ha="center", va="bottom", fontsize=11, fontweight="bold")
    ax.text((n_mod + 2 * n_mod - 1) / 2, top_y, "CDRSB", ha="center", va="bottom", fontsize=11, fontweight="bold")

    for i in range(vals.shape[0]):
        for j in range(vals.shape[1]):
            v = vals[i, j]
            p = p_df.iloc[i, j]
            if np.isnan(v):
                txt = "NA"
            else:
                star = ""
                if pd.notna(p):
                    if p < 0.001:
                        star = "***"
                    elif p < 0.01:
                        star = "**"
                    elif p < 0.05:
                        star = "*"
                txt = f"{v:.2f}{star}"
            text_color = "white" if np.isfinite(v) and abs(v) >= 0.20 else "black"
            ax.text(j, i, txt, ha="center", va="center", fontsize=9, color=text_color)

    cbar = fig.colorbar(im, ax=ax, fraction=0.035, pad=0.02)
    cbar.set_label("Spearman correlation" if cfg["correlation_method"] == "spearman" else "Pearson correlation")

    ax.set_title("Classification–Regression convergence", fontsize=13, fontweight="bold", pad=28)
    fig.tight_layout()

    out_png = os.path.join(outdir, "heatmap_convergence_mainpaper.png")
    fig.savefig(out_png, dpi=cfg["dpi"], bbox_inches="tight")
    plt.close(fig)
    saved.append(out_png)

    for score in score_order:
        sub_cols = [(score, modality_display[m]) for m in modality_order]
        sub_vals = val_df[sub_cols].to_numpy(dtype=float)
        sub_p = p_df[sub_cols]

        fig, ax = plt.subplots(figsize=(5.8, 4.6))
        im = ax.imshow(sub_vals, cmap="RdBu_r", vmin=-0.4, vmax=0.4, aspect="auto")

        ax.set_yticks(np.arange(len(val_df.index)))
        ax.set_yticklabels(val_df.index, fontsize=10)
        ax.set_xticks(np.arange(len(sub_cols)))
        ax.set_xticklabels([c[1] for c in sub_cols], fontsize=10)

        ax.set_xticks(np.arange(-.5, sub_vals.shape[1], 1), minor=True)
        ax.set_yticks(np.arange(-.5, sub_vals.shape[0], 1), minor=True)
        ax.grid(which="minor", color="white", linestyle="-", linewidth=1.5)
        ax.tick_params(which="minor", bottom=False, left=False)

        for i in range(sub_vals.shape[0]):
            for j in range(sub_vals.shape[1]):
                v = sub_vals[i, j]
                p = sub_p.iloc[i, j]
                if np.isnan(v):
                    txt = "NA"
                else:
                    star = ""
                    if pd.notna(p):
                        if p < 0.001:
                            star = "***"
                        elif p < 0.01:
                            star = "**"
                        elif p < 0.05:
                            star = "*"
                    txt = f"{v:.2f}{star}"
                text_color = "white" if np.isfinite(v) and abs(v) >= 0.20 else "black"
                ax.text(j, i, txt, ha="center", va="center", fontsize=9, color=text_color)

        cbar = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.03)
        cbar.set_label("Spearman correlation" if cfg["correlation_method"] == "spearman" else "Pearson correlation")

        ax.set_title(f"Convergence heatmap - {score}", fontsize=12, fontweight="bold")
        fig.tight_layout()

        out_png = os.path.join(outdir, f"heatmap_convergence_{score}.png")
        fig.savefig(out_png, dpi=cfg["dpi"], bbox_inches="tight")
        plt.close(fig)
        saved.append(out_png)

    return saved


# =============================================================================
# MAIN
# =============================================================================

def main(cfg):
    outdir = ensure_dir(cfg["output_dir"])

    print("=" * 80)
    print("OVERLAP ANALYSIS - MAIN PAPER")
    print("=" * 80)

    print("\nLoading classification tables...")
    class_tables = load_classification_tables(cfg["classification_pattern"])
    print(f"Loaded {len(class_tables)} classification tables")

    print("\nLoading regression tables...")
    reg_tables = load_regression_tables(cfg["regression_pattern"])
    print(f"Loaded {len(reg_tables)} regression tables")

    print("\nLoaded classification keys:")
    for k in sorted(class_tables.keys()):
        print(" -", repr(k))

    print("\nLoaded regression keys:")
    for k in sorted(reg_tables.keys()):
        print(" -", repr(k))

    print("\nMain paper classification tasks:")
    for p in cfg["main_pairs"]:
        print(" -", p)

    print("\nRegression scores:")
    for s in cfg["main_scores"]:
        print(" -", s)

    missing_pairs = [p for p in cfg["main_pairs"] if p not in class_tables]
    missing_scores = [s for s in cfg["main_scores"] if s not in reg_tables]

    if missing_pairs:
        print("\nDEBUG WARNING - missing classification pairs:")
        for p in missing_pairs:
            print(" -", p)

    if missing_scores:
        print("\nDEBUG WARNING - missing regression scores:")
        for s in missing_scores:
            print(" -", s)

    atlas_img = None
    node_to_label = None
    node_to_name = None

    if HAVE_NILEARN and cfg.get("atlas_nii") and os.path.exists(cfg["atlas_nii"]):
        atlas_img = nib.load(cfg["atlas_nii"])
        all_nodes = []
        for d in list(class_tables.values()) + list(reg_tables.values()):
            all_nodes.extend(d["node"].tolist())
        max_node = max(all_nodes) if all_nodes else 415
        node_to_label, node_to_name = load_mapping(
            cfg["node_label_mapping"],
            cfg["label_base"],
            max_node=max_node
        )
        print("\nAtlas loaded")
    else:
        print("\nAtlas not available or nilearn missing: brain plots skipped")

    print("\nComputing correlations...")
    corr_df = compute_correlations(class_tables, reg_tables, cfg)
    corr_csv = os.path.join(outdir, "correlation_results_mainpaper.csv")
    corr_xlsx = os.path.join(outdir, "correlation_results_mainpaper.xlsx")
    corr_df.to_csv(corr_csv, index=False)
    corr_df.to_excel(corr_xlsx, index=False)
    print("Saved:", corr_csv)

    print("\nGenerating heatmaps...")
    heatmap_files = plot_correlation_heatmaps(corr_df, outdir, cfg)
    print(f"Heatmaps created: {len(heatmap_files)}")

    print("\nGenerating four-panel brain plots...")
    all_rows = []
    saved_figures = []

    for pair in cfg["main_pairs"]:
        if pair not in class_tables:
            print(f"WARNING: classification pair not found: {pair}")
            continue

        class_df = class_tables[pair]

        for score in cfg["main_scores"]:
            if score not in reg_tables:
                print(f"WARNING: regression score not found: {score}")
                continue

            regr_df = reg_tables[score]

            sets_dict = get_sets_for_four_panel_plot(
                class_df=class_df,
                regr_df=regr_df,
                thr=cfg["threshold_overlap"]
            )

            set_sizes = {k: len(v) for k, v in sets_dict.items()}
            print(f"[OVERLAP] {pair} | {score} -> {set_sizes}")

            for set_name, nodes in sets_dict.items():
                for n in sorted(nodes):
                    all_rows.append({
                        "pair": pair,
                        "score": score,
                        "set": set_name,
                        "node": int(n),
                        "label": node_to_name.get(int(n)) if node_to_name else None,
                        "threshold": cfg["threshold_overlap"],
                    })

            if atlas_img is not None and node_to_label is not None:
                out_png = os.path.join(
                    outdir,
                    f"OVERLAP_MAINPAPER_{slugify(pair)}_{slugify(score)}_thr{str(cfg['threshold_overlap']).replace('.', 'p')}.png"
                )
                saved = _plot_four_panel_brain(
                    atlas_img=atlas_img,
                    node_to_label=node_to_label,
                    sets_dict=sets_dict,
                    title=f"{pair} | {score} | thr={cfg['threshold_overlap']}",
                    out_png=out_png,
                    cfg=cfg
                )
                saved_figures.extend(saved)

    overlap_df = pd.DataFrame(all_rows)
    overlap_csv = os.path.join(outdir, "overlap_nodes_mainpaper.csv")
    overlap_xlsx = os.path.join(outdir, "overlap_nodes_mainpaper.xlsx")
    overlap_df.to_csv(overlap_csv, index=False)
    overlap_df.to_excel(overlap_xlsx, index=False)

    manifest = {
        "config": cfg,
        "figures": saved_figures + heatmap_files,
        "brain_figures": saved_figures,
        "heatmaps": heatmap_files,
        "n_figures": len(saved_figures) + len(heatmap_files),
        "n_brain_figures": len(saved_figures),
        "n_heatmaps": len(heatmap_files),
        "n_overlap_rows": int(len(overlap_df)),
        "n_correlations": int(len(corr_df)),
        "correlation_csv": corr_csv,
        "overlap_csv": overlap_csv,
    }

    with open(os.path.join(outdir, "manifest_mainpaper.json"), "w", encoding="utf-8") as f:
        json.dump(manifest, f, indent=2)

    print("\nDone.")
    print("Output dir:", outdir)
    print("Correlation rows:", len(corr_df))
    print("Overlap rows:", len(overlap_df))
    print("Brain figures:", len(saved_figures))
    print("Heatmaps:", len(heatmap_files))

    return corr_df, overlap_df


if __name__ == "__main__":
    main(CONFIG)