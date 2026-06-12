#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Plot ROI importance maps for paper.

For each classification task (pair), create a figure with 4 subplots:
1) structural (modalità structural)
2) functional (modalità functional)
3) structural+functional (ROI insieme con colori distinti per structural / functional / overlap)
4) multilayer (ROI ML colorate in base all'overlap con structural e functional)

Also export, for each pair, a complete CSV with unified importance values
for each modality and for the structural / functional subcomponents derived from the
structural+functional modality.

Filename suffix: _for_paper
"""

import os
import re
import json
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

try:
    import nibabel as nib
    from nilearn import image
    from nilearn.plotting import plot_stat_map, plot_glass_brain
    HAVE_NILEARN = True
except Exception:
    HAVE_NILEARN = False

from matplotlib.colors import ListedColormap, BoundaryNorm
from matplotlib.patches import Patch

CONFIG = {
    "fi_csv": "/workspace/project/shared_data",
    "harmonization": "combat",
    "pairs": ["CN+_vs_Dementia+", "CN+_vs_MCI+", "Dementia+_vs_MCI+", "CN-_vs_Dementia+", "CN-_vs_MCI+", "CN+_vs_CN-"],
    "score_threshold": 1,
    "normalization": "zscore",
    "atlas_nii": "/workspace/project/shared_data",
    "node_label_mapping": "/workspace/project/shared_data",
    "label_base": "zero",
    "output_dir": None,
    "dpi": 400,
    "transparent_brain": "glass",
    "glass_display_mode": "lzry",
    "glass_plot_abs": False,
    "bg_img": None,
    "color_func": "#3b82f6",
    "color_struct": "#ef4444",
    "color_both": "#10b981",
    "color_ml_only": "#a855f7",
}

RE_NODE = re.compile(r"(?:^|_)(\d+)$")


def _slug_pair(pair: str) -> str:
    s = str(pair).strip()
    s = s.replace('+', 'plus').replace('-', 'minus')
    s = re.sub(r'\s*vs\s*', '_vs_', s, flags=re.IGNORECASE)
    s = re.sub(r'[^A-Za-z0-9]+', '_', s)
    s = re.sub(r'_+', '_', s).strip('_')
    return s


def ensure_outdir(path):
    os.makedirs(path, exist_ok=True)
    return path


def _extract_node(feat: str):
    m = RE_NODE.search(str(feat))
    if not m:
        return None
    try:
        return int(m.group(1))
    except Exception:
        return None


def _split_sf(feat: str):
    node = _extract_node(feat)
    if node is None:
        return None, None, None
    base = RE_NODE.sub("", str(feat))
    if base.startswith("func_"):
        return "func", base[len("func_"):], node
    if base.startswith("struct_"):
        return "struct", base[len("struct_"):], node
    return None, base, node


def _normalize(values: dict, how: str) -> dict:
    if how == "none":
        return dict(values)
    arr = np.array(list(values.values()), dtype=float)
    if arr.size == 0:
        return {}
    if how == "zscore":
        mu, sd = float(np.mean(arr)), float(np.std(arr))
        if sd == 0:
            return {n: 0.0 for n in values}
        return {n: (v - mu) / sd for n, v in values.items()}
    if how == "minmax":
        mn, mx = float(np.min(arr)), float(np.max(arr))
        if mx == mn:
            return {n: 0.0 for n in values}
        return {n: (v - mn) / (mx - mn) for n, v in values.items()}
    raise ValueError(f"Unsupported normalization: {how}")


def _fuse_metrics(metric_to_values: dict, normalization: str = "zscore") -> dict:
    if not metric_to_values:
        return {}
    normed = {m: _normalize(d, normalization) for m, d in metric_to_values.items()}
    nodes = set().union(*[set(d.keys()) for d in normed.values()]) if normed else set()
    fused = {}
    for n in nodes:
        vals = [d[n] for d in normed.values() if n in d]
        if vals:
            fused[n] = float(np.mean(vals))
    return fused


def _aggregate_by_mode(df: pd.DataFrame, mode: str, harm: str, pair: str) -> dict:
    sub = df[(df["mode"] == mode) & (df["harmonization"] == harm) & (df["pair"] == pair)]
    if sub.empty:
        return {}
    agg = sub.groupby("feature")["importance_mean"].mean().reset_index(name="imp")
    grouped = {}
    if mode == "structural+functional":
        for _, r in agg.iterrows():
            feat, val = r["feature"], float(r["imp"])
            side, base, node = _split_sf(feat)
            if node is None:
                continue
            if side not in {"func", "struct"}:
                side = "struct"
            metric_dict = grouped.setdefault(base, {"func": {}, "struct": {}})
            metric_dict[side][node] = max(val, metric_dict[side].get(node, -np.inf))
    else:
        for _, r in agg.iterrows():
            feat, val = r["feature"], float(r["imp"])
            node = _extract_node(feat)
            if node is None:
                continue
            base = RE_NODE.sub("", str(feat))
            d = grouped.setdefault(base, {})
            d[node] = max(val, d.get(node, -np.inf))
    return grouped


def _unified_from_sf(grouped_sf: dict, normalization: str):
    func_metrics, struct_metrics = {}, {}
    for metric, sides in grouped_sf.items():
        if sides.get("func"):
            func_metrics[metric] = sides["func"]
        if sides.get("struct"):
            struct_metrics[metric] = sides["struct"]
    return _fuse_metrics(func_metrics, normalization), _fuse_metrics(struct_metrics, normalization)


def load_mapping(mapping_csv: str, label_base: str, max_node: int):
    if mapping_csv is None:
        start = 1 if label_base == "one" else 0
        node_to_label = {i: i for i in range(start, max_node + 1)}
        node_to_name = {i: None for i in range(start, max_node + 1)}
        return node_to_label, node_to_name
    m = pd.read_csv(mapping_csv)
    if not {"node", "atlas_label"} <= set(m.columns):
        raise ValueError("node_label_mapping must contain node, atlas_label")
    node_to_label = dict(zip(m["node"].astype(int), m["atlas_label"].astype(int)))
    node_to_name = dict(zip(m["node"].astype(int), m["label"].astype(str))) if "label" in m.columns else {int(n): None for n in m["node"].astype(int)}
    return node_to_label, node_to_name


def build_roi_img(atlas_img, node_to_label: dict, values_per_node: dict):
    data = atlas_img.get_fdata()
    out = np.zeros_like(data, dtype=float)
    for node, val in values_per_node.items():
        lab = node_to_label.get(node, None)
        if lab is None:
            continue
        out[data == lab] = val
    return image.new_img_like(atlas_img, out)


def _mono_cmap(hex_color: str):
    import matplotlib.colors as mcolors
    base = mcolors.to_rgb(hex_color)
    return ListedColormap([(1, 1, 1, 0), base])


def _plot_panel_glass(ax, img, cmap, vmax, cfg):
    if img is None:
        ax.axis("off")
        return
    plot_glass_brain(
        img,
        display_mode=cfg.get("glass_display_mode", "lyr"),
        colorbar=False,
        plot_abs=cfg.get("glass_plot_abs", False),
        vmax=vmax,
        cmap=cmap,
        alpha=1,
        axes=ax,
    )


def _plot_panel_stat(ax, img, cmap, vmin, vmax, cfg):
    if img is None:
        ax.axis("off")
        return
    plot_stat_map(
        img,
        display_mode="ortho",
        cut_coords=4,
        axes=ax,
        colorbar=False,
        annotate=False,
        draw_cross=False,
        bg_img=None,
        black_bg=(cfg.get("transparent_brain") == "overlay"),
        dim=0,
        vmin=vmin,
        vmax=vmax,
        cmap=cmap,
    )


def _plot_four_panel(atlas_img, node_to_label, sets_dict, title, out_png, cfg):
    struct_nodes = sets_dict["structural_only"]
    func_nodes = sets_dict["functional_only"]
    sf_struct = sets_dict["sf_struct"]
    sf_func = sets_dict["sf_func"]
    sf_both = sets_dict["sf_both"]
    ml_f = sets_dict["ml_f"]
    ml_s = sets_dict["ml_s"]
    ml_fs = sets_dict["ml_fs"]
    ml_none = sets_dict["ml_none"]

    if HAVE_NILEARN and atlas_img is not None and node_to_label is not None:
        struct_img = build_roi_img(atlas_img, node_to_label, {n: 1.0 for n in struct_nodes}) if struct_nodes else None
        func_img = build_roi_img(atlas_img, node_to_label, {n: 1.0 for n in func_nodes}) if func_nodes else None

        sf_values = {}
        for n in sf_func:
            sf_values[n] = 1.0
        for n in sf_struct:
            sf_values[n] = 2.0
        for n in sf_both:
            sf_values[n] = 3.0
        sf_img = build_roi_img(atlas_img, node_to_label, sf_values) if sf_values else None

        ml_values = {}
        for n in ml_f:
            ml_values[n] = 1.0
        for n in ml_s:
            ml_values[n] = 2.0
        for n in ml_fs:
            ml_values[n] = 3.0
        for n in ml_none:
            ml_values[n] = 4.0
        ml_img = build_roi_img(atlas_img, node_to_label, ml_values) if ml_values else None

        cmap_tri = ListedColormap([cfg["color_func"], cfg["color_struct"], cfg["color_both"]])
        cmap_quad = ListedColormap([cfg["color_func"], cfg["color_struct"], cfg["color_both"], cfg["color_ml_only"]])

        if cfg.get("transparent_brain") == "glass":
            fig = plt.figure(figsize=(10, 8.5))
            gs = fig.add_gridspec(
                3, 2,
                height_ratios=[1, 1, 0.16],
                width_ratios=[1, 1],
                left=0.03, right=0.97, top=0.98, bottom=0.03,
                wspace=0.04, hspace=0.06
            )

            ax1 = fig.add_subplot(gs[0, 0])
            ax2 = fig.add_subplot(gs[0, 1])
            ax3 = fig.add_subplot(gs[1, 0])
            ax4 = fig.add_subplot(gs[1, 1])
            axL = fig.add_subplot(gs[2, :])
            axL.axis("off")

            _plot_panel_glass(ax1, struct_img, _mono_cmap(cfg["color_struct"]), 1.0, cfg)
            _plot_panel_glass(ax2, func_img, _mono_cmap(cfg["color_func"]), 1.0, cfg)
            _plot_panel_glass(ax3, sf_img, cmap_tri, 3.5, cfg)
            _plot_panel_glass(ax4, ml_img, cmap_quad, 4.5, cfg)

            legend_handles = [
                Patch(facecolor=cfg["color_struct"], label="Structural"),
                Patch(facecolor=cfg["color_func"], label="Functional"),
                Patch(facecolor=cfg["color_both"], label="Both"),
                Patch(facecolor=cfg["color_ml_only"], label="MUL-only"),
            ]
            axL.legend(handles=legend_handles, loc="center", frameon=False, ncol=4)
            fig.savefig(out_png, dpi=cfg["dpi"], bbox_inches="tight", pad_inches=0.01)
            plt.close(fig)
            return [out_png]
        else:
            fig = plt.figure(figsize=(10, 10))
            gs = fig.add_gridspec(5, 1, height_ratios=[1, 1, 1, 1, 0.16], left=0.03, right=0.97, top=0.97, bottom=0.03, hspace=0.08)
            ax1 = fig.add_subplot(gs[0, 0])
            ax2 = fig.add_subplot(gs[1, 0])
            ax3 = fig.add_subplot(gs[2, 0])
            ax4 = fig.add_subplot(gs[3, 0])
            axL = fig.add_subplot(gs[4, 0])
            axL.axis("off")

            _plot_panel_stat(ax1, struct_img, _mono_cmap(cfg["color_struct"]), 0.0, 1.0, cfg)
            _plot_panel_stat(ax2, func_img, _mono_cmap(cfg["color_func"]), 0.0, 1.0, cfg)
            _plot_panel_stat(ax3, sf_img, cmap_tri, 0.5, 3.5, cfg)
            _plot_panel_stat(ax4, ml_img, cmap_quad, 0.5, 4.5, cfg)

            legend_handles = [
                Patch(facecolor=cfg["color_struct"], label="Structural"),
                Patch(facecolor=cfg["color_func"], label="Functional"),
                Patch(facecolor=cfg["color_both"], label="Both"),
                Patch(facecolor=cfg["color_ml_only"], label="ML-only"),
            ]
            axL.legend(handles=legend_handles, loc="center", frameon=False, ncol=4)
            fig.savefig(out_png, dpi=cfg["dpi"], bbox_inches="tight", pad_inches=0.01)
            plt.close(fig)
            return [out_png]

    def _bar(ax, nodes, color, ttl):
        if nodes:
            xs = sorted(nodes)
            ax.bar(xs, [1] * len(xs), color=color)
            ax.set_title(ttl)
            ax.set_xlabel("node")
            ax.set_yticks([])
        else:
            ax.set_title(ttl + " (none)")
            ax.axis("off")

    fig = plt.figure(figsize=(14, 8))
    gs = fig.add_gridspec(2, 2, wspace=0.25, hspace=0.35)
    ax1 = fig.add_subplot(gs[0, 0])
    ax2 = fig.add_subplot(gs[0, 1])
    ax3 = fig.add_subplot(gs[1, 0])
    ax4 = fig.add_subplot(gs[1, 1])

    _bar(ax1, struct_nodes, cfg["color_struct"], "structural")
    _bar(ax2, func_nodes, cfg["color_func"], "functional")

    for lab, color, nodes in [("functional", cfg["color_func"], sf_func), ("structural", cfg["color_struct"], sf_struct), ("both", cfg["color_both"], sf_both)]:
        if nodes:
            ax3.bar(sorted(nodes), [1] * len(nodes), color=color, label=lab)
    ax3.legend(frameon=False)
    ax3.set_title("structural+functional")
    ax3.set_xlabel("node")
    ax3.set_yticks([])

    for lab, color, nodes in [("F", cfg["color_func"], ml_f), ("S", cfg["color_struct"], ml_s), ("FS", cfg["color_both"], ml_fs), ("ML-only", cfg["color_ml_only"], ml_none)]:
        if nodes:
            ax4.bar(sorted(nodes), [1] * len(nodes), color=color, label=lab)
    ax4.legend(frameon=False)
    ax4.set_title("multilayer")
    ax4.set_xlabel("node")
    ax4.set_yticks([])

    fig.savefig(out_png, dpi=cfg["dpi"], bbox_inches="tight")
    plt.close(fig)
    return [out_png]


def _export_importance_csv(pair, outdir, node_to_name, harmonization, thr, normalization, modalities):
    rows = []
    all_nodes = set()
    for _, d in modalities.items():
        all_nodes |= set(d.keys())
    for node in sorted(all_nodes):
        row = {
            "pair": pair,
            "harmonization": harmonization,
            "node": int(node),
            "label": node_to_name.get(int(node)) if node_to_name else None,
            "threshold": float(thr),
            "normalization": normalization,
        }
        for name, d in modalities.items():
            val = d.get(node, np.nan)
            row[f"importance_{name}"] = val
            row[f"selected_{name}"] = bool((not pd.isna(val)) and (val >= thr))
        rows.append(row)
    out_csv = os.path.join(outdir, f"ROI_importance_all_modalities_{_slug_pair(pair)}_for_paper.csv")
    pd.DataFrame(rows).to_csv(out_csv, index=False)
    return out_csv


def main(cfg):
    if cfg.get("output_dir") is None:
        outdir = os.path.join(os.path.dirname(cfg["fi_csv"]), "OVERLAP_PLOTS_for_paper")
    else:
        outdir = cfg["output_dir"]
    ensure_outdir(outdir)

    df = pd.read_csv(cfg["fi_csv"])
    required = {"mode", "harmonization", "pair", "feature", "importance_mean"}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(f"Missing columns: {missing}")

    df = df[df["harmonization"] == cfg["harmonization"]].copy()
    if df.empty:
        raise ValueError("No rows after harmonization filtering")

    pairs = cfg["pairs"] if cfg.get("pairs") else sorted(df["pair"].dropna().unique().tolist())

    atlas_img = None
    node_to_label = None
    node_to_name = None
    if HAVE_NILEARN and cfg.get("atlas_nii") and os.path.exists(cfg["atlas_nii"]):
        atlas_img = nib.load(cfg["atlas_nii"])
        all_nodes = [n for n in df["feature"].dropna().apply(_extract_node).tolist() if n is not None]
        node_to_label, node_to_name = load_mapping(cfg["node_label_mapping"], cfg["label_base"], max(all_nodes) if all_nodes else 0)

    saved_all = []
    csv_all = []

    for pair in pairs:
        grouped_struct = _aggregate_by_mode(df, "structural", cfg["harmonization"], pair)
        grouped_func = _aggregate_by_mode(df, "functional", cfg["harmonization"], pair)
        grouped_sf = _aggregate_by_mode(df, "structural+functional", cfg["harmonization"], pair)
        grouped_ml = _aggregate_by_mode(df, "multilayer", cfg["harmonization"], pair)

        struct_u = _fuse_metrics(grouped_struct, cfg["normalization"]) if grouped_struct else {}
        func_u = _fuse_metrics(grouped_func, cfg["normalization"]) if grouped_func else {}
        sf_func_u, sf_struct_u = _unified_from_sf(grouped_sf, cfg["normalization"]) if grouped_sf else ({}, {})
        ml_u = _fuse_metrics(grouped_ml, cfg["normalization"]) if grouped_ml else {}

        thr = float(cfg["score_threshold"]) if cfg["score_threshold"] is not None else 0.0

        csv_path = _export_importance_csv(
            pair,
            outdir,
            node_to_name,
            cfg["harmonization"],
            thr,
            cfg["normalization"],
            {
                "structural": struct_u,
                "functional": func_u,
                "sf_structural": sf_struct_u,
                "sf_functional": sf_func_u,
                "multilayer": ml_u,
            },
        )
        csv_all.append(csv_path)

        structural_only = {n for n, s in struct_u.items() if s >= thr}
        functional_only = {n for n, s in func_u.items() if s >= thr}
        sf_struct_sel = {n for n, s in sf_struct_u.items() if s >= thr}
        sf_func_sel = {n for n, s in sf_func_u.items() if s >= thr}
        ml_sel = {n for n, s in ml_u.items() if s >= thr}

        sf_both = sf_struct_sel & sf_func_sel
        sf_struct = sf_struct_sel - sf_both
        sf_func = sf_func_sel - sf_both

        ml_fs = ml_sel & sf_struct_sel & sf_func_sel
        ml_f = (ml_sel & sf_func_sel) - sf_struct_sel
        ml_s = (ml_sel & sf_struct_sel) - sf_func_sel
        ml_none = ml_sel - (ml_fs | ml_f | ml_s)

        sets_dict = {
            "structural_only": structural_only,
            "functional_only": functional_only,
            "sf_struct": sf_struct,
            "sf_func": sf_func,
            "sf_both": sf_both,
            "ml_f": ml_f,
            "ml_s": ml_s,
            "ml_fs": ml_fs,
            "ml_none": ml_none,
        }

        out_png = os.path.join(outdir, f"OVERLAP_{cfg['harmonization']}_{_slug_pair(pair)}_for_paper.png")
        saved = _plot_four_panel(atlas_img, node_to_label, sets_dict, pair, out_png, cfg)
        saved_all.extend(saved)

        overlap_rows = []
        for lab, nodes in sets_dict.items():
            for n in sorted(nodes):
                overlap_rows.append({
                    "pair": pair,
                    "set": lab,
                    "node": int(n),
                    "label": node_to_name.get(int(n)) if node_to_name else None,
                })
        if overlap_rows:
            pd.DataFrame(overlap_rows).to_csv(
                os.path.join(outdir, f"OVERLAP_nodes_{cfg['harmonization']}_{_slug_pair(pair)}_for_paper.csv"),
                index=False,
            )

    manifest = {"inputs": cfg, "figures": saved_all, "csv_files": csv_all}
    with open(os.path.join(outdir, "manifest_overlap_for_paper.json"), "w", encoding="utf-8") as f:
        json.dump(manifest, f, indent=2)

    print(f"[OK] Saved {len(saved_all)} figures and {len(csv_all)} importance CSV files in: {outdir}")


if __name__ == "__main__":
    for k in ("fi_csv", "harmonization", "score_threshold"):
        if CONFIG.get(k, None) in (None, ""):
            raise SystemExit(f"CONFIG[{k!r}] is required")
    main(CONFIG)
