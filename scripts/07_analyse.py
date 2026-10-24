"""Step 7: compare docking scores and MM/GBSA against experiment; make tables + figures.

Works after docking alone (MM/GBSA columns are left empty until step 6 has run).
Correlations carry bootstrap 95% CIs because n = 13 is small and the CIs are wide;
saying so honestly is part of the result.
"""
import argparse
import json

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats

from common import (BENCH, FIGURES, RESULTS, WORK, ki_to_dg, load_affinities, load_config,
                    ligand_names, stage_dir)

SURFACE, INK, INK_2, GRID = "#fcfcfb", "#0b0b0b", "#52514e", "#e4e3df"
COLOR = {"vina_score": "#2a78d6", "dG_mmgbsa": "#eb6834", "dG_mmgbsa_fixed": "#1baf7a"}
LABEL = {"vina_score": "Vina docking score (kcal/mol)", "dG_mmgbsa": "MM/GBSA ΔG (kcal/mol)",
         "dG_mmgbsa_fixed": "MM/GBSA ΔG (kcal/mol)"}
TITLE = {"vina_score": "Docking (Vina)", "dG_mmgbsa": "MD + MM/GBSA (docked poses)",
         "dG_mmgbsa_fixed": "MD + MM/GBSA (docking failures re-run)"}

plt.rcParams.update({
    "figure.facecolor": SURFACE, "axes.facecolor": SURFACE, "savefig.facecolor": SURFACE,
    "axes.edgecolor": INK_2, "axes.labelcolor": INK, "text.color": INK,
    "xtick.color": INK_2, "ytick.color": INK_2, "axes.grid": True, "grid.color": GRID,
    "grid.linewidth": 0.8, "axes.axisbelow": True,
    "axes.spines.top": False, "axes.spines.right": False,
    "font.size": 10,
})


def correlation(x, y, n_boot=10000, seed=0):
    x, y = np.asarray(x), np.asarray(y)
    rng = np.random.default_rng(seed)
    boot_r, boot_rho = [], []
    for idx in rng.integers(0, len(x), size=(n_boot, len(x))):
        if np.ptp(x[idx]) == 0 or np.ptp(y[idx]) == 0:
            continue
        boot_r.append(stats.pearsonr(x[idx], y[idx])[0])
        boot_rho.append(stats.spearmanr(x[idx], y[idx])[0])
    lo_r, hi_r = np.percentile(boot_r, [2.5, 97.5])
    lo_s, hi_s = np.percentile(boot_rho, [2.5, 97.5])
    return {"n": len(x), "pearson_r": stats.pearsonr(x, y)[0], "r_ci_low": lo_r, "r_ci_high": hi_r,
            "spearman_rho": stats.spearmanr(x, y)[0], "rho_ci_low": lo_s, "rho_ci_high": hi_s,
            "kendall_tau": stats.kendalltau(x, y)[0]}


def collect_mmgbsa(names, mm_root, replicas):
    rows = []
    for name in names:
        vals = [json.loads(f.read_text())["dG_mmgbsa"]
                for rep in replicas if (f := mm_root / name / f"rep{rep}" / "result.json").exists()]
        if vals:
            sem = np.std(vals, ddof=1) / np.sqrt(len(vals)) if len(vals) > 1 else np.nan
            rows.append({"ligand": name, "dG_mmgbsa": np.mean(vals), "dG_mmgbsa_sem": sem,
                         "n_replicas": len(vals)})
    return pd.DataFrame(rows, columns=["ligand", "dG_mmgbsa", "dG_mmgbsa_sem", "n_replicas"])


def tyk2_residue_numbers():
    """gmx_MMPBSA numbers receptor residues 1..N in file order; map back to 4GIH numbering."""
    seen = []
    for line in (BENCH / "protein.pdb").read_text().splitlines():
        if line.startswith(("ATOM", "HETATM")) and line[17:20] != "HOH":
            key = (line[17:20], int(line[22:26]))
            if not seen or seen[-1] != key:
                seen.append(key)
    return {i: f"{name} {num}" for i, (name, num) in enumerate(seen, start=1)}


def read_decomp(csv_path):
    """Mean per-residue ΔTOTAL (kcal/mol) for receptor residues from FINAL_DECOMP_MMPBSA.csv."""
    lines = csv_path.read_text().splitlines()
    start = lines.index("DELTAS:") + 3  # skip the TDC title and the column header
    rows = []
    for line in lines[start:]:
        parts = line.split(",")
        if len(parts) < 8:
            break
        if parts[1].startswith("R:"):
            rows.append((int(parts[1].split(":")[-1]), float(parts[-1])))
    return pd.DataFrame(rows, columns=["res", "total"]).groupby("res")["total"].mean()


def decomposition_figure(names, mm_root, replicas, top_n=12):
    per_lig = {}
    for name in names:
        reps = [read_decomp(f) for rep in replicas
                if (f := mm_root / name / f"rep{rep}" / "FINAL_DECOMP_MMPBSA.csv").exists()]
        if reps:
            per_lig[name] = pd.concat(reps, axis=1).mean(axis=1)
    if not per_lig:
        return
    labels = tyk2_residue_numbers()
    fig, axes = plt.subplots(1, len(per_lig), figsize=(max(4.6 * len(per_lig), 7.5), 4.2),
                             squeeze=False)
    for ax, (name, s) in zip(axes[0], per_lig.items()):
        top = s.nsmallest(top_n).iloc[::-1]
        ax.barh([labels[i] for i in top.index], top.values, color=COLOR["dG_mmgbsa"], height=0.7,
                edgecolor=SURFACE, linewidth=2)
        ax.set_title(name, loc="left", fontsize=11, fontweight="bold")
        ax.set_xlabel("Contribution to ΔG (kcal/mol)")
        ax.tick_params(axis="y", labelsize=8)
        ax.grid(axis="y", visible=False)
    fig.suptitle(f"MM/GBSA decomposition: top {top_n} stabilising TYK2 residues",
                 x=0.01, ha="left", fontsize=12, fontweight="bold")
    fig.tight_layout()
    fig.savefig(FIGURES / "fig3_residue_decomposition.png", dpi=200)


def scatter_panel(ax, df, col, st):
    err = df.get(f"{col}_sem")
    if err is not None and err.notna().any():
        ax.errorbar(df["dG_exp"], df[col], yerr=err, fmt="none", ecolor=INK_2, elinewidth=1,
                    capsize=0, zorder=2)
    rerun = df["rerun"] if col == "dG_mmgbsa_fixed" else pd.Series(False, index=df.index)
    ax.scatter(df.loc[~rerun, "dG_exp"], df.loc[~rerun, col], s=64, color=COLOR[col],
               edgecolor=SURFACE, linewidth=2, zorder=3)
    if rerun.any():  # hollow markers = started from the benchmark pose instead of the docked one
        ax.scatter(df.loc[rerun, "dG_exp"], df.loc[rerun, col], s=64, facecolor=SURFACE,
                   edgecolor=COLOR[col], linewidth=2, zorder=3)
        ax.text(0.97, 0.03, "○ re-run from benchmark pose", transform=ax.transAxes, ha="right",
                va="bottom", fontsize=8, color=INK_2)
    slope, intercept = np.polyfit(df["dG_exp"], df[col], 1)
    xs = np.linspace(df["dG_exp"].min() - 0.2, df["dG_exp"].max() + 0.2, 2)
    ax.plot(xs, slope * xs + intercept, color=INK_2, linewidth=1, linestyle="--", zorder=1)
    for _, r in df.loc[[df["dG_exp"].idxmin(), df["dG_exp"].idxmax()]].iterrows():
        ax.annotate(r["ligand"], (r["dG_exp"], r[col]), xytext=(6, 4), textcoords="offset points",
                    fontsize=8, color=INK_2)
    ax.set_title(TITLE[col], loc="left", fontsize=11, fontweight="bold", pad=30)
    ax.text(0, 1.02,  # stats sit between title and plot so they never cover a point
            f"Pearson r = {st['pearson_r']:.2f}  [{st['r_ci_low']:.2f}, {st['r_ci_high']:.2f}]\n"
            f"Spearman ρ = {st['spearman_rho']:.2f}  [{st['rho_ci_low']:.2f}, {st['rho_ci_high']:.2f}]",
            transform=ax.transAxes, va="bottom", fontsize=8.5, color=INK_2)
    ax.set_xlabel("Experimental ΔG (kcal/mol)")
    ax.set_ylabel(LABEL[col])


def main():
    cfg = load_config()
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--quick", action="store_true", help="use the step-6 --quick results")
    args = p.parse_args()

    names = ligand_names(cfg)
    temp = cfg["temperature_K"]
    ki = load_affinities()
    df = pd.DataFrame({"ligand": names, "Ki_nM": [ki[n] * 1e9 for n in names],
                       "dG_exp": [ki_to_dg(ki[n], temp) for n in names]})
    df = df.merge(pd.read_csv(RESULTS / "docking.csv"), on="ligand", how="left")
    reps = range(1, cfg["md"]["replicas"] + 1)
    mm_root = WORK / "mmgbsa_quick" if args.quick else stage_dir("mmgbsa")
    df = df.merge(collect_mmgbsa(names, mm_root, reps), on="ligand", how="left")

    # Control: ligands docking got wrong, re-run from the benchmark pose (04-06 --pose reference).
    ref = collect_mmgbsa(names, stage_dir("mmgbsa", "reference"), reps)
    if not args.quick and len(ref):
        ref.columns = ["ligand", "dG_mmgbsa_ref", "dG_mmgbsa_ref_sem", "n_replicas_ref"]
        df = df.merge(ref, on="ligand", how="left")
        df["rerun"] = df["dG_mmgbsa_ref"].notna()
        df["dG_mmgbsa_fixed"] = df["dG_mmgbsa_ref"].where(df["rerun"], df["dG_mmgbsa"])
        df["dG_mmgbsa_fixed_sem"] = df["dG_mmgbsa_ref_sem"].where(df["rerun"], df["dG_mmgbsa_sem"])
    df.to_csv(RESULTS / "summary.csv", index=False, float_format="%.3f")

    methods = [c for c in ("vina_score", "dG_mmgbsa", "dG_mmgbsa_fixed")
               if c in df and df[c].notna().sum() >= 3]
    st = {}
    for col in methods:
        sub = df.dropna(subset=[col])
        st[col] = correlation(sub["dG_exp"], sub[col])
    pd.DataFrame([{"method": TITLE[c], **s} for c, s in st.items()]).to_csv(
        RESULTS / "stats.csv", index=False, float_format="%.3f")

    FIGURES.mkdir(exist_ok=True)
    if methods:
        fig, axes = plt.subplots(1, len(methods), figsize=(5.2 * len(methods), 4.4), squeeze=False)
        for ax, col in zip(axes[0], methods):
            scatter_panel(ax, df.dropna(subset=[col]), col, st[col])
        fig.suptitle(f"TYK2, {len(df)} ligands: predicted vs experimental binding free energy",
                     x=0.01, ha="left", fontsize=12, fontweight="bold")
        fig.tight_layout()
        fig.savefig(FIGURES / "fig1_correlation.png", dpi=200)
    else:
        print("fewer than 3 ligands with results: skipping correlation figure")

    fig, ax = plt.subplots(figsize=(7, 3.6))
    d = df.sort_values("top_pose_rmsd_A")
    ax.bar(d["ligand"], d["top_pose_rmsd_A"], color=COLOR["vina_score"], width=0.7,
           edgecolor=SURFACE, linewidth=2)
    ax.axhline(2.0, color=INK_2, linewidth=1, linestyle="--")
    ax.set_ylim(0, max(2.6, d["top_pose_rmsd_A"].max() * 1.1))
    ax.text(-0.4, 2.05, "2 Å success threshold", ha="left", va="bottom", fontsize=8, color=INK_2)
    ok = (d["top_pose_rmsd_A"] < 2).sum()
    ax.set_title(f"Docking pose accuracy: top pose within 2 Å for {ok}/{len(d)} ligands",
                 loc="left", fontsize=11, fontweight="bold")
    ax.set_ylabel("Heavy-atom RMSD to benchmark pose (Å)")
    ax.tick_params(axis="x", rotation=45, labelsize=8)
    ax.grid(axis="x", visible=False)
    fig.tight_layout()
    fig.savefig(FIGURES / "fig2_docking_pose_rmsd.png", dpi=200)

    decomposition_figure([n for n in cfg["mmgbsa"]["decomp_ligands"] if n in names], mm_root,
                         range(1, cfg["md"]["replicas"] + 1))

    print(df.to_string(index=False, float_format=lambda v: f"{v:.2f}"))
    for col, s in st.items():
        print(f"{TITLE[col]:16s} r = {s['pearson_r']:.2f} [{s['r_ci_low']:.2f}, {s['r_ci_high']:.2f}]"
              f"  rho = {s['spearman_rho']:.2f}")


if __name__ == "__main__":
    main()
