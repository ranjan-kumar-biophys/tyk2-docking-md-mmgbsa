"""Step 11: MD analysis figures from the step-9 per-replica results.

  fig4_md_stability.png   backbone RMSD and ligand RMSD to the benchmark pose over time
                          (correctly docked, misdocked, and benchmark-pose control runs)
  fig5_rmsf.png           Cα RMSF with the binding-site residues and kinase motifs marked
  fig6_fingerprint.png    contact occupancy per residue for all ligands (strongest at top),
                          with H-bond occupancy as dots: the series' binding fingerprint
  fig7_rose_<ligand>.png  polar "rose" chart of contacts and H-bonds for one ligand
  fig8_pca.png            free-energy landscape on the first two Cα principal components,
                          strongest vs weakest binder, on shared axes
"""
import json
import math

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from common import FIGURES, RESULTS, WORK, ki_to_dg, load_affinities, load_config, ligand_names

SURFACE, INK, INK_2, GRID = "#fcfcfb", "#0b0b0b", "#52514e", "#e4e3df"
BLUE, ORANGE, AQUA = "#2a78d6", "#eb6834", "#1baf7a"
MOTIFS = {"P-loop": (903, 911), "hinge": (978, 984), "DFG": (1041, 1043)}
plt.rcParams.update({
    "figure.facecolor": SURFACE, "axes.facecolor": SURFACE, "savefig.facecolor": SURFACE,
    "axes.edgecolor": INK_2, "axes.labelcolor": INK, "text.color": INK, "xtick.color": INK_2,
    "ytick.color": INK_2, "axes.grid": True, "grid.color": GRID, "grid.linewidth": 0.8,
    "axes.axisbelow": True, "axes.spines.top": False, "axes.spines.right": False, "font.size": 10,
})


def load_runs():
    """[(pose, ligand, rep, npz, occupancy json)] for every analysed replica."""
    runs = []
    for npz in sorted((WORK / "analysis").glob("*/*/rep*.npz")):
        pose, lig, rep = npz.parts[-3], npz.parts[-2], int(npz.stem[3:])
        runs.append((pose, lig, rep, np.load(npz), json.loads(npz.with_suffix(".json").read_text())))
    return runs


def misdocked():
    d = pd.read_csv(RESULTS / "docking.csv")
    return set(d.loc[d["top_pose_rmsd_A"] > 2, "ligand"])


def title(ax, text):
    ax.set_title(text, loc="left", fontsize=11, fontweight="bold")


def fig_stability(runs, bad):
    fig, (a, b) = plt.subplots(1, 2, figsize=(11, 4.2))
    for pose, lig, _, d, _ in runs:
        a.plot(d["time_ns"], d["rmsd_backbone"], color=INK_2, alpha=0.25, linewidth=0.8)
        color = AQUA if pose == "reference" else ORANGE if lig in bad else BLUE
        b.plot(d["time_ns"], d["rmsd_ligand_to_benchmark"], color=color, alpha=0.55, linewidth=1)
    title(a, "Protein backbone RMSD")
    a.set_xlabel("Time (ns)"); a.set_ylabel("RMSD (Å)")
    title(b, "Ligand RMSD to the benchmark pose")
    b.axhline(2.0, color=INK_2, linewidth=1, linestyle="--")
    b.set_xlabel("Time (ns)"); b.set_ylabel("RMSD (Å)")
    for color, text in ((BLUE, "docked correctly"), (ORANGE, "misdocked start"), (AQUA, "control: benchmark start")):
        b.plot([], [], color=color, linewidth=2, label=text)
    b.legend(frameon=False, fontsize=8.5, loc="center right")   # the empty 2.5-4.5 Å band
    fig.tight_layout()
    fig.savefig(FIGURES / "fig4_md_stability.png", dpi=200)


def fig_rmsf(runs, contacts):
    resid = runs[0][3]["ca_resid_4gih"]
    rmsf = np.array([d["rmsf_ca"] for pose, *_, d, _ in [(r[0], r[3], r[4]) for r in runs] if pose == "docked"])
    mean, sd = rmsf.mean(axis=0), rmsf.std(axis=0)
    fig, ax = plt.subplots(figsize=(11, 3.6))
    for name, (lo, hi) in MOTIFS.items():
        ax.axvspan(lo - 0.5, hi + 0.5, color=GRID, alpha=0.8, zorder=0)
        ax.text((lo + hi) / 2, mean.max() * 1.05, name, ha="center", fontsize=8.5, color=INK_2)
    ax.fill_between(resid, mean - sd, mean + sd, color=BLUE, alpha=0.18, linewidth=0)
    ax.plot(resid, mean, color=BLUE, linewidth=1.6)
    pocket = [int(k.split()[1]) for k, v in contacts.items() if v >= 0.5]
    on_pocket = np.isin(resid, pocket)
    ax.scatter(resid[on_pocket], mean[on_pocket], s=22, color=ORANGE, edgecolor=SURFACE, linewidth=1.5,
               zorder=3, label="binding-site residue (contact ≥ 50%)")
    title(ax, "Cα fluctuations (RMSF), mean ± SD over all docked-pose runs")
    ax.set_xlabel("TYK2 residue"); ax.set_ylabel("RMSF (Å)")
    ax.set_ylim(0, mean.max() * 1.15)
    ax.legend(frameon=False, fontsize=8.5, loc="upper right")
    fig.tight_layout()
    fig.savefig(FIGURES / "fig5_rmsf.png", dpi=200)


def occupancy_table(runs, kind, pose="docked"):
    rows = {}
    for p, lig, _, _, occ in runs:
        if p == pose:
            rows.setdefault(lig, []).append(occ[kind])
    return pd.DataFrame({lig: pd.DataFrame(v).mean() for lig, v in rows.items()}).T.fillna(0.0)


def fig_fingerprint(runs, order):
    contacts, hbonds = occupancy_table(runs, "contacts"), occupancy_table(runs, "hbonds")
    cols = [c for c in contacts.columns if contacts[c].mean() >= 0.3]
    cols = sorted(cols, key=lambda c: int(c.split()[1]))
    rows = [l for l in order if l in contacts.index]
    c = contacts.loc[rows, cols]
    fig, ax = plt.subplots(figsize=(max(6, 0.42 * len(cols) + 2), 0.38 * len(rows) + 1.8))
    im = ax.imshow(c.values, cmap="Blues", vmin=0, vmax=1, aspect="auto")
    h = hbonds.reindex(index=rows, columns=cols).fillna(0.0)
    yy, xx = np.nonzero(h.values >= 0.1)
    ax.scatter(xx, yy, s=10 + 70 * h.values[yy, xx], color=ORANGE, edgecolor=SURFACE, linewidth=1)
    ax.set_xticks(range(len(cols)), [x.title().replace(" ", "") for x in cols], rotation=60, fontsize=8)
    ax.set_yticks(range(len(rows)), [r.replace("lig_", "") for r in rows], fontsize=8.5)
    ax.grid(False)
    title(ax, "Binding fingerprint: contact occupancy (blue) and H-bond occupancy (orange dots)")
    fig.colorbar(im, ax=ax, fraction=0.025, pad=0.01, label="contact occupancy")
    ax.set_ylabel("ligand (strongest binder at top)")
    fig.tight_layout()
    fig.savefig(FIGURES / "fig6_fingerprint.png", dpi=200)
    return contacts.mean()


def fig_rose(runs, ligand):
    contacts, hbonds = occupancy_table(runs, "contacts"), occupancy_table(runs, "hbonds")
    if ligand not in contacts.index:
        return
    c = contacts.loc[ligand]
    c = c[c >= 0.25].sort_index(key=lambda s: s.str.split().str[1].astype(int))
    h = hbonds.loc[ligand].reindex(c.index).fillna(0.0) if ligand in hbonds.index else c * 0
    theta = np.linspace(0, 2 * math.pi, len(c), endpoint=False)
    width = 2 * math.pi / len(c) * 0.85
    fig = plt.figure(figsize=(6.2, 6.2))
    ax = fig.add_subplot(projection="polar")
    ax.set_theta_zero_location("N")   # sequence order runs clockwise from 12 o'clock
    ax.set_theta_direction(-1)
    ax.bar(theta, c.values, width=width, color=BLUE, alpha=0.35, edgecolor=SURFACE, linewidth=2, label="contact")
    ax.bar(theta, h.values, width=width * 0.55, color=ORANGE, edgecolor=SURFACE, linewidth=1.5, label="H-bond")
    ax.set_xticks(theta, [x.title().replace(" ", "") for x in c.index], fontsize=8)
    ax.set_yticks([0.25, 0.5, 0.75, 1.0], ["25%", "50%", "75%", "100%"], fontsize=7, color=INK_2)
    ax.set_ylim(0, 1.05)
    ax.set_title(f"{ligand}: who holds the ligand?", fontsize=11, fontweight="bold", pad=18)
    ax.legend(frameon=False, fontsize=8.5, loc="upper right", bbox_to_anchor=(1.12, 1.1))
    fig.tight_layout()
    fig.savefig(FIGURES / f"fig7_rose_{ligand}.png", dpi=200)


def fig_pca(runs, strong, weak):
    docked = [(lig, d["ca_coords"]) for pose, lig, _, d, _ in runs if pose == "docked"]
    X = np.concatenate([c.reshape(len(c), -1) for _, c in docked])
    mean = X.mean(axis=0)
    _, s, vt = np.linalg.svd(X - mean, full_matrices=False)
    var = s ** 2 / (s ** 2).sum()
    proj = {lig: (c.reshape(len(c), -1) - mean) @ vt[:2].T for lig, c in docked}
    allp = np.concatenate(list(proj.values()))
    lim = [(allp[:, i].min(), allp[:, i].max()) for i in (0, 1)]
    panels = [l for l in (strong, weak) if l in proj]
    fig, axes = plt.subplots(1, len(panels), figsize=(5.2 * len(panels), 4.4), squeeze=False)
    for ax, lig in zip(axes[0], panels):
        p = np.concatenate([v for k, v in proj.items() if k == lig])
        H, xe, ye = np.histogram2d(p[:, 0], p[:, 1], bins=40, range=lim)
        with np.errstate(divide="ignore"):
            G = -np.log(H / H.sum())
        G -= np.nanmin(G[np.isfinite(G)])
        im = ax.contourf(0.5 * (xe[1:] + xe[:-1]), 0.5 * (ye[1:] + ye[:-1]), G.T,
                         levels=np.linspace(0, 6, 13), cmap="viridis_r")
        title(ax, f"{lig} (ΔG$_{{exp}}$ rank {'strongest' if lig == strong else 'weakest'})")
        ax.set_xlabel(f"PC1 ({var[0]:.0%} of variance, Å)"); ax.set_ylabel(f"PC2 ({var[1]:.0%}, Å)")
        ax.set_xlim(*lim[0]); ax.set_ylim(*lim[1]); ax.grid(False)
    fig.colorbar(im, ax=axes[0].tolist(), fraction=0.03, pad=0.02, label="free energy (kT)")
    fig.savefig(FIGURES / "fig8_pca.png", dpi=200, bbox_inches="tight")


def main():
    cfg = load_config()
    runs = load_runs()
    if not runs:
        raise SystemExit("no step-9 results in work/analysis")
    ki = load_affinities()
    order = sorted(ligand_names(cfg), key=lambda n: ki_to_dg(ki[n], cfg["temperature_K"]))
    FIGURES.mkdir(exist_ok=True)
    fig_stability(runs, misdocked())
    mean_contacts = fig_fingerprint(runs, order)
    fig_rmsf(runs, mean_contacts.to_dict())
    fig_rose(runs, order[0])
    fig_pca(runs, order[0], order[-1])
    print(f"{len(runs)} replicas -> fig4-fig8 in {FIGURES}")


if __name__ == "__main__":
    main()
