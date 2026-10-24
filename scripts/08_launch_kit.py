"""Step 8: turn the final tables into the README results, a LinkedIn post and a case-study PDF.

Reads results/summary.csv + results/stats.csv (step 7) and writes:
  README.md                  answer line + results section (between the HTML markers)
  results/linkedin_post.md   post draft, numbers filled in, links left as [placeholders]
  results/case_study.pdf     one A4 page to attach to proposals and outreach emails
The wording follows the numbers: it does not assume MM/GBSA wins.
"""
import re
import textwrap

import matplotlib

matplotlib.use("Agg")
import matplotlib.image as mpimg
import matplotlib.pyplot as plt
import pandas as pd

from common import FIGURES, RESULTS, ROOT

INK, INK_2 = "#0b0b0b", "#52514e"
METHOD_KEYS = {"Docking (Vina)": "dock", "MD + MM/GBSA (docked poses)": "mm",
               "MD + MM/GBSA (docking failures re-run)": "fixed"}


def md_table(df):
    """Markdown table without the optional tabulate dependency."""
    cell = lambda v: "" if pd.isna(v) else f"{v:.2f}" if isinstance(v, float) else str(v)
    rows = ["| " + " | ".join(df.columns) + " |", "|" + "---|" * len(df.columns)]
    rows += ["| " + " | ".join(cell(v) for v in r) + " |" for r in df.itertuples(index=False)]
    return "\n".join(rows)


def fmt(s):
    return (f"{s.pearson_r:.2f} [{s.r_ci_low:.2f}, {s.r_ci_high:.2f}]",
            f"{s.spearman_rho:.2f} [{s.rho_ci_low:.2f}, {s.rho_ci_high:.2f}]")


def findings(summary, stats):
    st = {METHOD_KEYS[m]: row for m, row in stats.set_index("method").iterrows() if m in METHOD_KEYS}
    n = len(summary)
    ok = int((summary["top_pose_rmsd_A"] < 2).sum())
    bad = summary.loc[summary["top_pose_rmsd_A"] >= 2, "ligand"].tolist()
    dock, mm, fixed = st.get("dock"), st.get("mm"), st.get("fixed")

    lines = [f"Vina placed {ok}/{n} ligands within 2 Å of the benchmark pose."]
    if dock is not None and mm is not None:
        better = mm.pearson_r > dock.pearson_r
        overlap = mm.r_ci_low <= dock.r_ci_high and dock.r_ci_low <= mm.r_ci_high
        lines.append(
            f"MD + MM/GBSA from the docked poses gave Pearson r = {mm.pearson_r:.2f} against "
            f"experiment, versus r = {dock.pearson_r:.2f} for the docking score"
            + (" (the 95% CIs overlap, so with 13 compounds this is a trend, not a proof)." if overlap
               else (", a clear improvement." if better else ", so no improvement.")))
    if fixed is not None and mm is not None and bad:
        delta = fixed.pearson_r - mm.pearson_r
        verdict = ("so wrong docking poses were a large part of the MM/GBSA error." if delta > 0.1
                   else "so wrong docking poses cost some accuracy." if delta > 0.03
                   else "so wrong docking poses made little difference to the ranking.")
        lines.append(f"Re-running the {len(bad)} misdocked ligands from the benchmark pose changed "
                     f"r to {fixed.pearson_r:.2f} ({delta:+.2f}), {verdict}")
    return lines, st, ok, n, bad


def readme_block(summary, stats, lines, ok, n):
    rows = ["| Method | Pearson r [95% CI] | Spearman ρ [95% CI] |", "|---|---|---|"]
    rows += [f"| {r.method} | {fmt(r)[0]} | {fmt(r)[1]} |" for r in stats.itertuples()]
    cols = ["ligand", "Ki_nM", "dG_exp", "vina_score", "top_pose_rmsd_A", "dG_mmgbsa", "dG_mmgbsa_sem"]
    cols += [c for c in ("dG_mmgbsa_ref", "dG_mmgbsa_ref_sem") if c in summary]
    per_lig = md_table(summary[cols])
    return "\n".join([
        *[f"- {l}" for l in lines], "",
        *rows, "",
        f"Docking top-pose accuracy: {ok}/{n} within 2 Å (Figure 2).", "",
        "![MM/GBSA per-residue decomposition](figures/fig3_residue_decomposition.png)", "",
        "<details><summary>Per-ligand numbers (kcal/mol; Ki in nM)</summary>", "",
        per_lig, "", "</details>",
    ])


def replace_between(text, tag, new):
    pattern = re.compile(rf"(<!-- {tag}:START -->\n).*?(\n<!-- {tag}:END -->)", re.S)
    if not pattern.search(text):
        raise RuntimeError(f"README has no {tag} markers")
    return pattern.sub(lambda m: m.group(1) + new + m.group(2), text)


def linkedin_post(lines, st, ok, n):
    dock, mm, fixed = st.get("dock"), st.get("mm"), st.get("fixed")
    hook = (f"Docking score vs MD + MM/GBSA: r = {dock.pearson_r:.2f} vs {mm.pearson_r:.2f} "
            "on 13 TYK2 inhibitors. Here's what I learned.")
    body = [
        hook, "",
        "I ran a small benchmark on a question that comes up in almost every hit-to-lead project: "
        "is it worth running MD after docking, or is the docking score good enough to rank "
        "compounds?", "",
        "Setup (all open-source, fully reproducible):",
        "• 13 TYK2 kinase inhibitors with measured Ki (public OpenFF benchmark)",
        "• AutoDock Vina docking from fresh conformers",
        "• GROMACS MD, 3 independent 10 ns replicas per complex, GAFF2/AM1-BCC ligands",
        "• MM/GBSA (gmx_MMPBSA) + per-residue decomposition",
        "• Control: misdocked ligands re-run from the correct pose", "",
        "What came out:",
        *[f"→ {l}" for l in lines], "",
        "Takeaway: check your poses before you trust any rescoring. With 13 compounds the "
        "confidence intervals are wide, and I've reported them rather than hide them.", "",
        "Next: alchemical FEP on the same series, to see how much rigorous free energies add.", "",
        "Code + figures: [GitHub link]",
        "I take on MD / free-energy projects for biotech teams and labs. DM me if you need this "
        "for your series.", "",
        "#ComputationalChemistry #DrugDiscovery #MolecularDynamics #GROMACS #Biophysics",
    ]
    return "\n".join(body)


def case_study_pdf(lines, out):
    fig = plt.figure(figsize=(8.27, 11.69))  # A4 portrait; y runs top -> bottom via a cursor
    left, line_h = 0.07, 0.0165
    y = 0.965

    def text(s, size=9.5, bold=False, color=INK, indent=""):
        nonlocal y
        for i, part in enumerate(textwrap.wrap(s, 100) or [""]):
            fig.text(left, y, (indent if i else "") + part, fontsize=size, color=color, va="top",
                     fontweight="bold" if bold else "normal")
            y -= line_h * size / 9.5
        return y

    def heading(s):
        nonlocal y
        y -= 0.012
        text(s, size=11, bold=True)
        y -= 0.004

    def image(path, width):
        nonlocal y
        img = mpimg.imread(path)
        h = width * img.shape[0] / img.shape[1] * 8.27 / 11.69
        ax = fig.add_axes([left, y - h, width, h])
        ax.imshow(img)
        ax.axis("off")
        y -= h + 0.01

    text("Case study: does MD rescoring beat docking?", size=17, bold=True)
    y -= 0.008
    text("TYK2 kinase · 13 inhibitors with measured Ki · open-source pipeline", size=10, color=INK_2)
    heading("Question")
    text("For a congeneric inhibitor series, does MD + MM/GBSA rank compounds better than the "
         "docking score?")
    y -= 0.006
    image(FIGURES / "fig1_correlation.png", 0.86)
    heading("Results")
    for l in lines:
        text("• " + l, indent="   ")
        y -= 0.004
    if (FIGURES / "fig3_residue_decomposition.png").exists():
        y -= 0.006
        image(FIGURES / "fig3_residue_decomposition.png", 0.62)
    heading("Method")
    text("AutoDock Vina docking → GAFF2/AM1-BCC ligand parameters → GROMACS MD (ff14SB, TIP3P, "
         "0.15 M NaCl; 3 × 10 ns replicas per complex) → MM/GBSA (gmx_MMPBSA, igb=5, mbondi2) "
         "with per-residue decomposition. Bootstrap 95% confidence intervals on all correlations.")
    heading("What you get when I run this for your series")
    for o in ("Ranked compounds with error bars, not single numbers",
              "Binding-mode check that flags poses docking got wrong",
              "Per-residue maps of what drives binding",
              "Paper-ready methods, figures and all input files"):
        text("✔ " + o)

    fig.text(left, 0.035, "Ranjan Kumar · Computational biophysicist · ranjankumarbijarnia70@gmail.com · "
             "github.com/ranjan-kumar-biophys/tyk2-docking-md-mmgbsa", fontsize=9, color=INK_2)
    fig.savefig(out)
    plt.close(fig)


def main():
    summary = pd.read_csv(RESULTS / "summary.csv")
    stats = pd.read_csv(RESULTS / "stats.csv")
    lines, st, ok, n, bad = findings(summary, stats)

    readme = (ROOT / "README.md").read_text()
    answer = "**Answer:** " + " ".join(lines[1:] or lines) + " See [Results](#results)."
    readme = replace_between(readme, "ANSWER", answer)
    readme = replace_between(readme, "RESULTS", readme_block(summary, stats, lines, ok, n))
    (ROOT / "README.md").write_text(readme)

    if "dock" in st and "mm" in st:
        (RESULTS / "linkedin_post.md").write_text(linkedin_post(lines, st, ok, n) + "\n")
    case_study_pdf(lines, RESULTS / "case_study.pdf")
    print("\n".join(lines))
    print("wrote README.md, results/linkedin_post.md, results/case_study.pdf")


if __name__ == "__main__":
    main()
