"""Step 6: MM/GBSA (gmx_MMPBSA) on the second half of every production replica.

GB model igb=5 with 0.15 M salt; crystal waters and ions are treated as solvent and
stripped. Per-residue decomposition runs only for the ligands listed in config.yaml
(it is much slower). Each replica writes result.json, which step 7 collects.
"""
import argparse
import json
import re

from common import WORK, load_config, ligand_names, run, stage_dir

# PBRadii=3 is mbondi2, the radii igb=5 was parameterised with (gmx_MMPBSA defaults to mbondi3).
GB_INPUT = """&general
sys_name="{name}_rep{rep}",
startframe={start},
interval={interval},
temperature={temp},
PBRadii=3,
/
&gb
igb={igb}, saltcon={salt},
/
"""
DECOMP_INPUT = """&decomp
idecomp=2, dec_verbose=0, print_res="within 4",
/
"""


def parse_delta_total(dat_text):
    """Return (mean, sd, sem) of DELTA TOTAL from FINAL_RESULTS_MMPBSA.dat."""
    block = dat_text.split("Delta (Complex - Receptor - Ligand)")[-1]
    m = re.search(r"^\s*ΔTOTAL\s+(-?[\d.]+)\s+[\d.]+\s+([\d.]+)\s+[\d.]+\s+([\d.]+)", block, re.M)
    if not m:
        raise RuntimeError("could not find ΔTOTAL in gmx_MMPBSA output")
    return tuple(float(x) for x in m.groups())


def run_replica(name, rep, cfg, roots, frame_ps, discard_ns, decomp):
    sys_root, md_root, mm_root = roots
    sysdir = sys_root / name
    mddir = md_root / name / f"rep{rep}"
    d = mm_root / name / f"rep{rep}"
    d.mkdir(parents=True, exist_ok=True)
    if (d / "result.json").exists():
        return json.loads((d / "result.json").read_text())
    if not (mddir / "prod_center.xtc").exists():
        raise FileNotFoundError(f"{mddir}/prod_center.xtc missing; run step 5 first")

    mm = cfg["mmgbsa"]
    text = GB_INPUT.format(name=name, rep=rep, start=int(discard_ns * 1000 / frame_ps) + 1,
                           interval=mm["stride_frames"], temp=cfg["temperature_K"],
                           igb=mm["igb"], salt=mm["saltcon_M"])
    if decomp:
        text += DECOMP_INPUT
    (d / "mmgbsa.in").write_text(text)

    launcher = ["mpirun", "-np", mm["mpi_procs"]] if mm["mpi_procs"] > 1 else []
    cmd = [*launcher, "gmx_MMPBSA", "-O", "-i", "mmgbsa.in",
           "-cs", mddir / "prod.tpr", "-ci", sysdir / "index.ndx", "-cg", "1", "2",
           "-ct", mddir / "prod_center.xtc", "-cp", sysdir / "complex.top",
           "-o", "FINAL_RESULTS_MMPBSA.dat", "-eo", "FINAL_RESULTS_MMPBSA.csv", "-nogui"]
    if decomp:
        cmd += ["-do", "FINAL_DECOMP_MMPBSA.dat", "-deo", "FINAL_DECOMP_MMPBSA.csv"]
    run(cmd, cwd=d, log=d / "gmx_MMPBSA.log")

    mean, sd, sem = parse_delta_total((d / "FINAL_RESULTS_MMPBSA.dat").read_text())
    for tmp in d.glob("_GMXMMPBSA_*"):  # ~10 MB of per-frame scratch files per replica
        tmp.unlink()
    result = {"ligand": name, "replica": rep, "dG_mmgbsa": mean, "sd": sd, "sem": sem}
    (d / "result.json").write_text(json.dumps(result, indent=2))
    return result


def main():
    cfg = load_config()
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--ligands", nargs="+")
    p.add_argument("--replicas", nargs="+", type=int)
    p.add_argument("--pose", choices=["docked", "reference"], default=cfg["build"]["start_pose"])
    p.add_argument("--quick", action="store_true", help="analyse the step-5 --quick runs")
    args = p.parse_args()

    md = cfg["md"]
    roots = tuple(stage_dir(s, args.pose) for s in ("systems", "md", "mmgbsa"))
    frame_ps, discard_ns = md["xtc_interval_ps"], cfg["mmgbsa"]["discard_ns"]
    if args.quick:
        roots = (stage_dir("systems"), WORK / "md_quick", WORK / "mmgbsa_quick")
        frame_ps, discard_ns = 2, 0.05
    replicas = args.replicas or list(range(1, md["replicas"] + 1))
    for name in args.ligands or ligand_names(cfg):
        for rep in replicas:
            r = run_replica(name, rep, cfg, roots, frame_ps, discard_ns,
                            decomp=name in cfg["mmgbsa"]["decomp_ligands"])
            print(f"{name:12s} rep{rep}  ΔG(MM/GBSA) = {r['dG_mmgbsa']:7.2f} ± {r['sem']:.2f} kcal/mol")


if __name__ == "__main__":
    main()
