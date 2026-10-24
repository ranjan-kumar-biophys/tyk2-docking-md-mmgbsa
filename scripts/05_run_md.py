"""Step 5: run EM -> NVT -> NPT -> production for every ligand and replica (GROMACS).

EM is done once per ligand; each replica then gets its own velocity seed. Finished
stages are skipped and production resumes from its checkpoint, so the script can simply
be re-run after a crash or a killed job. Use --ligands / --replicas / --gpu-id to split
the work across GPUs, --pose reference for the benchmark-pose control runs, and --quick
for a short smoke test (written to work/md_quick).
"""
import argparse
import os
import re
import shlex
from string import Template

from common import MDP, WORK, load_config, ligand_names, run, stage_dir

GMX = os.environ.get("GMX", "gmx")


def render(stage, out, **values):
    out.write_text(Template((MDP / f"{stage}.mdp").read_text()).substitute(**values))


def grompp_mdrun(d, name, coords, sysdir, mdrun_flags, restraint=None, checkpoint=None):
    if (d / f"{name}.gro").exists():
        print(f"  {d.name}/{name}: done, skipping")
        return
    cmd = [GMX, "grompp", "-f", f"{name}.mdp", "-c", coords, "-p", sysdir / "complex.top",
           "-n", sysdir / "index.ndx", "-o", f"{name}.tpr", "-maxwarn", "0"]
    if restraint:
        cmd += ["-r", restraint]
    if checkpoint:
        cmd += ["-t", checkpoint]
    if not (d / f"{name}.tpr").exists():
        run(cmd, cwd=d, log=d / f"{name}_grompp.log")
    run([GMX, "mdrun", "-deffnm", name, "-cpi", f"{name}.cpt", *mdrun_flags], cwd=d,
        log=d / f"{name}_mdrun.log")
    perf = re.search(r"Performance:\s+([\d.]+)", (d / f"{name}.log").read_text())
    if perf:
        print(f"  {d.name}/{name}: {perf.group(1)} ns/day")


def em_flags(flags):
    """Steepest descent supports only nonbonded work on the GPU: drop -update/-pme/-bonded."""
    out, skip = [], False
    for f in flags:
        if skip:
            skip = False
        elif f in ("-update", "-pme", "-bonded"):
            skip = True
        else:
            out.append(f)
    return out


def run_ligand(name, replicas, md, temp, mdrun_flags, root, pose, em_only=False):
    sysdir = stage_dir("systems", pose) / name
    lig_root = root / name
    em = lig_root / "em"
    em.mkdir(parents=True, exist_ok=True)
    (em / "em.mdp").write_text((MDP / "em.mdp").read_text())
    grompp_mdrun(em, "em", sysdir / "complex.gro", sysdir, em_flags(mdrun_flags))
    if em_only:
        return

    dt = md["dt_ps"]
    for rep in replicas:
        d = lig_root / f"rep{rep}"
        d.mkdir(exist_ok=True)
        common = dict(DT=dt, TEMP=temp)
        render("nvt", d / "nvt.mdp", NSTEPS=round(md["nvt_ps"] / dt), SEED=1000 + rep, **common)
        render("npt", d / "npt.mdp", NSTEPS=round(md["npt_ps"] / dt), **common)
        render("prod", d / "prod.mdp", NSTEPS=round(md["production_ns"] * 1000 / dt),
               NSTXTC=round(md["xtc_interval_ps"] / dt), **common)
        grompp_mdrun(d, "nvt", em / "em.gro", sysdir, mdrun_flags, restraint=em / "em.gro")
        grompp_mdrun(d, "npt", "nvt.gro", sysdir, mdrun_flags, restraint="nvt.gro",
                     checkpoint="nvt.cpt")
        grompp_mdrun(d, "prod", "npt.gro", sysdir, mdrun_flags, checkpoint="npt.cpt")
        if not (d / "prod_center.xtc").exists():
            run([GMX, "trjconv", "-s", "prod.tpr", "-f", "prod.xtc", "-n", sysdir / "index.ndx",
                 "-o", "prod_center.xtc", "-pbc", "mol", "-center"], cwd=d,
                stdin="Protein\nSystem\n", log=d / "trjconv.log")


def main():
    cfg = load_config()
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--ligands", nargs="+", help="subset of ligands (default: config)")
    p.add_argument("--replicas", nargs="+", type=int, help="replica numbers, e.g. 1 2 3")
    p.add_argument("--gpu-id", help="pass -gpu_id to mdrun")
    p.add_argument("--extra-flags", default="", help='more mdrun flags, e.g. "-ntomp 12 -pinoffset 12"')
    p.add_argument("--em-only", action="store_true",
                   help="only minimise (lets parallel workers split replicas without racing on EM)")
    p.add_argument("--pose", choices=["docked", "reference"], default=cfg["build"]["start_pose"])
    p.add_argument("--quick", action="store_true", help="10 ps NVT, 20 ps NPT, 100 ps production")
    args = p.parse_args()

    md = dict(cfg["md"])
    root = stage_dir("md", args.pose)
    if args.quick:
        md.update(nvt_ps=10, npt_ps=20, production_ns=0.1, xtc_interval_ps=2)
        root = WORK / "md_quick"
    flags = shlex.split(md.get("mdrun_flags") or "") + shlex.split(args.extra_flags)
    if args.gpu_id is not None:
        flags += ["-gpu_id", args.gpu_id]
    names = args.ligands or ligand_names(cfg)
    replicas = args.replicas or list(range(1, md["replicas"] + 1))

    for name in names:
        print(f"== {name}")
        run_ligand(name, replicas, md, cfg["temperature_K"], flags, root, args.pose, args.em_only)


if __name__ == "__main__":
    main()
