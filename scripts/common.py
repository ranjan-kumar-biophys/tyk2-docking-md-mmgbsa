"""Shared paths, config loading and helpers used by every pipeline step."""
import math
import shlex
import subprocess
from pathlib import Path

import yaml

ROOT = Path(__file__).resolve().parents[1]
BENCH = ROOT / "inputs" / "benchmark"
WORK = ROOT / "work"
RESULTS = ROOT / "results"
FIGURES = ROOT / "figures"
MDP = ROOT / "mdp"

R_KCAL = 0.0019872043  # gas constant, kcal/(mol K)


def stage_dir(stage, pose="docked"):
    """work/<stage> for the main (docked-pose) pipeline, work/<stage>_<pose> for control runs."""
    return WORK / (stage if pose == "docked" else f"{stage}_{pose}")


def load_config():
    with open(ROOT / "config.yaml") as fh:
        return yaml.safe_load(fh)


def load_affinities():
    """Benchmark ligands.yml -> {ligand: Ki in molar}."""
    with open(BENCH / "ligands.yml") as fh:
        data = yaml.safe_load(fh)
    scale = {"uM": 1e-6, "nM": 1e-9}
    out = {}
    for name, entry in data.items():
        m = entry["measurement"]
        if m["type"] != "ki":
            raise ValueError(f"{name}: expected Ki, got {m['type']}")
        out[name] = m["value"] * scale[m["unit"]]
    return out


def ki_to_dg(ki_molar, temperature):
    return R_KCAL * temperature * math.log(ki_molar)


def ligand_names(cfg):
    names = sorted(load_affinities())
    if cfg["ligands"] == "all":
        return names
    unknown = set(cfg["ligands"]) - set(names)
    if unknown:
        raise ValueError(f"unknown ligands in config: {sorted(unknown)}")
    return list(cfg["ligands"])


def run(cmd, cwd=None, stdin=None, log=None):
    """Run a command, echo it, and fail loudly with the tail of its output."""
    if isinstance(cmd, str):
        cmd = shlex.split(cmd)
    print("  $", " ".join(str(c) for c in cmd), flush=True)
    proc = subprocess.run(
        [str(c) for c in cmd], cwd=cwd, input=stdin, text=True,
        stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
    )
    if log:
        Path(log).write_text(proc.stdout)
    if proc.returncode != 0:
        tail = "\n".join(proc.stdout.splitlines()[-30:])
        raise RuntimeError(f"command failed ({proc.returncode}): {cmd[0]}\n{tail}")
    return proc.stdout
