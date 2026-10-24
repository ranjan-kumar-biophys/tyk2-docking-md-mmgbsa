"""Step 4: build solvated protein-ligand systems and convert them to GROMACS.

ff14SB (protein) + GAFF2/AM1-BCC (ligand) + TIP3P + Joung-Cheatham ions, truncated
octahedron, neutralised and brought to the salt concentration in config.yaml. tleap
builds the system, ParmEd writes GROMACS files, and this script adds heavy-atom
position restraints (POSRES) for equilibration plus an index file whose group order is
fixed: 0 System, 1 Protein, 2 LIG, 3 Protein_LIG, 4 Water_and_ions.

--pose reference builds from the benchmark pose into work/systems_reference (control
runs for ligands that docking placed wrongly).
"""
import argparse
import json
import re

import numpy as np
import parmed as pmd
from rdkit import Chem

from common import WORK, load_config, ligand_names, run, stage_dir

SOLVENT = {"WAT", "Na+", "Cl-"}
LEAP_HEAD = """source leaprc.protein.ff14SB
source leaprc.gaff2
source leaprc.water.tip3p
loadamberparams lig.frcmod
LIG = loadmol2 lig_pose.mol2
REC = loadpdb {receptor}
complex = combine {{REC LIG}}
solvateoct complex TIP3PBOX {padding}
"""


def write_pose_mol2(template, pose, ref, out):
    """Copy pose coordinates onto the parameterised mol2 (same atom order as ref)."""
    lines = template.read_text().splitlines()
    i0 = lines.index("@<TRIPOS>ATOM") + 1
    i1 = next(i for i in range(i0, len(lines)) if lines[i].startswith("@<TRIPOS>"))
    mol2_xyz = np.array([[float(x) for x in l.split()[2:5]] for l in lines[i0:i1]])
    if not np.allclose(mol2_xyz, ref.GetConformer().GetPositions(), atol=0.01):
        raise RuntimeError(f"{template}: atom order differs from the reference SDF")
    for k, xyz in enumerate(pose.GetConformer().GetPositions()):
        f = lines[i0 + k].split()
        lines[i0 + k] = (f"{int(f[0]):7d} {f[1]:<8s}{xyz[0]:10.4f}{xyz[1]:10.4f}{xyz[2]:10.4f} "
                         f"{f[5]:<6s}{f[6]:>4s} {f[7]:<6s}{float(f[8]):10.6f}")
    out.write_text("\n".join(lines) + "\n")


def tleap(d, script, tag):
    (d / f"{tag}.in").write_text(script)
    log = run(["tleap", "-f", f"{tag}.in"], cwd=d, log=d / f"{tag}.log")
    bad = [l for l in log.splitlines() if "FATAL" in l or "Created a new atom" in l
           or l.startswith("Error")]
    if bad:
        raise RuntimeError(f"tleap problems in {d / tag}.log:\n" + "\n".join(bad[:10]))
    return log


def add_posres(top_path, heavy_by_mol, fc=1000):
    """Append an #ifdef POSRES block to the end of each named moleculetype."""
    out, current = [], None

    def flush():
        if current in heavy_by_mol:
            out.extend(["", "#ifdef POSRES", "[ position_restraints ]"])
            out.extend(f"{i:6d}     1  {fc}  {fc}  {fc}" for i in heavy_by_mol[current])
            out.append("#endif")

    lines = top_path.read_text().splitlines()
    for i, line in enumerate(lines):
        if re.match(r"\s*\[\s*(moleculetype|system)\s*\]", line):
            flush()
            current = None
            if "moleculetype" in line:
                nxt = next(l for l in lines[i + 1:] if l.strip() and not l.lstrip().startswith(";"))
                current = nxt.split()[0]
        out.append(line)
    top_path.write_text("\n".join(out) + "\n")


def molecule_heavy_atoms(top_path):
    """{moleculetype name: [heavy atom indices]} for the protein and the ligand."""
    result, current, section, atoms = {}, None, None, []
    for line in top_path.read_text().splitlines() + ["[ moleculetype ]"]:
        s = line.split(";")[0].strip()
        m = re.match(r"\[\s*(\w+)\s*\]", s)
        if m:
            if section == "atoms" and current:
                resnames = {a[3] for a in atoms}
                if "LIG" in resnames or not resnames & SOLVENT and len(atoms) > 100:
                    result[current] = [int(a[0]) for a in atoms if float(a[7]) > 1.5]
            section = m.group(1)
            if section == "moleculetype":
                current, atoms = None, []
            continue
        if not s or s.startswith("#"):  # preprocessor lines, e.g. #ifdef FLEXIBLE in WAT
            continue
        if section == "moleculetype" and current is None:
            current = s.split()[0]
        elif section == "atoms":
            atoms.append(s.split())
    return result


def write_index(struct, path):
    groups = {"System": [], "Protein": [], "LIG": [], "Protein_LIG": [], "Water_and_ions": []}
    for i, atom in enumerate(struct.atoms, start=1):
        res = atom.residue.name
        groups["System"].append(i)
        if res in SOLVENT:
            groups["Water_and_ions"].append(i)
        else:
            groups["LIG" if res == "LIG" else "Protein"].append(i)
            groups["Protein_LIG"].append(i)
    with open(path, "w") as fh:
        for name, idx in groups.items():
            fh.write(f"[ {name} ]\n")
            for k in range(0, len(idx), 15):
                fh.write(" ".join(f"{x:6d}" for x in idx[k:k + 15]) + "\n")
    return {k: len(v) for k, v in groups.items()}


def build(name, cfg, pose):
    d = stage_dir("systems", pose) / name
    d.mkdir(parents=True, exist_ok=True)
    lig_dir = WORK / "ligands" / name
    ref = Chem.MolFromMolFile(str(lig_dir / "ref.sdf"), removeHs=False)
    pose_file = WORK / "docking" / name / "top_pose.sdf" if pose == "docked" else lig_dir / "ref.sdf"
    pose_mol = Chem.MolFromMolFile(str(pose_file), removeHs=False)
    write_pose_mol2(lig_dir / "lig.mol2", pose_mol, ref, d / "lig_pose.mol2")
    (d / "lig.frcmod").write_text((lig_dir / "lig.frcmod").read_text())

    head = LEAP_HEAD.format(receptor=WORK / "receptor" / "receptor_tleap.pdb",
                            padding=cfg["build"]["water_padding_A"])
    log = tleap(d, head + "charge complex\nquit\n", "leap_probe")
    n_wat = int(re.search(r"Added (\d+) residues", log).group(1))
    charge = round(float(re.search(r"Total unperturbed charge:\s*(-?[\d.]+)", log).group(1)))
    n_salt = round(cfg["build"]["salt_M"] * n_wat / 55.5)
    n_na, n_cl = n_salt + max(0, -charge), n_salt + max(0, charge)
    tleap(d, head + f"addionsrand complex Na+ {n_na} Cl- {n_cl}\n"
                    "saveamberparm complex complex.prmtop complex.inpcrd\n"
                    "savepdb complex complex_leap.pdb\nquit\n", "leap_build")

    struct = pmd.load_file(str(d / "complex.prmtop"), xyz=str(d / "complex.inpcrd"))
    struct.save(str(d / "complex.top"), format="gromacs", overwrite=True)
    struct.save(str(d / "complex.gro"), overwrite=True)
    heavy = molecule_heavy_atoms(d / "complex.top")
    if len(heavy) != 2:
        raise RuntimeError(f"{name}: expected protein + ligand moleculetypes, found {list(heavy)}")
    add_posres(d / "complex.top", heavy)
    counts = write_index(struct, d / "index.ndx")

    info = {"start_pose": pose, "solute_charge": charge,
            "waters": n_wat, "Na+": n_na, "Cl-": n_cl, "atoms": len(struct.atoms),
            "groups": counts}
    (d / "system_info.json").write_text(json.dumps(info, indent=2))
    return info


def main():
    cfg = load_config()
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--ligands", nargs="+", help="subset of ligands (default: config)")
    p.add_argument("--pose", choices=["docked", "reference"], default=cfg["build"]["start_pose"])
    args = p.parse_args()
    for name in args.ligands or ligand_names(cfg):
        info = build(name, cfg, args.pose)
        print(f"{name:12s} {info['atoms']} atoms, {info['waters']} waters, "
              f"{info['Na+']} Na+ / {info['Cl-']} Cl-, solute charge {info['solute_charge']}")


if __name__ == "__main__":
    main()
