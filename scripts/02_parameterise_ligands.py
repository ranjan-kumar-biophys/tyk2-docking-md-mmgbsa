"""Step 2: GAFF2 + AM1-BCC parameters for every ligand (antechamber / parmchk2).

Charges come from the benchmark 3D pose with its explicit hydrogens, so protonation
states match the benchmark. Atom order in lig.mol2 = order in the benchmark SDF; later
steps rely on that to copy docked coordinates onto the parameterised molecule.

antechamber's rounded charges can miss the formal charge by ~0.001 e, which makes
GROMACS refuse the system (net charge under PME), so the residual is spread back
over the atoms in exact 1e-6 e steps.
"""
import os
from concurrent.futures import ThreadPoolExecutor

from rdkit import Chem

from common import BENCH, WORK, load_config, ligand_names, run


def fix_total_charge(mol2, net_charge):
    """Make the mol2 partial charges sum exactly to net_charge; return the correction."""
    lines = mol2.read_text().splitlines()
    i0 = lines.index("@<TRIPOS>ATOM") + 1
    i1 = next(i for i in range(i0, len(lines)) if lines[i].startswith("@<TRIPOS>"))
    micro = [round(float(l.split()[8]) * 1e6) for l in lines[i0:i1]]
    diff = net_charge * 1_000_000 - sum(micro)
    n = len(micro)
    micro = [q + diff // n + (1 if k < diff % n else 0) for k, q in enumerate(micro)]
    for k, i in enumerate(range(i0, i1)):
        head = lines[i].rsplit(None, 1)[0]
        lines[i] = f"{head} {micro[k] / 1e6:10.6f}"
    mol2.write_text("\n".join(lines) + "\n")
    return diff / 1e6


def parameterise(name, mol):
    d = WORK / "ligands" / name
    d.mkdir(parents=True, exist_ok=True)
    charge = Chem.GetFormalCharge(mol)
    if not (d / "lig.frcmod").exists():
        with Chem.SDWriter(str(d / "ref.sdf")) as w:
            w.write(mol)
        run(["antechamber", "-i", "ref.sdf", "-fi", "mdl", "-o", "lig.mol2", "-fo", "mol2",
             "-c", "bcc", "-at", "gaff2", "-nc", charge, "-rn", "LIG", "-s", "2", "-pf", "y"],
            cwd=d, log=d / "antechamber.log")
        run(["parmchk2", "-i", "lig.mol2", "-f", "mol2", "-o", "lig.frcmod", "-s", "gaff2"], cwd=d)
    fixed = fix_total_charge(d / "lig.mol2", charge)
    if "ATTN" in (d / "lig.frcmod").read_text():
        return name, "WARNING: parmchk2 flagged missing parameters (see lig.frcmod)"
    return name, f"ok (net charge {charge}, rounding fix {fixed:+.6f} e)"


def main():
    cfg = load_config()
    wanted = set(ligand_names(cfg))
    mols = {m.GetProp("_Name"): m
            for m in Chem.SDMolSupplier(str(BENCH / "ligands.sdf"), removeHs=False)}
    with ThreadPoolExecutor(max_workers=os.cpu_count()) as pool:
        for name, status in pool.map(lambda n: parameterise(n, mols[n]), sorted(wanted)):
            print(f"{name:12s} {status}")


if __name__ == "__main__":
    main()
