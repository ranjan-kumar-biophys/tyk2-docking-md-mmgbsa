"""Step 1: prepare the TYK2 receptor for docking (Vina) and for system building (tleap).

The benchmark PDB (4GIH, already protonated and capped) is used as-is, with two fixes:
  * HIS -> HID/HIE/HIP from the hydrogens that are present, so tleap keeps the same tautomers
  * hydrogens are stripped for tleap (it rebuilds them with ff14SB names)
The two crystal waters (~3 A from the ligand) are kept for MD but not used in docking.
"""
from collections import defaultdict

from common import BENCH, WORK, run


def residue_key(line):
    return line[21], int(line[22:26])


def is_hydrogen(line):
    element = line[76:78].strip()
    return element == "H" if element else line[12:16].strip().startswith("H")


def main():
    out = WORK / "receptor"
    out.mkdir(parents=True, exist_ok=True)

    lines = [l for l in (BENCH / "protein.pdb").read_text().splitlines()
             if l.startswith(("ATOM", "HETATM"))]

    his_h = defaultdict(set)
    for l in lines:
        if l[17:20] == "HIS":
            his_h[residue_key(l)].add(l[12:16].strip())

    def his_state(key):
        names = his_h[key]
        if {"HD1", "HE2"} <= names:
            return "HIP"
        return "HID" if "HD1" in names else "HIE"

    protein, waters = [], []
    for l in lines:
        (waters if l[17:20] == "HOH" else protein).append(l)

    # Docking receptor: protein only, original names, all hydrogens.
    dock_pdb = out / "receptor_dock.pdb"
    dock_pdb.write_text("\n".join(protein) + "\nTER\nEND\n")
    run(["obabel", "-ipdb", dock_pdb, "-xr", "-opdbqt", "-O", out / "receptor_dock.pdbqt"])

    # tleap receptor: Amber residue names, heavy atoms only, waters as separate WAT residues.
    amber = []
    for l in protein:
        if is_hydrogen(l):
            continue
        if l[17:20] == "HIS":
            l = l[:17] + his_state(residue_key(l)) + l[20:]
        amber.append(l)
    amber.append("TER")
    for l in waters:
        if not is_hydrogen(l):
            amber += ["HETATM" + l[6:17] + "WAT" + l[20:], "TER"]
    (out / "receptor_tleap.pdb").write_text("\n".join(amber) + "\nEND\n")

    states = defaultdict(list)
    for key in his_h:
        states[his_state(key)].append(key[1])
    print("HIS states:", {k: sorted(v) for k, v in states.items()})
    print(f"wrote {out}/receptor_dock.pdbqt and receptor_tleap.pdb ({len(waters) // 3} crystal waters)")


if __name__ == "__main__":
    main()
