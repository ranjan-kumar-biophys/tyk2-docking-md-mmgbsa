"""Step 3: dock every ligand into TYK2 with AutoDock Vina and score pose accuracy.

To be fair, docking starts from a fresh RDKit conformer (ETKDG + MMFF), not from the
benchmark pose. The box is centred on the reference ligands. Poses are written back in
the benchmark atom order so step 4 can drop them straight onto the GAFF2 topology.
"""
import numpy as np
import pandas as pd
from meeko import MoleculePreparation, PDBQTMolecule, PDBQTWriterLegacy, RDKitMolCreate
from rdkit import Chem
from rdkit.Chem import AllChem, rdMolAlign
from vina import Vina

from common import BENCH, RESULTS, WORK, load_config, ligand_names


def heavy_coords(mol):
    pos = mol.GetConformer().GetPositions()
    return pos[[a.GetIdx() for a in mol.GetAtoms() if a.GetAtomicNum() > 1]]


def dock_one(name, ref, receptor, center, cfg):
    d = WORK / "docking" / name
    d.mkdir(parents=True, exist_ok=True)

    start = Chem.Mol(ref)
    start.RemoveAllConformers()
    AllChem.EmbedMolecule(start, randomSeed=cfg["seed"])
    AllChem.MMFFOptimizeMolecule(start)
    setup = MoleculePreparation().prepare(start)[0]
    pdbqt, ok, err = PDBQTWriterLegacy.write_string(setup)
    if not ok:
        raise RuntimeError(f"{name}: meeko failed: {err}")

    v = Vina(sf_name="vina", seed=cfg["seed"], verbosity=0)
    v.set_receptor(str(receptor))
    v.set_ligand_from_string(pdbqt)
    v.compute_vina_maps(center=list(center), box_size=[cfg["box_size_A"]] * 3)
    v.dock(exhaustiveness=cfg["exhaustiveness"], n_poses=cfg["n_poses"])
    n = cfg["n_poses"]
    out_pdbqt = v.poses(n_poses=n)
    scores = v.energies(n_poses=n)[:, 0]
    (d / "poses.pdbqt").write_text(out_pdbqt)

    docked = RDKitMolCreate.from_pdbqt_mol(PDBQTMolecule(out_pdbqt, skip_typing=True))[0]
    match = docked.GetSubstructMatch(ref)
    if len(match) != ref.GetNumAtoms():
        raise RuntimeError(f"{name}: docked molecule does not map onto the reference")
    docked = Chem.RenumberAtoms(docked, list(match))  # now in benchmark atom order

    ref_heavy = Chem.RemoveHs(ref)
    rows = []
    with Chem.SDWriter(str(d / "poses.sdf")) as w:
        for i, conf in enumerate(docked.GetConformers()):
            pose = Chem.Mol(docked, confId=conf.GetId())
            rmsd = rdMolAlign.CalcRMS(Chem.RemoveHs(pose), ref_heavy)
            pose.SetProp("_Name", f"{name}_pose{i + 1}")
            pose.SetProp("vina_score", f"{scores[i]:.3f}")
            pose.SetProp("rmsd_to_ref_A", f"{rmsd:.2f}")
            w.write(pose)
            rows.append((scores[i], rmsd))
            if i == 0:
                with Chem.SDWriter(str(d / "top_pose.sdf")) as top:
                    top.write(pose)

    best = int(np.argmin([r for _, r in rows]))
    return {"ligand": name, "vina_score": rows[0][0], "top_pose_rmsd_A": rows[0][1],
            "best_pose_rmsd_A": rows[best][1], "best_pose_rank": best + 1}


def main():
    cfg = load_config()
    names = ligand_names(cfg)
    refs = {m.GetProp("_Name"): m
            for m in Chem.SDMolSupplier(str(BENCH / "ligands.sdf"), removeHs=False)}
    center = np.vstack([heavy_coords(m) for m in refs.values()]).mean(axis=0)
    receptor = WORK / "receptor" / "receptor_dock.pdbqt"
    print("box centre (A):", np.round(center, 2), "size:", cfg["docking"]["box_size_A"])

    rows = []
    for name in names:
        row = dock_one(name, refs[name], receptor, center, cfg["docking"])
        print(f"{name:12s} vina {row['vina_score']:6.2f}  top-pose RMSD {row['top_pose_rmsd_A']:5.2f} A")
        rows.append(row)

    RESULTS.mkdir(exist_ok=True)
    df = pd.DataFrame(rows)
    df.to_csv(RESULTS / "docking.csv", index=False, float_format="%.3f")
    ok = (df["top_pose_rmsd_A"] < 2.0).sum()
    print(f"top pose within 2 A of the benchmark pose: {ok}/{len(df)}")


if __name__ == "__main__":
    main()
