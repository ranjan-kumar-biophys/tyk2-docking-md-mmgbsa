"""Step 9: per-replica MD analysis (MDAnalysis) for the figures in step 11.

For every production replica: protein backbone RMSD, symmetry-aware ligand RMSD to the
benchmark (crystal-like) pose, Cα RMSF, radius of gyration, ligand–residue H-bond and contact
occupancies, and strided Cα coordinates for a shared PCA. Frames are superposed on the
benchmark protein's Cα atoms, so ligand RMSD is measured in the crystal frame: it shows
whether MD keeps, repairs or loses a docked pose. Only protein + ligand are held in memory.
Writes work/analysis/<pose>/<ligand>/rep<N>.npz (+ .json with occupancies).
"""
import argparse
import json
import warnings

import MDAnalysis as mda
import numpy as np
from MDAnalysis.analysis import align, rms
from MDAnalysis.analysis.base import AnalysisFromFunction
from MDAnalysis.analysis.hydrogenbonds import HydrogenBondAnalysis
from MDAnalysis.coordinates.memory import MemoryReader
from MDAnalysis.lib.distances import capped_distance
from rdkit import Chem

from common import BENCH, WORK, load_config, ligand_names, stage_dir

warnings.filterwarnings("ignore", category=UserWarning)
CONTACT_A = 4.0      # heavy-atom contact cutoff
PCA_STRIDE = 5       # every 5th frame (50 ps) for the shared PCA


def load_solute(sysdir, xtc):
    """Protein + ligand of one replica, all frames, in memory (water and ions dropped)."""
    u = mda.Universe(str(sysdir / "complex.prmtop"), str(xtc))
    solute = u.select_atoms("protein or resname LIG or resname ACE NME")
    coords = AnalysisFromFunction(lambda ag: ag.positions.copy(), solute).run().results.timeseries
    m = mda.Merge(solute)
    m.load_new(coords, format=MemoryReader, dt=u.trajectory.dt)
    return m


def crystal_ca():
    """Cα atoms of the benchmark (crystal-frame) protein."""
    return mda.Universe(str(BENCH / "protein.pdb")).select_atoms("name CA")


def analyse_replica(name, pose, rep, cfg):
    sysdir = stage_dir("systems", pose) / name
    mddir = stage_dir("md", pose) / name / f"rep{rep}"
    out = WORK / "analysis" / pose / name
    out.mkdir(parents=True, exist_ok=True)
    m = load_solute(sysdir, mddir / "prod_center.xtc")

    ca = m.select_atoms("name CA")
    ref_ca = crystal_ca()
    if ca.n_atoms != ref_ca.n_atoms:
        raise RuntimeError(f"{name}: {ca.n_atoms} Cα in MD vs {ref_ca.n_atoms} in the benchmark")
    # superpose every frame on the crystal Cα (in place), so all positions are in that frame
    align.AlignTraj(m, ref_ca.universe, select=("name CA", "name CA"), in_memory=True).run()

    lig = m.select_atoms("resname LIG and not element H")
    ref_mol = Chem.MolFromMolFile(str(WORK / "ligands" / name / "ref.sdf"))  # heavy atoms only
    ref_lig = ref_mol.GetConformer().GetPositions()
    if len(ref_lig) != lig.n_atoms:
        raise RuntimeError(f"{name}: ligand heavy-atom count mismatch")
    # symmetry-equivalent atom orders (e.g. a flipped 2,6-dichlorophenyl ring is the same pose)
    perms = [list(p) for p in ref_mol.GetSubstructMatches(ref_mol, uniquify=False, maxMatches=1000)]

    bb = m.select_atoms("backbone")
    m.trajectory[0]
    bb0 = bb.positions.copy()
    t, rmsd_bb, rmsd_lig, rg = [], [], [], []
    protein = m.select_atoms("protein")
    heavy_prot = m.select_atoms("protein and not element H")
    contact_frames = {}
    ca_stride = []
    for ts in m.trajectory:
        t.append(ts.time / 1000.0)  # ns
        rmsd_bb.append(rms.rmsd(bb.positions, bb0, superposition=False))
        pos = lig.positions
        rmsd_lig.append(min(np.sqrt(((pos[p] - ref_lig) ** 2).sum(axis=1).mean()) for p in perms))
        rg.append(protein.radius_of_gyration())
        pairs = capped_distance(lig.positions, heavy_prot.positions, CONTACT_A,
                                box=None, return_distances=False)
        for resid in {heavy_prot[j].resid for _, j in pairs}:
            contact_frames[resid] = contact_frames.get(resid, 0) + 1
        if ts.frame % PCA_STRIDE == 0:
            ca_stride.append(ca.positions.copy())
    n = len(t)

    rmsf = rms.RMSF(ca).run().results.rmsf  # after the crystal superposition above

    hb = HydrogenBondAnalysis(m, between=["resname LIG", "protein"], d_a_cutoff=3.5,
                              d_h_a_angle_cutoff=135)
    hb.hydrogens_sel = hb.guess_hydrogens("resname LIG or protein")
    hb.acceptors_sel = hb.guess_acceptors("resname LIG or protein")
    hb.run()
    hb_frames = {}
    for frame, d, h, a, *_ in hb.results.hbonds:
        d, a = int(d), int(a)
        partner = m.atoms[d] if m.atoms[d].resname != "LIG" else m.atoms[a]
        hb_frames.setdefault(partner.resid, set()).add(int(frame))

    resname = {r.resid: f"{r.resname}{r.resid}" for r in m.residues}
    # resids here are 1..N in file order; map back to 4GIH numbering (ACE = 889)
    offset = m.residues[0].resid
    to_4gih = lambda r: f"{resname[r][:3]} {r - offset + 889}"
    occupancy = {
        "contacts": {to_4gih(r): c / n for r, c in sorted(contact_frames.items())},
        "hbonds": {to_4gih(r): len(f) / n for r, f in sorted(hb_frames.items())},
    }
    np.savez_compressed(out / f"rep{rep}.npz", time_ns=np.array(t), rmsd_backbone=np.array(rmsd_bb),
                        rmsd_ligand_to_benchmark=np.array(rmsd_lig), rg=np.array(rg), rmsf_ca=rmsf,
                        ca_resid_4gih=np.array([r - offset + 889 for r in ca.resids]),
                        ca_coords=np.array(ca_stride, dtype=np.float32))
    (out / f"rep{rep}.json").write_text(json.dumps(occupancy, indent=1))
    return {"frames": n, "final_lig_rmsd": rmsd_lig[-1], "mean_bb_rmsd": float(np.mean(rmsd_bb)),
            "top_hbond": max(occupancy["hbonds"].items(), key=lambda kv: kv[1], default=("none", 0))}


def main():
    cfg = load_config()
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--ligands", nargs="+")
    p.add_argument("--replicas", nargs="+", type=int)
    p.add_argument("--pose", choices=["docked", "reference"], default=cfg["build"]["start_pose"])
    args = p.parse_args()
    replicas = args.replicas or list(range(1, cfg["md"]["replicas"] + 1))
    for name in args.ligands or ligand_names(cfg):
        for rep in replicas:
            if not (stage_dir("md", args.pose) / name / f"rep{rep}" / "prod_center.xtc").exists():
                print(f"{name} rep{rep}: no trajectory, skipping")
                continue
            if (WORK / "analysis" / args.pose / name / f"rep{rep}.npz").exists():
                continue
            r = analyse_replica(name, args.pose, rep, cfg)
            print(f"{name:12s} {args.pose:9s} rep{rep}: {r['frames']} frames, backbone RMSD "
                  f"{r['mean_bb_rmsd']:.2f} Å, final ligand RMSD to benchmark {r['final_lig_rmsd']:.2f} Å, "
                  f"top H-bond {r['top_hbond'][0]} ({r['top_hbond'][1]:.0%})")


if __name__ == "__main__":
    main()
