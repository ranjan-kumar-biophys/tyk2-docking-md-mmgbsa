"""Step 12 (part 1): write MD frames for the movie (run in the pipeline env).

  python scripts/12_movie.py lig_jmc_23 1                  -> work/movie/lig_jmc_23_rep1.pdb
  pymol -cq scripts/12_movie_render.py -- lig_jmc_23 1     (part 2: PDB -> PNGs -> MP4 + GIF)

Takes every STRIDE-th frame of one production replica, keeps protein + ligand, and
superposes them on the crystal Cα so the camera can stay still while the molecule moves.
"""
import sys

import MDAnalysis as mda
from MDAnalysis.analysis import align
from MDAnalysis.analysis.base import AnalysisFromFunction
from MDAnalysis.coordinates.memory import MemoryReader

from common import BENCH, WORK, stage_dir

STRIDE = 5   # 1001 frames x 10 ps -> 201 frames, 8 s at 25 fps


def main():
    name, rep = sys.argv[1], int(sys.argv[2])
    pose = sys.argv[3] if len(sys.argv) > 3 else "docked"
    sysdir, mddir = stage_dir("systems", pose) / name, stage_dir("md", pose) / name / f"rep{rep}"
    u = mda.Universe(str(sysdir / "complex.prmtop"), str(mddir / "prod_center.xtc"))
    solute = u.select_atoms("protein or resname LIG or resname ACE NME")
    coords = AnalysisFromFunction(lambda ag: ag.positions.copy(), solute).run(step=STRIDE).results.timeseries
    m = mda.Merge(solute)
    m.load_new(coords, format=MemoryReader)
    ref = mda.Universe(str(BENCH / "protein.pdb"))
    align.AlignTraj(m, ref, select=("name CA", "name CA"), in_memory=True).run()

    out = WORK / "movie" / f"{name}_rep{rep}.pdb"
    out.parent.mkdir(parents=True, exist_ok=True)
    with mda.Writer(str(out), m.atoms.n_atoms, multiframe=True) as w:
        for _ in m.trajectory:
            w.write(m.atoms)
    print(f"wrote {out} ({m.trajectory.n_frames} frames, {m.atoms.n_atoms} atoms)")


if __name__ == "__main__":
    main()
