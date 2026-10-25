"""Step 12 (part 2): render the MD movie with PyMOL and encode it with ffmpeg.

  pymol -cq scripts/12_movie_render.py -- lig_jmc_23 1      (run from the repo root)

Loads work/movie/<ligand>_rep<N>.pdb (from part 1), smooths the frames, assigns secondary
structure, shows the ligand and key pocket residues over the trajectory, turns the camera
slowly, ray-traces each frame and writes figures/movie_<ligand>.mp4 (for posts) and
figures/movie_<ligand>.gif (plays in the README).
"""
import glob
import os
import shutil
import subprocess
import sys

from pymol import cmd

ROOT = os.getcwd()
name, rep = sys.argv[1], sys.argv[2]
pdb = os.path.join(ROOT, "work", "movie", f"{name}_rep{rep}.pdb")
frames = os.path.join(ROOT, "work", "movie", f"frames_{name}")
FIG = os.path.join(ROOT, "figures")
KEY = "903+911+928+978+979+980+981+1030+1041"
ORANGE, PROTEIN_GREY = "0xEB6834", "0xCFD6DE"
W, H, FPS = 1600, 1200, 25

cmd.load(pdb, "md")
cmd.alter("md", "resi=str(int(resi) + 888)")          # MDAnalysis numbers 1..N; 4GIH starts at 889
cmd.sort("md")
cmd.smooth("md", 2, 3)                                 # gentle smoothing for a calm movie
cmd.dss("md")                                          # MD frames carry no HELIX/SHEET records
cmd.remove("md and hydro and not (neighbor elem N+O)")
cmd.hide("everything")
cmd.bg_color("white")
for k, v in {"ray_trace_mode": 1, "ray_trace_gain": 0.6, "antialias": 2, "ambient": 0.35,
             "specular": 0.15, "ray_shadows": 0, "depth_cue": 0, "cartoon_fancy_helices": 1,
             "cartoon_transparency": 0.3, "stick_radius": 0.14, "dash_color": "0x0B0B0B",
             "dash_radius": 0.07, "dash_gap": 0.22}.items():
    cmd.set(k, v)
cmd.show("cartoon", "md and polymer")
cmd.color(PROTEIN_GREY, "md and polymer")
cmd.set("cartoon_color", PROTEIN_GREY, "md")
cmd.select("lig", "md and resn LIG")
cmd.show("sticks", "lig")
cmd.set("stick_radius", 0.26, "lig")
cmd.color(ORANGE, "lig and elem C")
cmd.select("key", f"md and polymer and resi {KEY} and not name N+C+O")
cmd.show("sticks", "key or (md and resi 979-981 and name N+C+O+CA)")
cmd.color("0x3A3A38", "(key or (md and resi 979-981)) and elem C")
cmd.util.cnc("md")

n = cmd.count_states("md")
cmd.orient("lig")
cmd.zoom("lig", 6)
cmd.turn("x", 60)                                      # same view as the hero render
cmd.mset(f"1 -{n}")
cmd.mview("store", 1)
cmd.turn("y", 35)
cmd.mview("store", n)
cmd.mview("interpolate")

shutil.rmtree(frames, ignore_errors=True)
os.makedirs(frames)
cmd.viewport(W, H)
cmd.set("ray_trace_frames", 1)
cmd.mpng(os.path.join(frames, "f"), width=W, height=H)
pngs = sorted(glob.glob(os.path.join(frames, "f*.png")))
print(f"rendered {len(pngs)} frames")

mp4, gif = os.path.join(FIG, f"movie_{name}.mp4"), os.path.join(FIG, f"movie_{name}.gif")
pattern = os.path.join(frames, "f%04d.png")
subprocess.run(["ffmpeg", "-y", "-loglevel", "error", "-framerate", str(FPS), "-i", pattern,
                "-c:v", "libx264", "-pix_fmt", "yuv420p", "-crf", "18", mp4], check=True)
subprocess.run(["ffmpeg", "-y", "-loglevel", "error", "-framerate", str(FPS), "-i", pattern, "-vf",
                "fps=15,scale=720:-1:flags=lanczos,split[a][b];[a]palettegen=max_colors=128[p];[b][p]paletteuse",
                gif], check=True)
print("wrote", mp4, "and", gif)
