"""Step 10: publication-quality structure renders with PyMOL (headless, ray-traced).

  pymol -cq scripts/10_render.py            (run from the repo root)

Writes to figures/:
  render_hero.png            strongest binder (lig_jmc_23) in the TYK2 pocket, hinge H-bonds shown
  render_pose_<ligand>.png   docked pose (blue) over the benchmark pose (orange), one per ligand
                             listed in POSE_LIGANDS: a docking success and a docking failure
White background, outlined ray tracing, 2400 x 1800 px at 300 dpi.
"""
import json
import os

from pymol import cmd

ROOT = os.getcwd()
BENCH = os.path.join(ROOT, "inputs", "benchmark")
WORK = os.path.join(ROOT, "work")
FIG = os.path.join(ROOT, "figures")
POSE_LIGANDS = ["lig_jmc_23", "lig_ejm_46"]   # 0.33 Å success vs ~7.5 Å flip in Project 1
ORANGE, BLUE, PROTEIN_GREY = "0xEB6834", "0x2A78D6", "0xCFD6DE"
# Residues that dominate binding in the MM/GBSA decomposition: P-loop Leu903, Val911, Ala928,
# gatekeeper Met978, hinge Glu979-Tyr980-Val981, Leu1030, DFG Asp1041.
KEY = "903+911+928+978+979+980+981+1030+1041"
LABELS = {"978": "Met978 (gatekeeper)", "981": "Val981 (hinge)", "1041": "Asp1041 (DFG)"}


def style_scene():
    cmd.bg_color("white")
    cmd.set("ray_opaque_background", 1)
    cmd.set("ray_trace_mode", 1)          # shaded + black outline
    cmd.set("ray_trace_gain", 0.6)
    cmd.set("antialias", 2)
    cmd.set("ambient_occlusion_mode", 1)
    cmd.set("ambient", 0.35)
    cmd.set("specular", 0.15)
    cmd.set("ray_shadows", 0)
    cmd.set("depth_cue", 0)
    cmd.set("cartoon_fancy_helices", 1)
    cmd.set("cartoon_transparency", 0.3)
    cmd.set("stick_radius", 0.14)
    cmd.set("surface_quality", 1)
    cmd.set("transparency", 0.88)
    cmd.set("dash_color", "0x0B0B0B")
    cmd.set("dash_gap", 0.22)
    cmd.set("dash_radius", 0.07)
    cmd.set("label_size", -1.05)          # world units: scales with the scene
    cmd.set("label_color", "0x0B0B0B")
    cmd.set("label_font_id", 7)
    cmd.set("label_outline_color", "white")
    cmd.set("float_labels", 1)
    cmd.set("label_position", (0, 1.2, 6))   # lift labels toward the camera so they draw in front


def load_protein():
    cmd.load(os.path.join(BENCH, "protein.pdb"), "prot")
    cmd.remove("prot and (resn HOH or hydro)")
    cmd.hide("everything", "prot")
    cmd.show("cartoon", "prot")
    cmd.color(PROTEIN_GREY, "prot")
    cmd.set("cartoon_color", PROTEIN_GREY, "prot")   # stays grey when residue atoms are recoloured


def style_ligand(obj, color):
    cmd.remove(f"{obj} and hydro and not (neighbor elem N+O)")   # keep polar H only
    cmd.show("sticks", obj)
    cmd.set("stick_radius", 0.26, obj)
    cmd.color(color, f"{obj} and elem C")
    cmd.util.cnc(obj)


def show_pocket(ligand_obj, label=True):
    cmd.select("key", f"prot and resi {KEY} and not name N+C+O")
    cmd.show("sticks", "key and not hydro")
    cmd.show("sticks", "prot and resi 979-981 and name N+C+O+CA")   # hinge backbone carries the H-bonds
    cmd.color("0x3A3A38", "key and elem C")
    cmd.color("0x3A3A38", "prot and resi 979-981 and elem C")
    cmd.util.cnc("prot and resi " + KEY)
    # no pocket surface: in ray tracing a translucent surface washes out the residues behind it
    cmd.distance("hbonds", ligand_obj, "prot and resi 978-982", mode=2, cutoff=3.4)
    cmd.hide("labels", "hbonds")
    if label:
        for resi, text in LABELS.items():
            cmd.label(f"prot and resi {resi} and name CB", f'"{text}"')


# Camera: orient on the ligand, then turn. TURN_X/TURN_Y and WIDTH can be overridden from the
# environment to try angles with quick small previews before the full-size render.
TURN_X = float(os.environ.get("TURN_X", 60))   # +60: view into the cleft, nothing in front of the ligand
TURN_Y = float(os.environ.get("TURN_Y", 0))
WIDTH = int(os.environ.get("WIDTH", 2400))


def render(path, ligand_obj):
    cmd.orient(ligand_obj)
    cmd.zoom(ligand_obj, 5)
    cmd.turn("x", TURN_X)
    cmd.turn("y", TURN_Y)
    cmd.ray(WIDTH, WIDTH * 3 // 4)
    cmd.png(path, dpi=300)
    print("wrote", path)


def hero():
    cmd.reinitialize()
    style_scene()
    load_protein()
    cmd.load(os.path.join(WORK, "ligands", "lig_jmc_23", "ref.sdf"), "lig")
    style_ligand("lig", ORANGE)
    show_pocket("lig")
    render(os.path.join(FIG, "render_hero.png"), "lig")


def pose_overlay(name):
    cmd.reinitialize()
    style_scene()
    load_protein()
    cmd.load(os.path.join(WORK, "ligands", name, "ref.sdf"), "benchmark")
    cmd.load(os.path.join(WORK, "docking", name, "top_pose.sdf"), "docked")
    for obj, color in (("benchmark", ORANGE), ("docked", BLUE)):
        style_ligand(obj, color)
        cmd.remove(f"{obj} and hydro")
    show_pocket("benchmark", label=False)
    rmsd = cmd.rms_cur("docked", "benchmark", matchmaker=-1)  # same atom order, no fitting
    render(os.path.join(FIG, f"render_pose_{name}.png"), "benchmark or docked")
    return rmsd


os.makedirs(FIG, exist_ok=True)
hero()
if os.environ.get("HERO_ONLY"):
    raise SystemExit
rmsds = {name: pose_overlay(name) for name in POSE_LIGANDS}
with open(os.path.join(FIG, "render_pose_rmsd.json"), "w") as fh:
    json.dump(rmsds, fh, indent=1)
print("pose RMSD (Å, docked vs benchmark):", {k: round(v, 2) for k, v in rmsds.items()})
