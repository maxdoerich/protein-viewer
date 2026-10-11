"""
Gas vesicle animation  --  GvpA monomer -> shell lattice -> GvpC cage
======================================================================

Builds the complete scene (structures, materials, lights, camera, animation,
render settings) for Blender 5.1 + Molecular Nodes 4.5.x.

HOW TO RUN
  1. Keep the folder "data" next to this script.
  2. Blender -> "Scripting" workspace -> Text editor -> Open -> this file -> Run Script (Alt+P / Play button).
     (Takes ~1 min. The current scene is cleared first.)
  3. Press Space in the viewport to preview, or Render -> Render Animation (Ctrl+F12).

All knobs you might want to change are in the SETTINGS block below.

Data sources
  PDB 7R1C (B. megaterium GvpA2, Huber et al. 2023, Cell), EMD-14238 (local wall map),
  full-vesicle pseudo-atomic model (Zenodo 10.5281/zenodo.6458345),
  PDB 8GBS (A. flos-aquae GvpA/GvpC integrative model, Dutka et al. 2023, Structure).
"""

import bpy, os, sys, math, importlib
import numpy as np
from mathutils import Matrix, Vector, Quaternion

# =============================== SETTINGS ====================================
RESOLUTION = (1920, 1080)
FPS = 30
ENGINE = "EEVEE"            # "EEVEE" (fast, minutes-hours) or "CYCLES" (nicer light, much slower)
EEVEE_SAMPLES = 64
CYCLES_SAMPLES = 96
USE_GPU_FOR_CYCLES = True
OUTPUT_PATH = "//render/gas_vesicle"   # "//" = folder of the saved .blend file
OUTPUT_FORMAT = "MP4"       # "MP4" (one video file) or "PNG" (image sequence, safest for long renders)
PREVIEW = False             # True: 50 % resolution + fewer samples for a quick test render

BACKGROUND = (1.0, 1.0, 1.0)     # white
COLORS = {
    "shell":       "#8FA3B8",   # GvpA in the shell
    "half_a":      "#2F5F98",   # highlight: helical half 1
    "half_b":      "#3E9C8F",   # highlight: helical half 2
    "gvpc":        "#EE8A2A",   # GvpC
    "hydrophobic": "#E3A92B",   # hydropathy colouring of the single GvpA
    "neutral":     "#E6E1D8",
    "hydrophilic": "#3F7FC1",
    "helix":       "#C8553D",   # cartoon colours
    "sheet":       "#E9B949",
    "loop":        "#9AA5B1",
    "density":     "#6F86A6",
}
DENSITY_THRESHOLD = 0.8          # EMD-14238 contour level (map units; ~2 sigma)
SHELL_SURFACE_QUALITY = 2        # 1-3; lower = lighter geometry for the 1730 shell copies

# Timeline (frames at 30 fps). Shift numbers here to re-time the story.
T = dict(
    end=1800,
    zoom_start=210, dissolve=(300, 390), hero_in=392,
    cartoon_in=(425, 455), surface_out=(430, 460),
    surface_in=(565, 595), cartoon_out=(585, 610),
    hero_to_shell_colour=(760, 800),
    patch_start=790, patch_step=7,
    density_in=(955, 990), density_out=(1045, 1075),
    morph=(1070, 1100), swap=1100, regrow=(1100, 1270),
    highlight_a=(1270, 1330), highlight_b=(1310, 1370), highlight_fade=(1395, 1425),
    gvpc=(1420, 1630),
)
DUR = dict(shell_in=15, shell_out=10, patch_in=22, gvpc_in=32)
# ==============================================================================


# ----------------------------------------------------------------------------- helpers
def find_data_dir():
    cands = []
    try:
        cands.append(os.path.dirname(os.path.abspath(__file__)))
    except NameError:
        pass
    for t in bpy.data.texts:
        if t.filepath:
            cands.append(os.path.dirname(bpy.path.abspath(t.filepath)))
    if bpy.data.filepath:
        cands.append(os.path.dirname(bpy.data.filepath))
    cands.append(os.getcwd())
    for c in cands:
        d = os.path.join(c, "data")
        if os.path.exists(os.path.join(d, "instances.npz")):
            return d
    raise FileNotFoundError("Could not find the 'data' folder. Keep it next to build_gv_animation.py "
                            "and open the script from disk in Blender's text editor.")


def import_mn():
    for name in ("bl_ext.blender_org.molecularnodes", "bl_ext.user_default.molecularnodes", "molecularnodes"):
        try:
            return importlib.import_module(name)
        except ImportError:
            pass
    for name in list(sys.modules):
        if name.endswith(".molecularnodes"):
            return sys.modules[name]
    raise ImportError("Molecular Nodes is not installed/enabled: Edit > Preferences > Get Extensions > 'Molecular Nodes'.")


def hex2rgb(h):
    h = h.lstrip("#")
    srgb = [int(h[i:i + 2], 16) / 255 for i in (0, 2, 4)]
    return tuple(c / 12.92 if c <= 0.04045 else ((c + 0.055) / 1.055) ** 2.4 for c in srgb)  # -> linear


def rgba(h):
    return (*hex2rgb(h), 1.0)


def sock(node, name, output=False):
    socks = node.outputs if output else node.inputs
    for s in socks:
        if s.name == name and s.enabled:
            return s
    for s in socks:
        if s.name == name:
            return s
    raise KeyError(f"{node.bl_idname} has no socket {name}")


def clear_scene():
    for ob in list(bpy.data.objects):
        bpy.data.objects.remove(ob, do_unlink=True)
    for c in list(bpy.data.collections):
        bpy.data.collections.remove(c)
    for block in (bpy.data.meshes, bpy.data.materials, bpy.data.cameras, bpy.data.lights,
                  bpy.data.volumes, bpy.data.actions, bpy.data.worlds):
        for b in list(block):
            if b.users == 0 or block in (bpy.data.actions,):
                try:
                    block.remove(b)
                except Exception:
                    pass
    for ng in list(bpy.data.node_groups):
        if ng.name.startswith("GV_"):
            bpy.data.node_groups.remove(ng)
    try:
        bpy.data.orphans_purge(do_recursive=True)
    except Exception:
        pass
    try:                                    # forget molecules of a previous run
        bpy.context.scene.MNSession.prune()
    except Exception:
        pass


def new_collection(name):
    c = bpy.data.collections.new(name)
    bpy.context.scene.collection.children.link(c)
    return c


def move_to(ob, coll):
    for c in list(ob.users_collection):
        c.objects.unlink(ob)
    coll.objects.link(ob)


def set_attr(mesh, name, values, atype, domain="POINT"):
    if name in mesh.attributes:
        mesh.attributes.remove(mesh.attributes[name])
    a = mesh.attributes.new(name, atype, domain)
    key = {"FLOAT": "value", "INT": "value", "FLOAT_VECTOR": "vector", "QUATERNION": "value",
           "FLOAT_COLOR": "color"}[atype]
    a.data.foreach_set(key, np.ascontiguousarray(values, dtype=np.float32 if atype != "INT" else np.int32).ravel())


def mat_to_quat_arr(Ms):
    return np.array([list(Matrix(M[:3, :3].tolist()).to_quaternion()) for M in Ms])  # w, x, y, z


# ----------------------------------------------------------------------------- setup
DATA = find_data_dir()
mn = import_mn()
if not hasattr(bpy.types.Scene, "MNSession"):
    mn.register()

clear_scene()
scene = bpy.context.scene
scene.frame_start, scene.frame_end = 1, T["end"]
scene.render.fps = FPS
scene.frame_set(1)
S = 0.01                                   # Molecular Nodes world scale: 1 Angstrom = 0.01 Blender units
inst = np.load(os.path.join(DATA, "instances.npz"))

col_src = new_collection("GV_sources (hidden)")
col_main = new_collection("GV_scene")

# root: vesicle axis (model z) -> Blender X, hero subunit faces the camera (-Y)
root = bpy.data.objects.new("GV_root", None)
col_main.objects.link(root)
P = Matrix(((0, 0, 1), (1, 0, 0), (0, 1, 0)))
root.matrix_world = P.to_4x4()
Pn = np.array(P)


def world(v_A):
    """model-frame Angstrom coordinates -> Blender world coordinates"""
    return Vector((Pn @ (np.asarray(v_A) * S)).tolist())


def scaled(M):
    M = np.array(M, dtype=float).copy()
    M[:3, 3] *= S
    return Matrix(M.tolist())


# ----------------------------------------------------------------------------- materials
def principled(name, color=None, rough=0.55, alpha_default=1.0, blended=False):
    m = bpy.data.materials.new(name)
    if m.node_tree is None:
        m.use_nodes = True
    nt = m.node_tree
    for n in list(nt.nodes):
        nt.nodes.remove(n)
    out = nt.nodes.new("ShaderNodeOutputMaterial"); out.location = (500, 0)
    bsdf = nt.nodes.new("ShaderNodeBsdfPrincipled"); bsdf.location = (200, 0); bsdf.name = "BSDF"
    sock(bsdf, "Roughness").default_value = rough
    try:
        sock(bsdf, "Specular IOR Level").default_value = 0.35
    except KeyError:
        pass
    if color:
        sock(bsdf, "Base Color").default_value = rgba(color)
    sock(bsdf, "Alpha").default_value = alpha_default
    nt.links.new(sock(bsdf, "BSDF", True), sock(out, "Surface"))
    try:
        m.surface_render_method = "BLENDED" if blended else "DITHERED"
    except Exception:
        pass
    return m, nt, bsdf


def value_node(nt, name, val, loc):
    v = nt.nodes.new("ShaderNodeValue"); v.name = v.label = name; v.location = loc
    v.outputs[0].default_value = val
    return v


def key_node_value(mat, node_name, pairs):
    s = mat.node_tree.nodes[node_name].outputs[0]
    for f, v in pairs:
        s.default_value = v
        s.keyframe_insert("default_value", frame=f)


def smath(nt, op, a, b=None, loc=(0, 0)):
    n = nt.nodes.new("ShaderNodeMath"); n.operation = op; n.location = loc
    for i, x in enumerate((a, b)):
        if x is None:
            continue
        if isinstance(x, (int, float)):
            n.inputs[i].default_value = x
        else:
            nt.links.new(x, n.inputs[i])
    return n.outputs[0]


def mix_rgb(nt, fac, a, b, loc=(0, 0)):
    n = nt.nodes.new("ShaderNodeMix"); n.data_type = "RGBA"; n.location = loc
    for s_, x in ((sock(n, "Factor"), fac), (sock(n, "A"), a), (sock(n, "B"), b)):
        if isinstance(x, tuple):
            s_.default_value = x
        elif isinstance(x, (int, float)):
            s_.default_value = x
        else:
            nt.links.new(x, s_)
    return sock(n, "Result", True)


# shell material: colour driven by per-instance attributes (filament id / position)
hero = int(inst["hero"])
fil_len = inst["fil_len"]
hero_pos = float(inst["shell_pos"][hero])
m_shell, nt, bsdf = principled("GV_Shell", COLORS["shell"], rough=0.5)
a_fil = nt.nodes.new("ShaderNodeAttribute"); a_fil.attribute_type = "INSTANCER"; a_fil.attribute_name = "fil"; a_fil.location = (-1400, 200)
a_pos = nt.nodes.new("ShaderNodeAttribute"); a_pos.attribute_type = "INSTANCER"; a_pos.attribute_name = "fpos"; a_pos.location = (-1400, -100)
LA = value_node(nt, "L_A", 0.0, (-1400, -350)); LB = value_node(nt, "L_B", 0.0, (-1400, -450)); HF = value_node(nt, "HL_fade", 1.0, (-1400, -550))
fil, fpos = sock(a_fil, "Factor", True), sock(a_pos, "Factor", True)
isA = smath(nt, "LESS_THAN", smath(nt, "ABSOLUTE", smath(nt, "SUBTRACT", fil, 0.0)), 0.5)
isB = smath(nt, "LESS_THAN", smath(nt, "ABSOLUTE", smath(nt, "SUBTRACT", fil, 1.0)), 0.5)
inA = smath(nt, "LESS_THAN", smath(nt, "ABSOLUTE", smath(nt, "SUBTRACT", fpos, hero_pos)), LA.outputs[0])
inB = smath(nt, "GREATER_THAN", fpos, smath(nt, "SUBTRACT", float(fil_len[1] - 1) + 0.5, LB.outputs[0]))
hA = smath(nt, "MULTIPLY", smath(nt, "MULTIPLY", isA, inA), HF.outputs[0])
hB = smath(nt, "MULTIPLY", smath(nt, "MULTIPLY", isB, inB), HF.outputs[0])
c1 = mix_rgb(nt, hA, rgba(COLORS["shell"]), rgba(COLORS["half_a"]), (-300, 100))
c2 = mix_rgb(nt, hB, c1, rgba(COLORS["half_b"]), (-100, 100))
nt.links.new(c2, sock(bsdf, "Base Color"))
SA = value_node(nt, "alpha", 1.0, (-300, -300)); nt.links.new(SA.outputs[0], sock(bsdf, "Alpha"))
key_node_value(m_shell, "L_A", [(T["highlight_a"][0], 0.0), ((T["highlight_a"][0] + T["highlight_a"][1]) // 2, 60.0),
                                (T["highlight_a"][1], float(fil_len.max()) + 5)])   # slow first turn, then sweep
key_node_value(m_shell, "L_B", [(T["highlight_b"][0], 0.0), (T["highlight_b"][1], float(fil_len.max()) + 5)])
key_node_value(m_shell, "HL_fade", [(T["highlight_fade"][0], 1.0), (T["highlight_fade"][1], 0.0)])

m_gvpc, _, _ = principled("GV_GvpC", COLORS["gvpc"], rough=0.45)

# hero surface: hydropathy colours (atom attribute "Color") blended with the shell colour
m_hsurf, nt, bsdf = principled("GV_HeroSurface", None, rough=0.5)
a_col = nt.nodes.new("ShaderNodeAttribute"); a_col.attribute_name = "Color"; a_col.location = (-700, 100)
HM = value_node(nt, "hydro_mix", 0.0, (-700, -150)); HA = value_node(nt, "alpha", 1.0, (-700, -300))
nt.links.new(mix_rgb(nt, HM.outputs[0], rgba(COLORS["shell"]), sock(a_col, "Color", True)), sock(bsdf, "Base Color"))
nt.links.new(HA.outputs[0], sock(bsdf, "Alpha"))

m_hcart, nt, bsdf = principled("GV_HeroCartoon", None, rough=0.4, alpha_default=0.0)
a_col = nt.nodes.new("ShaderNodeAttribute"); a_col.attribute_name = "Color"; a_col.location = (-400, 100)
CA_ = value_node(nt, "alpha", 0.0, (-400, -200))
nt.links.new(sock(a_col, "Color", True), sock(bsdf, "Base Color"))
nt.links.new(CA_.outputs[0], sock(bsdf, "Alpha"))

m_dens, nt, bsdf = principled("GV_Density", COLORS["density"], rough=0.6, alpha_default=0.0, blended=True)
DA = value_node(nt, "alpha", 0.0, (-400, -200))
nt.links.new(DA.outputs[0], sock(bsdf, "Alpha"))
try:
    m_dens.use_backface_culling = False
except Exception:
    pass

# ----------------------------------------------------------------------------- molecules
def load_mol(fname, name, coll):
    mol = mn.Molecule.load(os.path.join(DATA, fname), name=name)
    move_to(mol.object, coll)
    return mol


def mesh_attr(obj, name, comps=1):
    a = obj.data.attributes[name]
    arr = np.zeros(len(a.data) * comps, dtype=np.float32)
    a.data.foreach_get("value" if comps == 1 else "vector", arr)
    return arr.reshape(-1, comps) if comps > 1 else arr


def set_color_attr(obj, cols):
    me = obj.data
    if "Color" in me.attributes:
        me.attributes.remove(me.attributes["Color"])
    a = me.attributes.new("Color", "FLOAT_COLOR", "POINT")
    a.data.foreach_set("color", np.ascontiguousarray(cols, dtype=np.float32).ravel())


# (1) GvpA used for the 1730 shell copies
src_shell = load_mol("gvpa_monomer.cif", "GvpA_shell_source", col_src)
src_shell.add_style(mn.StyleSurface(quality=SHELL_SURFACE_QUALITY), color=None, material=m_shell)

# (2) GvpC repeat used for the cage
src_gvpc = load_mol("gvpc_repeat.cif", "GvpC_source", col_src)
src_gvpc.add_style(mn.StyleSurface(quality=2), color=None, material=m_gvpc)

# (3) the single "hero" GvpA: surface coloured by hydropathy + cartoon coloured by secondary structure
hero_M = scaled(inst["hero_M"])
hs = load_mol("gvpa_monomer.cif", "GvpA_single_surface", col_main)
kd = mesh_attr(hs.object, "b_factor")
t = np.clip(kd / 4.5, -1, 1)[:, None]
c_pho, c_neu, c_phi = (np.array(hex2rgb(COLORS[k])) for k in ("hydrophobic", "neutral", "hydrophilic"))
cols = np.where(t > 0, c_neu + t * (c_pho - c_neu), c_neu + (-t) * (c_phi - c_neu))
set_color_attr(hs.object, np.hstack([cols, np.ones((len(cols), 1))]))
hs.add_style(mn.StyleSurface(quality=3), color=None, material=m_hsurf)

hc = load_mol("gvpa_monomer.cif", "GvpA_single_cartoon", col_main)
if "sec_struct" in hc.object.data.attributes:
    ss = mesh_attr(hc.object, "sec_struct").astype(int)
    pal = {1: hex2rgb(COLORS["helix"]), 2: hex2rgb(COLORS["sheet"])}
    cols = np.array([pal.get(s_, hex2rgb(COLORS["loop"])) for s_ in ss])
    cartoon = mn.StyleCartoon(quality=3, peptide_rounded=True, peptide_width=2.6, peptide_thickness=0.9, peptide_loop_radius=0.45)
else:
    rid = mesh_attr(hc.object, "res_id"); x = (rid - rid.min()) / max(1, np.ptp(rid))
    cols = np.stack([0.2 + 0.7 * x, 0.35 + 0.3 * np.sin(np.pi * x), 0.9 - 0.7 * x], 1)
    cartoon = mn.StyleCartoon(quality=3, peptide_dssp=True, peptide_rounded=True, peptide_width=2.6, peptide_thickness=0.9, peptide_loop_radius=0.45)
set_color_attr(hc.object, np.hstack([cols, np.ones((len(cols), 1))]))
hc.add_style(cartoon, color=None, material=m_hcart)

for ob in (hs.object, hc.object):
    ob.parent = root
    ob.matrix_basis = hero_M

# (4) EMD-14238 density, placed with the same transform that maps 7R1C onto the hero subunit
dens = mn.entities.density.io.load(os.path.join(DATA, "emd_14238.mrc"), style="density_surface", overwrite=True)
move_to(dens.object, col_main)
dens.object.name = "EMD-14238"
for n in dens.tree.nodes:
    if n.bl_idname == "GeometryNodeGroup" and "Density" in n.node_tree.name:
        for s_ in n.inputs:
            if s_.name == "Threshold":
                s_.default_value = DENSITY_THRESHOLD
            elif s_.name == "Hide Dust":
                s_.default_value = 5.0          # removes isolated noise specks
            elif s_.name == "Material":
                s_.default_value = m_dens
dens.object.parent = root
dens.object.matrix_basis = scaled(inst["map_M"])

# hide the sources (they are only used as instance templates)
for ob in (src_shell.object, src_gvpc.object):
    ob.hide_render = True
    ob.location = (0, 0, 0)
col_src.hide_render = False

# ----------------------------------------------------------------------------- instancer node group
def build_instancer_group():
    ng = bpy.data.node_groups.new("GV_Instancer", "GeometryNodeTree")
    itf = ng.interface
    itf.new_socket("Geometry", in_out="INPUT", socket_type="NodeSocketGeometry")
    itf.new_socket("Source", in_out="INPUT", socket_type="NodeSocketObject")
    s_m = itf.new_socket("Morph", in_out="INPUT", socket_type="NodeSocketFloat"); s_m.min_value, s_m.max_value = 0, 1
    s_i = itf.new_socket("Duration In", in_out="INPUT", socket_type="NodeSocketFloat"); s_i.default_value = 15
    s_o = itf.new_socket("Duration Out", in_out="INPUT", socket_type="NodeSocketFloat"); s_o.default_value = 10
    itf.new_socket("Geometry", in_out="OUTPUT", socket_type="NodeSocketGeometry")
    N = ng.nodes; L = ng.links
    gin = N.new("NodeGroupInput"); gin.location = (-1600, 0)
    gout = N.new("NodeGroupOutput"); gout.location = (900, 0)
    time = N.new("GeometryNodeInputSceneTime"); time.location = (-1600, 400)
    frame = sock(time, "Frame", True)

    def attr(name, dtype, y):
        n = N.new("GeometryNodeInputNamedAttribute"); n.data_type = dtype; n.location = (-1400, y)
        sock(n, "Name").default_value = name
        return sock(n, "Attribute", True)

    def gmath(op, a, b=None, y=0, x=-1000):
        n = N.new("ShaderNodeMath"); n.operation = op; n.location = (x, y)
        for i, v in enumerate((a, b)):
            if v is None:
                continue
            if isinstance(v, (int, float)):
                n.inputs[i].default_value = v
            else:
                L.new(v, n.inputs[i])
        return n.outputs[0]

    def ramp(start, dur, y):
        n = N.new("ShaderNodeMapRange"); n.data_type = "FLOAT"; n.interpolation_type = "SMOOTHSTEP"; n.clamp = True
        n.location = (-800, y)
        L.new(frame, sock(n, "Value")); L.new(start, sock(n, "From Min"))
        L.new(gmath("ADD", start, dur, y), sock(n, "From Max"))
        return sock(n, "Result", True)

    din, dout = sock(gin, "Duration In", True), sock(gin, "Duration Out", True)
    a0, v0, a1, v1 = (attr(k, "FLOAT", 300 - 90 * i) for i, k in enumerate(("a0", "v0", "a1", "v1")))
    g0, h0 = ramp(a0, din, 400), ramp(v0, dout, 300)
    g1, h1 = ramp(a1, din, 200), ramp(v1, dout, 100)
    w0 = gmath("MULTIPLY", g0, gmath("SUBTRACT", 1.0, h0, 400, -600), 400, -450)
    w1 = gmath("MULTIPLY", g1, gmath("SUBTRACT", 1.0, h1, 200, -600), 200, -450)
    s = gmath("MAXIMUM", w0, w1, 300, -300)
    visible = gmath("GREATER_THAN", s, 0.0005, 300, -150)

    # position = mix(t0, t1, morph) + fly * (1-g1)^2 * [window 0 is over]
    t0, t1, fly = attr("t0", "FLOAT_VECTOR", -200), attr("t1", "FLOAT_VECTOR", -300), attr("fly", "FLOAT_VECTOR", -400)
    q0, q1 = attr("q0", "QUATERNION", -500), attr("q1", "QUATERNION", -600)
    morph = sock(gin, "Morph", True)
    mixv = N.new("ShaderNodeMix"); mixv.data_type = "VECTOR"; mixv.location = (-800, -250)
    L.new(morph, sock(mixv, "Factor")); L.new(t0, sock(mixv, "A")); L.new(t1, sock(mixv, "B"))
    after0 = gmath("GREATER_THAN", frame, gmath("ADD", v0, dout, -350, -1000), -350, -800)
    remain = gmath("SUBTRACT", 1.0, g1, -450, -1000)
    fly_f = gmath("MULTIPLY", gmath("MULTIPLY", remain, remain, -450, -800), after0, -450, -600)
    sc = N.new("ShaderNodeVectorMath"); sc.operation = "SCALE"; sc.location = (-400, -400)
    L.new(fly, sc.inputs[0]); L.new(fly_f, sock(sc, "Scale"))
    add = N.new("ShaderNodeVectorMath"); add.operation = "ADD"; add.location = (-200, -300)
    L.new(sock(mixv, "Result", True), add.inputs[0]); L.new(sc.outputs[0], add.inputs[1])
    mixq = N.new("ShaderNodeMix"); mixq.data_type = "ROTATION"; mixq.location = (-400, -600)
    L.new(morph, sock(mixq, "Factor")); L.new(q0, sock(mixq, "A")); L.new(q1, sock(mixq, "B"))

    setpos = N.new("GeometryNodeSetPosition"); setpos.location = (0, 0)
    L.new(sock(gin, "Geometry", True), sock(setpos, "Geometry")); L.new(add.outputs[0], sock(setpos, "Position"))
    oinfo = N.new("GeometryNodeObjectInfo"); oinfo.location = (0, 300); oinfo.transform_space = "ORIGINAL"
    L.new(sock(gin, "Source", True), sock(oinfo, "Object")); sock(oinfo, "As Instance").default_value = True
    iop = N.new("GeometryNodeInstanceOnPoints"); iop.location = (400, 0)
    L.new(sock(setpos, "Geometry", True), sock(iop, "Points"))
    L.new(visible, sock(iop, "Selection"))
    L.new(sock(oinfo, "Geometry", True), sock(iop, "Instance"))
    L.new(sock(mixq, "Result", True), sock(iop, "Rotation"))
    L.new(s, sock(iop, "Scale"))
    L.new(sock(iop, "Instances", True), sock(gout, "Geometry"))
    return ng


GN = build_instancer_group()


def make_instancer(name, source, Ms, a0, v0, a1, v1, fly, extra=None, Ms1=None, dur_in=15, dur_out=10):
    n = len(Ms)
    me = bpy.data.meshes.new(name)
    tr = np.array([M[:3, 3] for M in Ms]) * S
    me.vertices.add(n)
    me.vertices.foreach_set("co", tr.astype(np.float32).ravel())
    ob = bpy.data.objects.new(name, me)
    col_main.objects.link(ob)
    ob.parent = root
    q = mat_to_quat_arr(Ms)
    Ms1 = Ms if Ms1 is None else Ms1
    set_attr(me, "t0", tr, "FLOAT_VECTOR"); set_attr(me, "t1", np.array([M[:3, 3] for M in Ms1]) * S, "FLOAT_VECTOR")
    set_attr(me, "q0", q, "QUATERNION"); set_attr(me, "q1", mat_to_quat_arr(Ms1), "QUATERNION")
    for k, v in (("a0", a0), ("v0", v0), ("a1", a1), ("v1", v1)):
        set_attr(me, k, np.broadcast_to(np.asarray(v, float), (n,)), "FLOAT")
    set_attr(me, "fly", fly * S, "FLOAT_VECTOR")
    for k, v in (extra or {}).items():
        set_attr(me, k, np.asarray(v, float), "FLOAT")
    mod = ob.modifiers.new("GV_Instancer", "NODES"); mod.node_group = GN
    ids = {it.name: it.identifier for it in GN.interface.items_tree if getattr(it, "in_out", "") == "INPUT"}
    mod[ids["Source"]] = source
    mod[ids["Duration In"]] = float(dur_in)
    mod[ids["Duration Out"]] = float(dur_out)
    mod[ids["Morph"]] = 0.0
    return ob, mod, ids


NEVER, ALWAYS = 1.0e7, -1.0e7

# ---- the full shell (1730 GvpA)
shell_M = inst["shell_M"]
cent = shell_M[:, :3, 3]
d_hero = np.linalg.norm(cent - cent[hero], axis=1)
rank = np.argsort(np.argsort(d_hero)) / (len(d_hero) - 1)          # 0 = hero ... 1 = farthest
v0 = T["dissolve"][0] + (1 - rank) * (T["dissolve"][1] - T["dissolve"][0])   # far ones vanish first
v0[hero] = T["hero_in"]
a1 = T["regrow"][0] + rank * (T["regrow"][1] - T["regrow"][0])
patch_match = inst["patch_match"]
a1[patch_match] = T["swap"] - DUR["shell_in"]                         # fully grown exactly at the swap
a1[hero] = T["swap"] - DUR["shell_in"]
shell_ob, shell_mod, ids = make_instancer(
    "GV_shell_GvpA", src_shell.object, shell_M, ALWAYS, v0, a1, NEVER, np.zeros((len(shell_M), 3)),
    extra={"fil": inst["shell_fil"], "fpos": inst["shell_pos"]}, dur_in=DUR["shell_in"], dur_out=DUR["shell_out"])

# ---- the local patch (exact helical lattice of 7R1C inside EMD-14238), assembled around the hero
pn = inst["patch_n"]; pM = inst["patch_M"]
keep = pn != 0
pn, pM, pmatch = pn[keep], pM[keep], patch_match[keep]
rib = np.round(pn / 93.0); along = pn - 93 * rib
order = np.abs(along) + 6 * np.abs(rib)
pa1 = T["patch_start"] + (order - order.min()) * T["patch_step"]
pfly = (pM[:, :3, 3] - inst["hero_M"][:3, 3]) * 0.9
radial = pM[:, :3, 3].copy(); radial[:, 2] = 0
pfly += radial / np.linalg.norm(radial, axis=1)[:, None] * 25.0
patch_ob, patch_mod, _ = make_instancer(
    "GV_patch_GvpA", src_shell.object, pM, NEVER, ALWAYS, pa1, T["swap"], pfly,
    extra={"fil": np.full(len(pM), -1.0), "fpos": np.zeros(len(pM))}, Ms1=shell_M[pmatch],
    dur_in=DUR["patch_in"], dur_out=1)
for f, v in ((T["morph"][0], 0.0), (T["morph"][1], 1.0)):
    patch_mod[ids["Morph"]] = v
    patch_mod.keyframe_insert(f'["{ids["Morph"]}"]', frame=f)

# same patch as cartoons (secondary-structure colours), shown inside the density map
pcart_ob, pcart_mod, _ = make_instancer(
    "GV_patch_cartoon", hc.object, pM, NEVER, ALWAYS, pa1, T["swap"], pfly,
    extra={"fil": np.full(len(pM), -1.0), "fpos": np.zeros(len(pM))}, Ms1=shell_M[pmatch],
    dur_in=DUR["patch_in"], dur_out=1)

# ---- GvpC cage
gM = inst["gvpc_M"]; gmol = inst["gvpc_mol"]
nmol = gmol.max() + 1
mol_z = np.array([gM[gmol == m][:, 2, 3].mean() for m in range(nmol)])
rng = np.random.default_rng(7)
mol_order = np.argsort(np.argsort(mol_z + rng.normal(0, 120, nmol))) / max(1, nmol - 1)
ga1 = T["gvpc"][0] + mol_order[gmol] * (T["gvpc"][1] - T["gvpc"][0] - DUR["gvpc_in"])
rad = gM[:, :3, 3].copy(); rad[:, 2] = 0
gfly = rad / np.linalg.norm(rad, axis=1)[:, None] * 90.0
gvpc_ob, gvpc_mod, _ = make_instancer("GV_GvpC_cage", src_gvpc.object, gM, NEVER, ALWAYS, ga1, NEVER, gfly,
                                      dur_in=DUR["gvpc_in"], dur_out=10)

# ----------------------------------------------------------------------------- per-object visibility keys
def key_visibility(ob, spans):
    """spans: list of (first_visible_frame, last_visible_frame)"""
    frames = sorted({1, *[a for a, b in spans], *[b + 1 for a, b in spans]})
    for f in frames:
        vis = any(a <= f <= b for a, b in spans)
        ob.hide_render = ob.hide_viewport = not vis
        ob.keyframe_insert("hide_render", frame=f)
        ob.keyframe_insert("hide_viewport", frame=f)


key_visibility(hs.object, [(T["hero_in"], T["swap"])])
key_visibility(hc.object, [(T["cartoon_in"][0] - 1, T["cartoon_out"][1] + 1),
                           (T["density_in"][0] - 1, T["density_out"][1] + 1)])
key_visibility(pcart_ob, [(T["density_in"][0] - 1, T["density_out"][1] + 1)])
key_visibility(dens.object, [(T["density_in"][0] - 1, T["density_out"][1] + 1)])

din_, dout_ = T["density_in"], T["density_out"]
key_node_value(m_hsurf, "alpha", [(T["surface_out"][0], 1.0), (T["surface_out"][1], 0.0),
                                  (T["surface_in"][0], 0.0), (T["surface_in"][1], 1.0),
                                  (din_[0], 1.0), (din_[1], 0.0), (dout_[0], 0.0), (dout_[1], 1.0)])
# surfaces of the assembled patch give way to cartoons while the cryo-EM density is shown
key_node_value(m_shell, "alpha", [(din_[0], 1.0), (din_[1], 0.0), (dout_[0], 0.0), (dout_[1], 1.0)])
key_node_value(m_hsurf, "hydro_mix", [(T["surface_out"][0], 0.0), (T["surface_in"][0], 1.0),
                                      (T["hero_to_shell_colour"][0], 1.0), (T["hero_to_shell_colour"][1], 0.0)])
key_node_value(m_hcart, "alpha", [(T["cartoon_in"][0], 0.0), (T["cartoon_in"][1], 1.0),
                                  (T["cartoon_out"][0], 1.0), (T["cartoon_out"][1], 0.0),
                                  (din_[0], 0.0), (din_[1], 1.0), (dout_[0], 1.0), (dout_[1], 0.0)])
key_node_value(m_dens, "alpha", [(T["density_in"][0], 0.0), (T["density_in"][1], 0.55),
                                 (T["density_out"][0], 0.55), (T["density_out"][1], 0.0)])

# ----------------------------------------------------------------------------- world & lights
world_ = bpy.data.worlds.new("GV_World"); scene.world = world_
if world_.node_tree is None:
    world_.use_nodes = True
wn = world_.node_tree
for n in list(wn.nodes):
    wn.nodes.remove(n)
wout = wn.nodes.new("ShaderNodeOutputWorld"); wout.location = (400, 0)
lp = wn.nodes.new("ShaderNodeLightPath"); lp.location = (-400, 200)
bg_cam = wn.nodes.new("ShaderNodeBackground"); bg_cam.location = (-200, 0)
bg_cam.inputs[0].default_value = (*BACKGROUND, 1); bg_cam.inputs[1].default_value = 1.0
bg_light = wn.nodes.new("ShaderNodeBackground"); bg_light.location = (-200, -150)
bg_light.inputs[0].default_value = (1, 1, 1, 1); bg_light.inputs[1].default_value = 0.55
mixs = wn.nodes.new("ShaderNodeMixShader"); mixs.location = (150, 0)
wn.links.new(sock(lp, "Is Camera Ray", True), mixs.inputs[0])
wn.links.new(bg_light.outputs[0], mixs.inputs[1]); wn.links.new(bg_cam.outputs[0], mixs.inputs[2])
wn.links.new(mixs.outputs[0], wout.inputs[0])


def add_light(name, ltype, energy, rot_deg, size=None, color=(1, 1, 1)):
    ld = bpy.data.lights.new(name, ltype); ld.energy = energy; ld.color = color
    if ltype == "SUN":
        ld.angle = math.radians(size or 8)
    ob = bpy.data.objects.new(name, ld); col_main.objects.link(ob)
    ob.rotation_euler = [math.radians(a) for a in rot_deg]
    return ob


add_light("Key", "SUN", 3.2, (50, 0, -35), size=10, color=(1.0, 0.97, 0.93))
add_light("Fill", "SUN", 0.9, (65, 0, 140), size=30, color=(0.9, 0.95, 1.0))
add_light("Rim", "SUN", 1.2, (-60, 0, 20), size=15)

# ----------------------------------------------------------------------------- camera
cam_d = bpy.data.cameras.new("GV_Camera"); cam_d.lens = 50; cam_d.clip_start = 0.02; cam_d.clip_end = 500
cam = bpy.data.objects.new("GV_Camera", cam_d); col_main.objects.link(cam); scene.camera = cam
cam.rotation_mode = "QUATERNION"

O = Vector((0, 0, 0))
H = world(inst["hero_M"][:3, 3])
HP = world(inst["patch_M"][:, :3, 3].mean(0))                                            # centre of the local patch
n_out = world(np.array([inst["hero_M"][0, 3], inst["hero_M"][1, 3], 0.0])).normalized()   # outward normal at hero
ax = Vector((1, 0, 0))                                                                       # vesicle axis
up = n_out.cross(ax).normalized()
if up.z < 0:
    up = -up


def direction(az, el):
    az, el = math.radians(az), math.radians(el)
    return (math.cos(el) * (math.cos(az) * n_out + math.sin(az) * ax) + math.sin(el) * up).normalized()


# key poses: frame, target, azimuth, elevation, distance (Blender units; 1 BU = 100 A)
KEYS = [
    (1,    O,  -32, 22, 21.0),
    (T["zoom_start"], O,   -8, 12, 17.5),
    (T["dissolve"][0], H, 0, 12, 4.2),
    (T["hero_in"], H,       4, 10, 1.5),
    (470,  H,  25, 12, 1.35),
    (560,  H,  80,  6, 1.6),
    (650,  H, 135,  8, 1.65),
    (720,  H, 155, 12, 1.75),
    (800,  H,  60, 24, 2.9),
    (900,  HP, 22, 14, 5.4),
    (1070, HP, 12, 12, 5.6),
    (1300, O, -18, 14, 20.0),
    (1410, O,  -6, 14, 19.0),
    (1560, O+Vector((1.2, 0, 0)), 18, 22, 9.5),
    (1650, O+Vector((0.8, 0, 0)), 26, 20, 11.0),
    (T["end"], O, 40, 12, 19.5),
]


def pchip(xk, yk, x):
    xk, yk = np.asarray(xk, float), np.asarray(yk, float)
    h = np.diff(xk); d = np.diff(yk) / h
    m = np.zeros_like(yk)
    for i in range(1, len(yk) - 1):
        if d[i - 1] * d[i] > 0:
            w1, w2 = 2 * h[i] + h[i - 1], h[i] + 2 * h[i - 1]
            m[i] = (w1 + w2) / (w1 / d[i - 1] + w2 / d[i])
    i = np.clip(np.searchsorted(xk, x) - 1, 0, len(h) - 1)
    t_ = (x - xk[i]) / h[i]
    h00, h10, h01, h11 = 2*t_**3 - 3*t_**2 + 1, t_**3 - 2*t_**2 + t_, -2*t_**3 + 3*t_**2, t_**3 - t_**2
    return h00 * yk[i] + h10 * h[i] * m[i] + h01 * yk[i + 1] + h11 * h[i] * m[i + 1]


kf = [k[0] for k in KEYS]
frames = np.arange(1, T["end"] + 1)
tgt = np.stack([pchip(kf, [k[1][j] for k in KEYS], frames) for j in range(3)], 1)
az_ = pchip(kf, [k[2] for k in KEYS], frames)
el_ = pchip(kf, [k[3] for k in KEYS], frames)
dist = np.exp(pchip(kf, [math.log(k[4]) for k in KEYS], frames))
cam_d.dof.use_dof = True
cam_d.dof.aperture_fstop = 5.6
for i, f in enumerate(frames):
    tv = Vector(tgt[i].tolist())
    loc = tv + direction(az_[i], el_[i]) * dist[i]
    cam.location = loc
    cam.rotation_quaternion = (tv - loc).to_track_quat("-Z", "Y")
    cam.keyframe_insert("location", frame=int(f))
    cam.keyframe_insert("rotation_quaternion", frame=int(f))
    cam_d.dof.focus_distance = float(dist[i])
    cam_d.dof.keyframe_insert("focus_distance", frame=int(f))

# ----------------------------------------------------------------------------- render settings
r = scene.render
r.resolution_x, r.resolution_y = RESOLUTION
r.resolution_percentage = 50 if PREVIEW else 100
r.film_transparent = False
try:
    scene.view_settings.view_transform = "Standard"
    scene.view_settings.look = "None"
except Exception:
    pass
if ENGINE.upper() == "CYCLES":
    r.engine = "CYCLES"
    scene.cycles.samples = 16 if PREVIEW else CYCLES_SAMPLES
    scene.cycles.use_denoising = True
    if USE_GPU_FOR_CYCLES:
        scene.cycles.device = "GPU"
else:
    r.engine = "BLENDER_EEVEE"
    ee = scene.eevee
    ee.taa_render_samples = 16 if PREVIEW else EEVEE_SAMPLES
    for attr_, val in (("use_raytracing", True), ("use_shadows", True), ("fast_gi_method", "GLOBAL_ILLUMINATION")):
        try:
            setattr(ee, attr_, val)
        except Exception:
            pass
r.filepath = OUTPUT_PATH + ("_" if OUTPUT_FORMAT.upper() == "PNG" else ".mp4")
ims = r.image_settings
if OUTPUT_FORMAT.upper() == "MP4":
    try:
        ims.media_type = "VIDEO"
    except Exception:
        pass
    ims.file_format = "FFMPEG"
    r.ffmpeg.format = "MPEG4"
    r.ffmpeg.codec = "H264"
    r.ffmpeg.constant_rate_factor = "HIGH"
    r.ffmpeg.ffmpeg_preset = "GOOD"
else:
    try:
        ims.media_type = "IMAGE"
    except Exception:
        pass
    ims.file_format = "PNG"
    ims.color_mode = "RGB"

scene.frame_set(1)
print("Gas vesicle scene built:", len(shell_M), "GvpA,", len(pM), "patch copies,", len(gM), "GvpC repeats,",
      T["end"], "frames at", FPS, "fps. Output ->", r.filepath)
