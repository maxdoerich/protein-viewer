"""Build the data for the interactive Mol* viewer on the Model page from ../data/.

Usage:  python prep_viewer.py [out_dir]
  out_dir defaults to static/gv-viewer/ in the wiki repository.
Writes:
  gv-assembly.cif  GvpA (chain A) and one full-length GvpC (chain C) once each, plus every
                   copy of them as an mmCIF assembly operator. Mol* draws each copy as a
                   GPU instance of the same geometry, so the 1730-copy shell stays light.
  emd-14238.mrc    the EMD-14238 map from data/, cropped to the wall and stored as int8.
  gv-viewer.json   isovalue + camera frames for each step of the viewer.
Requires: numpy, gemmi, mrcfile

GvpC is the AlphaFold Server prediction data/gvpc_af.cif (job "Ser_GvpC_single"), subject to the
AlphaFold Server Output Terms of Use: see NOTICE.md. The 8GBS repeat (gvpc_repeat.cif) is used only
to say where on the shell each molecule sits.

Frame: everything is placed in the frame of the EMD-14238 map (= PDB 7R1C), so the map
needs no transform. In this frame the vesicle's long axis is z.
"""
import json, os, sys
import numpy as np, gemmi, mrcfile

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "..", "data")
REPO = os.path.abspath(os.path.join(HERE, "..", "..", ".."))
OUT = sys.argv[1] if len(sys.argv) > 1 else os.path.join(REPO, "static", "gv-viewer")
os.makedirs(OUT, exist_ok=True)

d = np.load(os.path.join(DATA, "instances.npz"))
F = np.linalg.inv(d["map_M"])                     # shell frame -> map frame


def place(Ms):
    return np.array([F @ M for M in Ms])


hero = int(d["hero"])
shell = place(d["shell_M"])                       # centred GvpA -> its place in the shell
hero_op = F @ d["hero_M"]
patch = place(d["patch_M"])
patch_n = d["patch_n"]
gvpc = place(d["gvpc_M"])                         # centred GvpC repeat -> along a rib

# --- which helical half each shell subunit belongs to (tips: nearest assigned neighbour)
fil = d["shell_fil"].copy()
cent = shell[:, :3, 3]
lone = np.where(fil < 0)[0]
done = np.where(fil >= 0)[0]
for k in lone:
    fil[k] = fil[done[np.argmin(np.linalg.norm(cent[done] - cent[k], axis=1))]]

# --- GvpC: the full-length AlphaFold model, placed as one rigid body on each 5-repeat molecule of the
# 8GBS-based placement (gvpc_mol groups the 280 repeats into 56 molecules). Kabsch fit of the model's five
# 32-residue stretches onto the five repeats; N to C follows the rib order. Molecules on the cones are
# left out: a straight helix does not follow a surface whose radius changes along it.
af = gemmi.read_structure(os.path.join(DATA, "gvpc_af.cif"))
af_ca = np.array([r["CA"][0].pos.tolist() for r in af[0][0]])
rep_ca = np.array([r["CA"][0].pos.tolist() for r in gemmi.read_structure(os.path.join(DATA, "gvpc_repeat.cif"))[0][0]])
STRETCH = [(0, 32), (33, 65), (66, 98), (99, 131), (131, 163)]
CYL_Z = 240.0                                     # the cylinder is |z| <= 240 in the shell frame (radius 183)
axis_pt, axis_dir = F[:3, 3], F[:3, 2] / np.linalg.norm(F[:3, 2])      # the vesicle axis in this frame


def radial(P):                                    # distance from the vesicle axis
    v = P - axis_pt
    return np.linalg.norm(v - np.outer(v @ axis_dir, axis_dir), axis=1)


def kabsch(P, Q):
    pc, qc = P.mean(0), Q.mean(0)
    U, _, Vt = np.linalg.svd((P - pc).T @ (Q - qc))
    R = Vt.T @ np.diag([1, 1, np.sign(np.linalg.det(Vt.T @ U.T))]) @ U.T
    return R, qc - R @ pc


gvpc_ops, fit_rms, dropped = [], [], 0
for m in sorted(set(d["gvpc_mol"].tolist())):
    ks = np.where(d["gvpc_mol"] == m)[0]
    assert len(ks) == 5, "a GvpC molecule is not five repeats"
    if np.abs(d["gvpc_M"][ks][:, 2, 3]).max() > CYL_Z:
        dropped += 1
        continue
    tgt = [(gvpc[k][:3, :3] @ rep_ca.T).T + gvpc[k][:3, 3] for k in ks]
    P = np.vstack([af_ca[a_:b_][:32] for a_, b_ in STRETCH])
    Q = np.vstack([t_[:32] for t_ in tgt])
    R, t = kabsch(P, Q)
    fit_rms.append(float(np.sqrt(((P @ R.T + t - Q) ** 2).sum(1).mean())))
    # a straight helix on a curved rib would dip into the wall in the middle: lift it out until it clears
    rod = af_ca @ R.T + t
    lift = max(0.0, radial(np.vstack(tgt)).min() - radial(rod).min())
    out_ = rod.mean(0) - axis_pt
    out_ -= (out_ @ axis_dir) * axis_dir
    t = t + lift * out_ / np.linalg.norm(out_)
    M = np.eye(4); M[:3, :3] = R; M[:3, 3] = t
    gvpc_ops.append(M)
gvpc_ops = np.array(gvpc_ops)
print("GvpC: %d molecules on the cylinder (%d on the cones left out), rigid fit rmsd mean %.1f max %.1f A"
      % (len(gvpc_ops), dropped, np.mean(fit_rms), np.max(fit_rms)))

# --- the patch: one rib = one turn of the GvpA helix: 360 / 3.874 deg = 92.9 subunits,
# 92.9 x 0.525 A = 48.8 A apart. rib 0 holds the hero (n = 0). In step 2 each copy drops
# in on its own (gv-viewer.js), in this order: the hero's rib from the hero outwards,
# then the rib before it, the one after, and the part of the next.
PER_TURN = 360 / 3.874
pc = patch[:, :3, 3]
rib = np.round(patch_n / PER_TURN).astype(int)
drop_order = sorted((i for i in range(len(patch_n)) if patch_n[i] != 0),
                    key=lambda i: (abs(rib[i]), rib[i], abs(patch_n[i] - rib[i] * PER_TURN)))

# --- the shell, for the growth animation (step 3). Each half is one helical filament
# of GvpA that grows at the central seam: every unit added there pushes the filament one
# step along its helix (it turns and moves out). Both filaments run from their tip
# (position 0) to the seam; the 13 subunits beyond each tip are there from the start.
# The filament is cut into blocks of BLOCK units: finished blocks move as one, the block
# being filled is drawn unit by unit (gv-viewer.js).
BLOCK = 32
spos = d["shell_pos"]
filament = {}                                         # half -> subunit indices, tip to seam
tips = {}                                             # half -> subunits outside the filament
for half in (0, 1):
    ks = np.where(fil == half)[0]
    on = ks[spos[ks] >= 0]
    filament[half] = on[np.argsort(spos[on])]
    tips[half] = ks[spos[ks] < 0]
axis_xy = (cent[:, :2].min(0) + cent[:, :2].max(0)) / 2   # the long axis runs along z through here

# --- one model: GvpA as chain A, the GvpC repeat as chain C -------------------------
a = gemmi.read_structure(os.path.join(DATA, "gvpa_monomer.cif"))
c = af
a.rename_chain(a[0][0].name, "A")
c.rename_chain(c[0][0].name, "C")
a[0].add_chain(c[0]["C"])
for r in a[0]["C"]:                               # B = pLDDT here, and GvpC carries no hydropathy: keep B neutral
    for at in r:
        at.b_iso = 0.0
for ch_, sub_ in (("A", "A"), ("C", "B")):          # two chains, two label_asym_ids: an assembly op names one,
    for r in a[0][ch_]:                              # and a shared label would put GvpC into every GvpA copy
        r.subchain = sub_
        r.entity_id = ""
a.setup_entities()
a.assign_label_seq_id()
# the model is one continuous helix (every i to i+3 CA distance < 6 A); say so, or Mol* draws it as coil
gc = a[0]["C"]
h = gemmi.Helix()
h.start = gemmi.AtomAddress("C", gc[0].seqid, gc[0].name, "CA")
h.end = gemmi.AtomAddress("C", gc[len(gc) - 1].seqid, gc[len(gc) - 1].name, "CA")
h.pdb_helix_class = gemmi.Helix.HelixClass.RAlpha
h.length = len(gc)
a.helices.append(h)
names = {"A": "GvpA (B. megaterium, PDB 7R1C)", "C": "GvpC (AlphaFold Server model Ser_GvpC_single, 163 residues)"}
sub_of = {ch.name: ch[0].subchain for ch in a[0]}
desc = {e.name: names[ch] for e in a.entities for ch, sub in sub_of.items() if sub in e.subchains}

doc = a.make_mmcif_document()
blk = doc.sole_block()
ent = [(i, t) for i, t in blk.find("_entity.", ["id", "type"])]
blk.find_mmcif_category("_entity.").erase()
loop = blk.init_mmcif_loop("_entity.", ["id", "type", "pdbx_description"])
for i, t in ent:
    loop.add_row([i, t, gemmi.cif.quote(desc.get(i, "?"))])
ids = blk.find_values("_struct_conf.id")
for i in range(len(ids)):
    ids[i] = f"HELX_P{i+1}"
for cat in ("_pdbx_struct_assembly.", "_pdbx_struct_assembly_gen.", "_pdbx_struct_oper_list.",
            "_cell.", "_symmetry."):
    blk.find_mmcif_category(cat).erase()
# label_asym_id of each chain (gemmi assigns them in setup_entities)
lab = {}
for row in blk.find("_atom_site.", ["auth_asym_id", "label_asym_id"]):
    lab.setdefault(row[0], row[1])

# every placement is written once; assemblies pick theirs by operator id
ops, gens, assemblies = [], [], []


def add_ops(mats):
    first = len(ops) + 1
    ops.extend(mats)
    return list(range(first, len(ops) + 1))


def expr(ids):
    """[1, 2, 3, 7, 9, 10] -> '1-3,7,9-10'"""
    out, i = [], 0
    while i < len(ids):
        j = i
        while j + 1 < len(ids) and ids[j + 1] == ids[j] + 1:
            j += 1
        out.append(str(ids[i]) if i == j else f"{ids[i]}-{ids[j]}")
        i = j + 1
    return ",".join(out)


def assembly(aid, details, op_ids, chain):
    gens.append((aid, "(" + expr(op_ids) + ")", lab[chain]))
    assemblies.append((aid, details, len(op_ids)))


hero_id = add_ops([hero_op])
shell_ids = add_ops(list(shell))
patch_ids = add_ops(list(patch))
gvpc_ids = add_ops(list(gvpc_ops))

assembly("hero", "the GvpA monomer in its place in the wall", hero_id, "A")
patch_list = []
for order, i in enumerate(drop_order):
    pid = f"p-{order:02d}"
    assembly(pid, f"patch copy n = {patch_n[i]} (rib {rib[i]}), on the EMD-14238 lattice", [patch_ids[i]], "A")
    patch_list.append({"id": pid, "rib": int(rib[i])})
assembly("mono", "one GvpA, centred (placed by the viewer)", add_ops([np.eye(4)]), "A")
growth = {}
for half in (0, 1):
    h_ = str(half + 1)
    assembly(f"h{h_}-tip", f"helical half {h_}: the tip, beyond the filament", [shell_ids[x] for x in tips[half]], "A")
    blocks = []
    for b0 in range(0, len(filament[half]), BLOCK):
        ks = [x for x in filament[half][b0:b0 + BLOCK] if x != hero]
        bid = f"h{h_}-b{b0 // BLOCK:02d}"
        assembly(bid, f"helical half {h_}: filament units {b0}-{b0 + BLOCK - 1} from the tip",
                 [shell_ids[x] for x in ks], "A")
        blocks.append(bid)
    growth[h_] = {"tip": f"h{h_}-tip", "blocks": blocks, "n": int(len(filament[half]))}
assembly("gvpc", "GvpC on the cylinder of the shell (illustrative coverage)", gvpc_ids, "C")

loop = blk.init_mmcif_loop("_pdbx_struct_assembly.", ["id", "details", "oligomeric_details", "oligomeric_count"])
for aid, det, n in assemblies:
    loop.add_row([aid, gemmi.cif.quote(det), "?", str(n)])
loop = blk.init_mmcif_loop("_pdbx_struct_assembly_gen.", ["assembly_id", "oper_expression", "asym_id_list"])
for aid, expr, asym_id in gens:
    loop.add_row([aid, gemmi.cif.quote(expr), asym_id])
cols = ["id", "type"] + [f"matrix[{i}][{j}]" for i in (1, 2, 3) for j in (1, 2, 3)] + [f"vector[{i}]" for i in (1, 2, 3)]
loop = blk.init_mmcif_loop("_pdbx_struct_oper_list.", cols)
for i, M in enumerate(ops, 1):
    loop.add_row([str(i), gemmi.cif.quote("general operation")]
                 + [f"{M[r, s]:.6f}" for r in range(3) for s in range(3)]
                 + [f"{M[r, 3]:.4f}" for r in range(3)])
doc.write_file(os.path.join(OUT, "gv-assembly.cif"), gemmi.cif.Style.Pdbx)
# chain C is derived from AlphaFold Server output: the terms ask for a conspicuous notice with every copy
with open(os.path.join(OUT, "gv-assembly.cif")) as f:
    body = f.read()
with open(os.path.join(OUT, "gv-assembly.cif"), "w") as f:
    f.write("# Chain C (GvpC) is derived from an AlphaFold Server prediction. This information is subject to the\n"
            "# AlphaFold Server Output Terms of Use found at alphafoldserver.com/output-terms (non-commercial use).\n" + body)

# --- the map: crop to the wall, quantise to int8 ----------------------------------
m = mrcfile.open(os.path.join(DATA, "emd_14238.mrc"))
vs = float(m.voxel_size.x)
org = np.array([m.header.origin.x, m.header.origin.y, m.header.origin.z], float)
den = m.data.astype(np.float32)                   # z, y, x
pos = np.argwhere(den > 0.3)
lo = np.maximum(pos.min(0) - 2, 0)
hi = np.minimum(pos.max(0) + 3, den.shape)
crop = den[lo[0]:hi[0], lo[1]:hi[1], lo[2]:hi[2]]
scale = 127.0 / float(crop.max())
q = np.clip(np.round(crop * scale), -128, 127).astype(np.int8)
with mrcfile.new(os.path.join(OUT, "emd-14238.mrc"), overwrite=True) as o:
    o.set_data(q)
    o.voxel_size = vs
    o.header.origin.x, o.header.origin.y, o.header.origin.z = org + lo[::-1] * vs
    o.header.nxstart = o.header.nystart = o.header.nzstart = 0
    o.update_header_stats()

# --- camera frames -----------------------------------------------------------------
z = np.array([0.0, 0.0, 1.0])
hc = hero_op[:3, 3]
out = np.array([hc[0], hc[1], 0.0]); out /= np.linalg.norm(out)   # radial: away from the axis
up = np.cross(out, z); up /= np.linalg.norm(up)                   # tangent: keeps the long axis horizontal


def sphere(P, pad):
    ctr = (P.min(0) + P.max(0)) / 2
    return ctr, float(np.linalg.norm(P - ctr, axis=1).max() + pad)


mono_xyz = np.array([at.pos.tolist() for r in a[0]["A"] for at in r])
mono_xyz = mono_xyz @ hero_op[:3, :3].T + hero_op[:3, 3]
pc_all = pc
allc = np.vstack([cent, gvpc_ops[:, :3, 3]])
lo_all, hi_all = allc.min(0) - 25, allc.max(0) + 25
frames = {
    # half-extents across the view (along `up`) and along the long axis (z): JS fits both
    "monomer": dict(target=sphere(mono_xyz, 2)[0].tolist(), radius=sphere(mono_xyz, 2)[1]),
    "rib": dict(target=sphere(pc_all, 0)[0].tolist(), radius=sphere(pc_all, 22)[1]),
    "vesicle": dict(target=((lo_all + hi_all) / 2).tolist(),
                    halfHeight=float((hi_all - lo_all)[:2].max() / 2), halfWidth=float((hi_all - lo_all)[2] / 2)),
}
# the 4-GvpA preview: four neighbours of rib 0 (n = -1, 0 = the hero, 1, 2), seen along the
# vesicle axis, so the wall is edge-on with the outside up (`out`) and the gas inside down
h4_i = [int(np.where(patch_n == n)[0][0]) for n in (-1, 1, 2)]
h4_ids = ["hero"] + [patch_list[drop_order.index(i)]["id"] for i in h4_i]
h4_ctr = np.mean([hc] + [patch[i][:3, 3] for i in h4_i], axis=0)
hydro4 = dict(ids=h4_ids, target=h4_ctr.tolist(), radius=42.0, dir=(-z).tolist(), up=out.tolist())

manifest = {
    "frame": "EMD-14238 / PDB 7R1C; vesicle axis = z",
    "view": {"dir": (-out).tolist(), "up": up.tolist()},
    "frames": frames,
    # 0.8 map units, the contour the animation uses (~2 sigma), in the int8 file's units
    "isovalue": round(0.8 * scale, 2),
    "heroHalf": int(fil[hero]) + 1,                # which half's colour the hero takes in the last step
    # step 2: the patch copies, in the order they drop in
    "patch": patch_list,
    # step 3: per half, the tip and the filament blocks (tip to seam); every filament
    # unit's placement is in gv-shell.json
    "growth": {"block": BLOCK, "axis": [round(float(v), 2) for v in axis_xy], "halves": growth,
               "hero": {"half": str(int(fil[hero]) + 1), "pos": int(np.where(filament[int(fil[hero])] == hero)[0][0])}},
    "hydro4": hydro4,
    "counts": {a_: n for a_, _, n in assemblies},
}
with open(os.path.join(OUT, "gv-viewer.json"), "w") as f:
    json.dump(manifest, f, indent=1)

# every filament unit's placement, tip to seam: 3x3 rotation (row-major) + translation
shell_json = {str(h + 1): [[round(float(v), 4) for v in np.r_[shell[k][:3, :3].ravel(), shell[k][:3, 3]]]
                           for k in filament[h]] for h in (0, 1)}
with open(os.path.join(OUT, "gv-shell.json"), "w") as f:
    json.dump(shell_json, f, separators=(",", ":"))

for fn in ("gv-assembly.cif", "emd-14238.mrc", "gv-viewer.json", "gv-shell.json"):
    print(f"{fn:18s} {os.path.getsize(os.path.join(OUT, fn)) / 1e3:8.1f} kB")
print("assemblies:", manifest["counts"], " isovalue:", manifest["isovalue"])
