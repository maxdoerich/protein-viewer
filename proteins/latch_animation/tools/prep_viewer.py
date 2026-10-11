"""Build the data for the interactive Mol* viewer of the latch (biosensor) mechanism.

Usage:
  python prep_viewer.py --trim <film_dir>   once: copy the inputs from the Drylab film folder
                                            (latch_anim/rank016/film/ in the Second Brain Vault) into ../data/
  python prep_viewer.py [out_dir]           build the viewer data from ../data/ and the gas vesicle data
                                            in ../../gv_animation/data/; out_dir defaults to
                                            static/latch-viewer/ in the wiki repository
Writes:
  latch-assembly.cif  every part of the construct once (binder, latch, target-binding domain, target, GvpA
                      pentamer) plus one GvpA monomer with an mmCIF assembly operator per copy in the wall.
                      Each part is its own assembly, so latch-viewer.js can move and fade it on its own.
  latch-viewer.json   hinge, fly-in vector, dock offset and camera frames.
Requires: numpy, gemmi

Frame: the film frame of the Drylab pipeline (latch_animation.py), i.e. the frame of the rank016 complex.
Everything in the film is a rigid body except the RFdiffusion trajectory, which the viewer does not use:
  - binder (112 aa, chain F) and the GvpA pentamer (chains A-E) stand still;
  - latch (120 aa, last frame of the trajectory = the closed design, chain L), target-binding domain
    (PD-L1 binder, 177 aa, chain P) and target (PD-L1, 211 aa, chain T) form ONE block that swings by
    `block_hinge` (axis, angle, centre) around the binder's C terminus; the target flies in along `fly_pdl1`;
  - the wall is the 1730-copy shell of the Huber 2023 model (Zenodo 10.5281/zenodo.6458345), moved into
    the film frame by fitting the pentamer onto five consecutive wall subunits (as latch_blender.py does).
The latch has no sequence (RFdiffusion backbone) and the target-binding domain is attached by a 3.8 A CA-CA
bridge: the geometry shows the design intent, it is not a prediction. See README.md.
"""
import json, os, sys
import numpy as np, gemmi

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "..", "data")
GV = os.path.join(HERE, "..", "..", "gv_animation", "data")
REPO = os.path.abspath(os.path.join(HERE, "..", "..", ".."))

FLY = 90.0             # A the target starts from, along place_domain.py's fly-in line (the Blender film flies 150 A)
AWAY = 65.0            # A the whole construct stands off the wall until the last step (as in the film)
GAP = 25.0             # A of room around the model for each camera frame
TRIM = 26              # PD-L1 binder residues 1-26 are left out: an N-terminal helix that runs 30 A away from the latch,
                       # so the folded domain would hang far from it. A drawn linker takes their place (see `linker`).
LINKER_STEP = 3.8      # A, CA-CA


# ------------------------------------------------------------------ trim: inputs -> ../data/
def pdb_model(path, idx):
    """the ATOM lines of the idx-th MODEL of a PDB file (the whole file if it has no MODEL records)"""
    out, n = [], -1
    for l in open(path):
        if l.startswith("MODEL"):
            n += 1
        elif l.startswith(("ATOM", "HETATM", "TER")) and n in (-1, idx):
            out.append(l)
    return "".join(out) + "END\n"


if len(sys.argv) > 2 and sys.argv[1] == "--trim":
    film = sys.argv[2]
    nt = json.load(open(os.path.join(film, "film_geometry.json")))["n_traj_frames"]
    os.makedirs(DATA, exist_ok=True)
    for src, dst, idx in (("film_binder", "binder", 0), ("film_gvpa", "gvpa", 0), ("film_latch", "latch", nt - 1),
                          ("film_pdl1binder", "tbd", 0), ("film_pdl1", "target", 0)):
        with open(os.path.join(DATA, dst + ".pdb"), "w") as f:
            f.write(pdb_model(os.path.join(film, src + ".pdb"), idx))
    geo = json.load(open(os.path.join(film, "film_geometry.json")))
    with open(os.path.join(DATA, "geometry.json"), "w") as f:
        json.dump(geo, f, indent=1)
    print("trimmed into", os.path.abspath(DATA))
    sys.exit(0)

OUT = sys.argv[1] if len(sys.argv) > 1 else os.path.join(REPO, "static", "latch-viewer")
os.makedirs(OUT, exist_ok=True)

geo = json.load(open(os.path.join(DATA, "geometry.json")))
bh = geo["block_hinge"]
hinge_c = np.array(bh["center"])
hinge_ax = np.array(bh["axis"]) / np.linalg.norm(bh["axis"])
hinge_deg = float(bh["deg"])
fly = np.array(geo["fly_pdl1"]) / np.linalg.norm(geo["fly_pdl1"])
iface = set(geo["interface_res"])


def kabsch(P, Q):
    pc, qc = P.mean(0), Q.mean(0)
    U, _, Vt = np.linalg.svd((P - pc).T @ (Q - qc))
    R = Vt.T @ np.diag([1, 1, np.sign(np.linalg.det(Vt.T @ U.T))]) @ U.T
    return R, qc - R @ pc


def rot(axis, deg):
    k = np.radians(deg)
    K = np.array([[0, -axis[2], axis[1]], [axis[2], 0, -axis[0]], [-axis[1], axis[0], 0]])
    return np.eye(3) + np.sin(k) * K + (1 - np.cos(k)) * K @ K


R_open = rot(hinge_ax, hinge_deg)
swing = lambda P: (P - hinge_c) @ R_open.T + hinge_c            # closed -> open, about the binder's C terminus


def xyz(chain):
    return np.array([a.pos.tolist() for r in chain for a in r])


# ------------------------------------------------------------------ one model, one chain per part
CHAINS = {  # chain, file, what it is
    "binder": ("F", "binder.pdb", "GvpA binder rank016 (112 aa, RFdiffusion + ProteinMPNN design)"),
    "latch": ("L", "latch.pdb", "Latch (120 aa, RFdiffusion backbone, poly-Gly: no sequence yet)"),
    "tbd": ("P", "tbd.pdb", "Target-binding domain: PD-L1 binder (residues 27-177)"),
    "target": ("T", "target.pdb", "Target: PD-L1 ectodomain (residues 18-228)"),
}
PLACE = os.path.join(DATA, "placement.json")      # written by place_domain.py: where the domain + PD-L1 sit relative to the binder
place = json.load(open(PLACE)) if os.path.exists(PLACE) else None


def move(chain, R, t):
    for r in chain:
        for a in r:
            a.pos = gemmi.Position(*(R @ np.array(a.pos.tolist()) + t))


st = None
for key, (ch, fn, _) in CHAINS.items():
    s = gemmi.read_structure(os.path.join(DATA, fn))
    s.rename_chain(s[0][0].name, ch)
    if key == "tbd":
        del s[0][ch][:TRIM]
    if place and key in ("tbd", "target"):
        move(s[0][ch], np.array(place["R"]), np.array(place["t"]))
    if st is None:
        st = s
    else:
        st[0].add_chain(s[0][ch])
penta = gemmi.read_structure(os.path.join(DATA, "gvpa.pdb"))
for ch in penta[0]:
    st[0].add_chain(ch)
mono = gemmi.read_structure(os.path.join(GV, "gvpa_monomer.cif"))
mono.rename_chain(mono[0][0].name, "G")
st[0].add_chain(mono[0]["G"])

# Bridge and linker. Binder, latch and domain are chains of their own (so each can move and fade alone), and a cartoon stops at the
# end of its chain: the 3.8 A CA-CA joins would show as gaps. Each bridge is a two-residue chain, the last residue before
# the join and the first after it, drawn with the part that arrives later (J with the latch, K with the domain).
def bridge(name, before, after):
    ch = gemmi.Chain(name)
    for res in (before[len(before) - 1], after[0]):
        ch.add_residue(res)
    for i, r in enumerate(ch, 1):
        r.seqid = gemmi.SeqId(i, " ")
    st[0].add_chain(ch)


def linker(name, before, after):
    """latch C terminus -> first residue of the domain: a straight, schematic Gly linker (CA every 3.8 A), drawn as a chain of
    its own that starts and ends on the two residues it joins. Straight because nothing is known about its fold: the real
    construct has a flexible linker here (the AlphaFold pilot used GGSGGS), which this stands for."""
    a, b = before[len(before) - 1], after[0]
    p0, p1 = np.array(a["CA"][0].pos.tolist()), np.array(b["CA"][0].pos.tolist())
    n = max(1, int(round(np.linalg.norm(p1 - p0) / LINKER_STEP)) - 1)       # residues in between
    u = (p1 - p0) / np.linalg.norm(p1 - p0)
    w = np.cross(u, [0, 0, 1.0]); w /= np.linalg.norm(w)
    ch = gemmi.Chain(name)
    ch.add_residue(a)
    for i in range(1, n + 1):
        q = p0 + (p1 - p0) * i / (n + 1)
        res = gemmi.Residue(); res.name = "GLY"
        for an, el, pos in (("N", "N", q - 1.2 * u), ("CA", "C", q), ("C", "C", q + 1.2 * u), ("O", "O", q + 1.2 * u + 1.2 * w)):
            at = gemmi.Atom(); at.name = an; at.element = gemmi.Element(el); at.pos = gemmi.Position(*pos)
            res.add_atom(at)
        ch.add_residue(res)
    ch.add_residue(b)
    for i, r in enumerate(ch, 1):
        r.seqid = gemmi.SeqId(i, " ")
    st[0].add_chain(ch)
    return n


bridge("J", st[0]["F"], st[0]["L"])                   # binder C terminus -> latch
n_link = linker("K", st[0]["L"], st[0]["P"])                   # latch C terminus -> target-binding domain
print("linker: %d residues between latch and domain" % n_link)

# Clash view: with the latch closed and PD-L1 on its domain, which residues of binder and PD-L1 overlap? (heavy atoms
# closer than CLASH_A; the binder file has no side chains, so this is a lower bound.)
CLASH_A = 2.8
from scipy.spatial import cKDTree
fa = [(r, a) for r in st[0]["F"] for a in r]
ta = [(r, a) for r in st[0]["T"] for a in r]
tf = cKDTree(np.array([a.pos.tolist() for _, a in fa]))
clash_f, clash_t, n_pairs = set(), set(), 0
for (r, a), hits in zip(ta, tf.query_ball_point(np.array([a.pos.tolist() for _, a in ta]), CLASH_A)):
    if hits:
        clash_t.add(r.seqid.num)
        clash_f.update(fa[h][0].seqid.num for h in hits)
        n_pairs += len(hits)
print("clash (closed, PD-L1 bound): %d PD-L1 residues x %d binder residues, %d atom pairs < %.1f A"
      % (len(clash_t), len(clash_f), n_pairs, CLASH_A))
for r in st[0]["F"]:                                  # B-factor: 1 = GvpA-binding site, 2 = in clash with PD-L1; the viewer colours by it
    for a in r:
        a.b_iso = 2.0 if r.seqid.num in clash_f else (1.0 if r.seqid.num in iface else 0.0)
for r in st[0]["T"]:                                  # PD-L1: 1 = in clash with the binder
    for a in r:
        a.b_iso = 1.0 if r.seqid.num in clash_t else 0.0
for name in ("L", "P", "A", "B", "C", "D", "E", "G", "J", "K"):
    for r in st[0][name]:
        for a in r:
            a.b_iso = 0.0
for ch in st[0]:
    for r in ch:
        r.subchain = ch.name
        r.entity_id = ""
st.setup_entities()
for e in st.entities:                                 # PDB files carry no SEQRES: without it label_seq_id stays empty
    ch = next(c for c in st[0] if c[0].subchain in e.subchains)     # and Mol* cannot chain the residues up
    e.full_sequence = [r.name for r in ch]
st.assign_label_seq_id()

names = {v[0]: v[2] for v in CHAINS.values()}
names.update({c: "GvpA, wall-binding site (pentamer)" for c in "ABCDE"})
names["G"] = "GvpA (B. megaterium, PDB 7R1C), one copy per wall subunit"
names["J"] = "Bridge: binder C terminus to latch"
names["K"] = "Linker: latch C terminus to target-binding domain (schematic, straight)"
sub_of = {ch.name: ch[0].subchain for ch in st[0]}
desc = {e.name: names[ch] for e in st.entities for ch, sub in sub_of.items() if sub in e.subchains}

doc = st.make_mmcif_document()
blk = doc.sole_block()
ent = [(i, t) for i, t in blk.find("_entity.", ["id", "type"])]
blk.find_mmcif_category("_entity.").erase()
loop = blk.init_mmcif_loop("_entity.", ["id", "type", "pdbx_description"])
for i, t in ent:
    loop.add_row([i, t, gemmi.cif.quote(desc.get(i, "?"))])
for cat in ("_pdbx_struct_assembly.", "_pdbx_struct_assembly_gen.", "_pdbx_struct_oper_list.", "_cell.", "_symmetry."):
    blk.find_mmcif_category(cat).erase()
lab = {}
for row in blk.find("_atom_site.", ["auth_asym_id", "label_asym_id"]):
    lab.setdefault(row[0], row[1])

# ------------------------------------------------------------------ the wall, in the film frame
d = np.load(os.path.join(GV, "instances.npz"))
cen, fil, pos, shell_M = d["shell_cent"], d["shell_fil"], d["shell_pos"], d["shell_M"]
pent = [xyz(st[0][c]) for c in "ABCDE"]
ours = np.array([p.mean(0) for p in pent])
cc = cen - cen.mean(0)
gv_ax0 = np.linalg.svd(cc, full_matrices=False)[2][0]
axial = np.abs(cc @ gv_ax0)
cands = []
for f in np.unique(fil):
    idx = np.where((fil == f) & (pos >= 0))[0]
    idx = idx[np.argsort(pos[idx])]
    for k in range(len(idx) - 4):
        run_ = idx[k:k + 5]
        for tgt in (cen[run_], cen[run_][::-1]):
            R, tr = kabsch(tgt, ours)
            rms = float(np.sqrt(((tgt @ R.T + tr - ours) ** 2).sum(1).mean()))
            cands.append((rms, axial[run_].mean(), R, tr, run_.copy()))
best = min(c[0] for c in cands)
# the film takes rmsd <= 1.6 x best (1.08 A) and lands 245 A from the equator, where the wall starts to curve into the
# cone; the viewer looks straight at the wall, so it accepts 1.7 A and gets a flat piece of the cylinder instead
good = sorted((c for c in cands if c[0] <= 1.7), key=lambda c: c[1])
rms, ax_d, R, tr, used = good[0]               # of the good fits, the one nearest the equator (the tips are cones)
print("wall: pentamer fitted onto 5 wall subunits, rmsd %.2f A (best %.2f), %.0f A from the equator" % (rms, best, ax_d))
K = np.eye(4); K[:3, :3] = R; K[:3, 3] = tr
keep = [i for i in range(len(shell_M)) if i not in set(used.tolist())]
wall_ops = [K @ shell_M[i] for i in keep]
wall_c = cen @ R.T + tr                       # centroids of all subunits, film frame
gv_pt = wall_c.mean(0)
gv_ax = np.linalg.svd(wall_c - gv_pt, full_matrices=False)[2][0]

# ------------------------------------------------------------------ assemblies
ops, gens, assemblies = [], [], []


def add_ops(mats):
    first = len(ops) + 1
    ops.extend(mats)
    return list(range(first, len(ops) + 1))


def expr(ids):
    out, i = [], 0
    while i < len(ids):
        j = i
        while j + 1 < len(ids) and ids[j + 1] == ids[j] + 1:
            j += 1
        out.append(str(ids[i]) if i == j else f"{ids[i]}-{ids[j]}")
        i = j + 1
    return ",".join(out)


def assembly(aid, details, op_ids, chains):
    gens.append((aid, "(" + expr(op_ids) + ")", ",".join(lab[c] for c in chains)))
    assemblies.append((aid, details, len(op_ids)))


one = add_ops([np.eye(4)])
assembly("binder", names["F"], one, "F")
assembly("latch", names["L"], one, "LJ")
assembly("tbd", names["P"], one, "PK")
assembly("target", names["T"], one, "T")
assembly("gvpa", "GvpA pentamer the binder docks onto", one, "ABCDE")
assembly("wall", "wall of the gas vesicle: %d GvpA copies" % len(wall_ops), add_ops(wall_ops), "G")

loop = blk.init_mmcif_loop("_pdbx_struct_assembly.", ["id", "details", "oligomeric_details", "oligomeric_count"])
for aid, det, n in assemblies:
    loop.add_row([aid, gemmi.cif.quote(det), "?", str(n)])
loop = blk.init_mmcif_loop("_pdbx_struct_assembly_gen.", ["assembly_id", "oper_expression", "asym_id_list"])
for aid, e, asym in gens:
    loop.add_row([aid, gemmi.cif.quote(e), asym])
cols = ["id", "type"] + [f"matrix[{i}][{j}]" for i in (1, 2, 3) for j in (1, 2, 3)] + [f"vector[{i}]" for i in (1, 2, 3)]
loop = blk.init_mmcif_loop("_pdbx_struct_oper_list.", cols)
for i, M in enumerate(ops, 1):
    loop.add_row([str(i), gemmi.cif.quote("general operation")]
                 + [f"{M[r, s]:.6f}" for r in range(3) for s in range(3)] + [f"{M[r, 3]:.4f}" for r in range(3)])
path = os.path.join(OUT, "latch-assembly.cif")
doc.write_file(path, gemmi.cif.Style.Pdbx)
with open(path) as f:
    body = f.read()
with open(path, "w") as f:
    # chains P and T most likely come from AlphaFold Server predictions (file names fold_pdl1_*): confirm, see README.md
    f.write("# Chains P and T (PD-L1 binder, PD-L1) are probably derived from AlphaFold Server predictions. If so, this\n"
            "# information is subject to the AlphaFold Server Output Terms of Use found at\n"
            "# alphafoldserver.com/output-terms (non-commercial use).\n" + body)

# ------------------------------------------------------------------ camera frames
P_binder, P_latch, P_tbd, P_target = (xyz(st[0][c]) for c in "FLPT")
P_pent = np.vstack(pent)
away_dir = P_binder.mean(0) - P_pent.mean(0)
away_dir /= np.linalg.norm(away_dir)
away = away_dir * AWAY
site = ours.mean(0)
out_ = site - gv_pt
out_ -= (out_ @ gv_ax) * gv_ax
out_ /= np.linalg.norm(out_)
print("away offset points %+.2f along the wall's outward normal (should be > 0)" % (away_dir @ out_))


def sphere(P):
    c = (P.min(0) + P.max(0)) / 2
    return c, float(np.linalg.norm(P - c, axis=1).max() + GAP)


def frame(P, dir_, up):
    c, r = sphere(P)
    return dict(target=c.round(2).tolist(), radius=round(r, 1), dir=dir_.round(4).tolist(), up=up.round(4).tolist())


# steps 1-4: look down the least-spread axis of everything that ever shows (closed and open), so the swing
# happens in the picture plane; the construct stands off by `away`
everything = np.vstack([P_binder, P_latch, swing(P_latch), P_tbd, swing(P_tbd), P_target, swing(P_target)])
Vt = np.linalg.svd(everything - everything.mean(0), full_matrices=False)[2]
view_dir, view_up = Vt[2], Vt[0]
sw_target = swing(P_target)
fly_dir = np.array(place["fly_closed"]) if place else fly
frames = {
    "binder": frame(P_binder + away, view_dir, view_up),
    "latch": frame(np.vstack([P_binder, P_latch]) + away, view_dir, view_up),
    "tbd": frame(np.vstack([P_binder, P_latch, P_tbd]) + away, view_dir, view_up),
    # step 4: the block swings open (no PD-L1 yet), then PD-L1 flies in to the open domain; the viewer flies in the closed
    # frame and swings the result, so the way in is swing(closed bound pose + fly * s)
    "mechanism": frame(np.vstack([P_binder, P_latch, swing(P_latch), P_tbd, swing(P_tbd), sw_target,
                                  swing(P_target + fly_dir * FLY * 0.5)]) + away, view_dir, view_up),
}
# the clash view: latch shut, PD-L1 on its domain. Look at the contact from the side: along the in-plane axis across which
# binder and PD-L1 spread least, with the line binder -> PD-L1 vertical, so the overlap sits in the middle of the picture
v = P_target.mean(0) - P_binder.mean(0)
v /= np.linalg.norm(v)
both = np.vstack([P_binder, P_target])
flat = both - np.outer(both @ v, v)
w = np.linalg.svd(flat - flat.mean(0), full_matrices=False)[2]
side = w[1] - (w[1] @ v) * v                                  # in-plane axis with the least spread
side /= np.linalg.norm(side)
frames["clash"] = frame(both + away, side, v)
# step 5: from outside the wall, looking radially in; the vesicle's long axis runs across the picture
docked = np.vstack([P_binder, swing(P_latch), swing(P_tbd), sw_target])
near = wall_c[np.linalg.norm(wall_c - site, axis=1) < 60]
frames["dock"] = frame(np.vstack([docked, docked + away, P_pent, near]), -out_, gv_ax)

manifest = {
    "frame": "film frame of latch_animation.py (rank016 complex); standing off by `away` until the last step",
    "hinge": {"center": hinge_c.round(3).tolist(), "axis": hinge_ax.round(4).tolist(), "deg": hinge_deg},
    "fly": {"dir": (np.array(place["fly_closed"]) if place else fly).round(4).tolist(), "dist": FLY},
    "clash": {"binder_residues": len(clash_f), "target_residues": len(clash_t), "atom_pairs": n_pairs, "max_A": CLASH_A},
    "away": away.round(2).tolist(),
    "frames": frames,
    "counts": {a: n for a, _, n in assemblies},
    "geometry": {**{k: bh[k] for k in ("iface_covered_closed", "iface_covered_open", "n_iface")},
                 "placement": place["metrics"] if place else None},
}
with open(os.path.join(OUT, "latch-viewer.json"), "w") as f:
    json.dump(manifest, f, indent=1)

for fn in ("latch-assembly.cif", "latch-viewer.json"):
    print(f"{fn:20s} {os.path.getsize(os.path.join(OUT, fn)) / 1e3:8.1f} kB")
print("assemblies:", manifest["counts"])
