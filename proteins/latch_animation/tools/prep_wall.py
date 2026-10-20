"""Turn the vesicle wall of the latch viewer into a plain 3D mesh: no atoms, no residues, nothing to click.

Run after prep_viewer.py (it reads and rewrites static/latch-viewer/ in place; running it twice changes nothing):

    python tools/prep_wall.py [static/latch-viewer]       # needs numpy + scipy + gemmi + scikit-image

- the wall (1725 copies of GvpA) plus the five GvpA the binder docks onto (the subunits prep_viewer.py left out of the wall,
  so there is no hole) become one density map -> one closed surface, written as latch-wall.ply.gz (the viewer unzips it);
- chains A-E and G, and the `gvpa` and `wall-*` assemblies with their operators, are taken out of latch-assembly.cif;
- latch-viewer.json: `wall` (list of assemblies) is replaced by `wallMesh` (the file name).
The surface is coarse on purpose (STEP / SIGMA): the wall is a backdrop, not a structure.
"""
import gzip
import json
import os
import sys

import gemmi
import numpy as np
from scipy.ndimage import gaussian_filter
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import connected_components
from skimage.measure import marching_cubes

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = sys.argv[1] if len(sys.argv) > 1 else os.path.join(HERE, "..", "..", "..", "static", "latch-viewer")
GV = os.path.join(HERE, "..", "..", "gv_animation", "data")
STEP = 2.0       # A, grid spacing of the density map = size of the smallest feature of the mesh
SIGMA = 3.0      # A, blur of the atoms
LEVEL = 0.5      # of the density inside the protein
GONE = set("ABCDEG")

cif_path = os.path.join(OUT, "latch-assembly.cif")
with open(cif_path) as f:
    header = [l for l in f if l.startswith("#") and not l.strip() == "#"][:3]
doc = gemmi.cif.read(cif_path)
blk = doc.sole_block()

# ---- the wall: every copy of G, plus the five subunits under the pentamer
st = gemmi.make_structure_from_block(blk)
chains = {c.name: c for c in st[0]}
if "G" not in chains:
    sys.exit("latch-assembly.cif has no wall any more: nothing to do (run prep_viewer.py first to start over)")


def xyz(chain):
    return np.array([a.pos.tolist() for r in chain for a in r])


G = xyz(chains["G"])
cols = [f"matrix[{i}][{j}]" for i in (1, 2, 3) for j in (1, 2, 3)] + [f"vector[{i}]" for i in (1, 2, 3)]
oper = {int(r[0]): np.array([float(x) for x in list(r)[1:]]) for r in blk.find("_pdbx_struct_oper_list.", ["id"] + cols)}
wall_ids = sorted(i for i in oper if i > 1)          # operator 1 is the identity: every other one is a wall subunit


def mat(v):
    M = np.eye(4)
    M[:3, :3] = v[:9].reshape(3, 3)
    M[:3, 3] = v[9:]
    return M


wall_ops = [mat(oper[i]) for i in wall_ids]

# the five subunits prep_viewer.py cut out (same fit as there: pentamer onto five consecutive wall subunits)
d = np.load(os.path.join(GV, "instances.npz"))
cen, fil, pos, shell_M = d["shell_cent"], d["shell_fil"], d["shell_pos"], d["shell_M"]
ours = np.array([xyz(chains[c]).mean(0) for c in "ABCDE"])


def kabsch(P, Q):
    pc, qc = P.mean(0), Q.mean(0)
    U, _, Vt = np.linalg.svd((P - pc).T @ (Q - qc))
    R = Vt.T @ np.diag([1, 1, np.sign(np.linalg.det(Vt.T @ U.T))]) @ U.T
    return R, qc - R @ pc


gv_ax0 = np.linalg.svd(cen - cen.mean(0), full_matrices=False)[2][0]
axial = np.abs((cen - cen.mean(0)) @ gv_ax0)
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
rms, _, R, tr, used = sorted((c for c in cands if c[0] <= 1.7), key=lambda c: c[1])[0]
K = np.eye(4)
K[:3, :3], K[:3, 3] = R, tr
keep = [i for i in range(len(shell_M)) if i not in set(used.tolist())]
mine = np.array([(K @ shell_M[i])[:3].ravel() for i in keep])
have = np.array([m[:3].ravel() for m in wall_ops])
gap = np.abs(mine[:, None] - have[None]).max(2).min(1)           # the file keeps 4 decimals
assert len(mine) == len(have) and gap.max() < 1e-2, "wall does not match prep_viewer.py's (%.4f)" % gap.max()
wall_ops += [K @ shell_M[i] for i in used]

# ---- density -> surface
P = np.vstack([G @ M[:3, :3].T + M[:3, 3] for M in wall_ops])
lo = P.min(0) - 3 * SIGMA
n = np.ceil((P.max(0) + 3 * SIGMA - lo) / STEP).astype(int) + 1
h, _ = np.histogramdd(P, bins=n, range=[(lo[i] - STEP / 2, lo[i] + (n[i] - 0.5) * STEP) for i in range(3)])
dens = gaussian_filter(h, SIGMA / STEP)
level = LEVEL * np.percentile(dens[dens > 0.1 * dens.max()], 75)
verts, faces, _, _ = marching_cubes(dens, level, spacing=(STEP, STEP, STEP))
verts += lo
# loose specks of the blur (a few dozen vertices each) are not wall
a, b = np.r_[faces[:, 0], faces[:, 1]], np.r_[faces[:, 1], faces[:, 2]]
_, comp = connected_components(coo_matrix((np.ones(len(a)), (a, b)), shape=(len(verts),) * 2), directed=False)
big = np.bincount(comp)[comp] > 2000
faces = faces[big[faces].all(1)]
used_v = np.unique(faces)
remap = np.zeros(len(verts), int)
remap[used_v] = np.arange(len(used_v))
verts, faces = verts[used_v], remap[faces]
# vertex normals from the faces (area weighted): Mol* would work them out too, but zero-area triangles, which marching cubes
# leaves where the density sits right on the level, give it NaN normals: dark lines across the wall
fn = np.cross(verts[faces[:, 1]] - verts[faces[:, 0]], verts[faces[:, 2]] - verts[faces[:, 0]])
normals = np.zeros_like(verts)
for k in range(3):
    np.add.at(normals, faces[:, k], fn)
normals /= np.maximum(np.linalg.norm(normals, axis=1), 1e-9)[:, None]
with gzip.open(os.path.join(OUT, "latch-wall.ply.gz"), "wt", compresslevel=9) as f:
    f.write("ply\nformat ascii 1.0\ncomment GvpA wall of the latch viewer (prep_wall.py): a surface, no atoms\n"
            f"element vertex {len(verts)}\nproperty float x\nproperty float y\nproperty float z\n"
            "property float nx\nproperty float ny\nproperty float nz\n"
            f"element face {len(faces)}\nproperty list uchar int vertex_indices\nend_header\n")
    for v, nn in zip(verts, normals):
        f.write("%.1f %.1f %.1f %.2f %.2f %.2f\n" % (*v, *nn))
    for t in faces:
        f.write("3 %d %d %d\n" % tuple(t))

# ---- take the wall out of the structure file and the manifest
GONE_ASM = {"gvpa"} | {r[0] for r in blk.find("_pdbx_struct_assembly.", ["id"]) if r[0].startswith("wall-")}


def rebuild(cat, keep_row):
    """rewrite one mmCIF category with only the rows keep_row accepts"""
    t = blk.find_mmcif_category(cat)
    tags = [x.split(".", 1)[1] for x in t.tags]
    rows = [list(r) for r in t if keep_row(dict(zip(tags, r)))]
    t.erase()
    loop = blk.init_mmcif_loop(cat, tags)
    for r in rows:
        loop.add_row(r)


rebuild("_entity_poly.", lambda r: r["entity_id"] not in GONE)
rebuild("_struct_asym.", lambda r: r["id"] not in GONE)
rebuild("_entity.", lambda r: r["id"] not in GONE)
rebuild("_atom_site.", lambda r: r["label_asym_id"] not in GONE)
rebuild("_pdbx_struct_assembly.", lambda r: r["id"] not in GONE_ASM)
rebuild("_pdbx_struct_assembly_gen.", lambda r: r["assembly_id"] not in GONE_ASM)
rebuild("_pdbx_struct_oper_list.", lambda r: int(r["id"]) not in wall_ids)
tmp = cif_path + ".tmp"
doc.write_file(tmp, gemmi.cif.Style.Pdbx)
with open(tmp) as f, open(cif_path, "w") as g:
    g.write("".join(header) + f.read())
os.remove(tmp)

mpath = os.path.join(OUT, "latch-viewer.json")
with open(mpath) as f:
    man = json.load(f)
man.pop("wall")
man["wallMesh"] = "latch-wall.ply.gz"
with open(mpath, "w") as f:
    json.dump(man, f, indent=1)

print("wall mesh: %d vertices, %d faces, %.1f MB; latch-assembly.cif %.0f kB"
      % (len(verts), len(faces), os.path.getsize(os.path.join(OUT, "latch-wall.ply.gz")) / 1e6, os.path.getsize(cif_path) / 1e3))
