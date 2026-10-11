"""Prepare all structural data for the gas-vesicle animation (not needed for rendering).

Usage:  python prep_data.py <raw_folder>
  <raw_folder> must contain: 7R1C.cif, 8GBS.cif, emd_14238.map.gz and
  fullGV.cif (= HuberEtAl_B_megaterium_GV_pseudoAtomicModel.cif from Zenodo 10.5281/zenodo.6458345)
Writes: ../data/  (gvpa_monomer.cif, gvpc_repeat.cif, emd_14238.mrc, instances.npz)
Requires: numpy, scipy, gemmi, mrcfile
"""
import os, sys
import numpy as np, gemmi, mrcfile
from scipy.spatial import cKDTree
RAW = (sys.argv[1] if len(sys.argv) > 1 else "raw").rstrip("/") + "/"
OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "data") + "/"
KD = dict(ILE=4.5,VAL=4.2,LEU=3.8,PHE=2.8,CYS=2.5,MET=1.9,ALA=1.8,GLY=-0.4,THR=-0.7,SER=-0.8,TRP=-0.9,
          TYR=-1.3,PRO=-1.6,HIS=-3.2,GLU=-3.5,GLN=-3.5,ASP=-3.5,ASN=-3.5,LYS=-3.9,ARG=-4.5)

def kabsch(P, Q):
    pc, qc = P.mean(0), Q.mean(0); H = (P-pc).T @ (Q-qc); U, S, Vt = np.linalg.svd(H)
    d = np.sign(np.linalg.det(Vt.T @ U.T)); R = Vt.T @ np.diag([1, 1, d]) @ U.T
    M = np.eye(4); M[:3, :3] = R; M[:3, 3] = qc - R @ pc; return M
def T(v): M = np.eye(4); M[:3, 3] = v; return M
def app(M, X): return X @ M[:3, :3].T + M[:3, 3]

# ---------- GvpA monomer (7R1C, chain N) -------------------------------------
s7 = gemmi.read_structure(RAW + "7R1C.cif"); s7.setup_entities()
ch7 = s7[0]["N"]
ref_CA = np.array([r["CA"][0].pos.tolist() for r in ch7])
c_mono = ref_CA.mean(0)
mono = gemmi.read_structure(RAW + "7R1C.cif"); mono.setup_entities()
for r in mono[0]["N"]:
    for a in r:
        a.pos = gemmi.Position(*(np.array(a.pos.tolist()) - c_mono)); a.b_iso = KD.get(r.name, 0.0)   # Kyte-Doolittle hydropathy as B-factor
doc = mono.make_mmcif_document()
ids = doc.sole_block().find_values("_struct_conf.id")
for i in range(len(ids)): ids[i] = f"HELX_P{i+1}"          # MN detects helices by the 'HELX' prefix
doc.write_file(OUT + "gvpa_monomer.cif")
# radial check: hydrophobic residues should point to the vesicle interior (smaller radius)
rad = {}
for r in ch7:
    sc = [a for a in r if a.name not in ("N", "CA", "C", "O")]
    if sc: rad[r.name + str(r.seqid.num)] = (np.hypot(*np.mean([a.pos.tolist() for a in sc], 0)[:2]), KD[r.name])
v = np.array(list(rad.values())); print("corr(radius, hydropathy) =", np.corrcoef(v[:, 0], v[:, 1])[0, 1].round(2), "(negative = hydrophobic inside)")

# ---------- full vesicle model (Zenodo) --------------------------------------
docf = gemmi.cif.read(RAW + "fullGV.cif"); CA = []; blk = []
for bi, b in enumerate(docf):
    st = gemmi.make_structure_from_block(b)
    for m in st: CA.append([r["CA"][0].pos.tolist() for r in m[0]]); blk.append(bi)
CA = np.array(CA); blk = np.array(blk); N = len(CA); cents = CA.mean(1)
Mfull = np.array([kabsch(ref_CA, CA[k]) for k in range(N)])      # 7R1C (uncentred) -> subunit k
shell_M = np.array([M @ T(c_mono) for M in Mfull])               # centred monomer -> subunit k

# filament tracing via the helical neighbour operator
th = np.deg2rad(-3.874)
def Hn(n):
    a = th*n; M = np.eye(4); M[:2, :2] = [[np.cos(a), -np.sin(a)], [np.sin(a), np.cos(a)]]; M[2, 3] = 0.525*n; return M
tree = cKDTree(cents); cm = np.append(c_mono, 1)
nx = np.full(N, -1); pv = np.full(N, -1)
for k in range(N):
    d, j = tree.query((Mfull[k] @ Hn(1) @ cm)[:3]);  nx[k] = j if d < 4 else -1
    d, j = tree.query((Mfull[k] @ Hn(-1) @ cm)[:3]); pv[k] = j if d < 4 else -1
for k in range(N):
    if nx[k] >= 0 and pv[nx[k]] != k: nx[k] = -1
for k in range(N):
    if pv[k] >= 0 and nx[pv[k]] != k: pv[k] = -1
fil = np.full(N, -1); fpos = np.full(N, -1); chains = []
for s in [k for k in range(N) if pv[k] < 0]:
    c = [s]
    while nx[c[-1]] >= 0: c.append(nx[c[-1]])
    if len(c) > 100:
        for p, k in enumerate(c): fil[k] = len(chains); fpos[k] = p
        chains.append(c)
print("filaments:", [len(c) for c in chains], "subunits outside filaments:", (fil < 0).sum())

# hero subunit: cylinder part, facing -x, away from the seam
az = np.degrees(np.arctan2(cents[:, 1], cents[:, 0]))
cand = np.where((blk == 0) & (np.abs(cents[:, 2]-160) < 25) & (np.abs(np.abs(az)-180) < 2.5))[0]
hero = int(cand[np.argmin(np.abs(np.abs(az[cand])-180))]); print("hero subunit", hero, cents[hero].round(1))
M_hero = Mfull[hero]

# ---------- local patch in the EMD-14238 frame (exact helical lattice of 7R1C) --
mp = mrcfile.open(RAW + "emd_14238.map.gz", permissive=True); vs = float(mp.voxel_size.x)
ns = np.array([int(mp.header.nxstart), int(mp.header.nystart), int(mp.header.nzstart)]); dat = mp.data.astype(np.float32)
lo = ns*vs; hi = (ns + np.array(dat.shape[::-1]) - 1)*vs
patch_n, patch_M, patch_match = [], [], []
for n in range(-400, 401):
    X = app(Hn(n), ref_CA); c = X.mean(0)
    if np.all(c > lo+8) and np.all(c < hi-8):
        patch_n.append(n); patch_M.append(M_hero @ Hn(n) @ T(c_mono))
        patch_match.append(int(tree.query(app(M_hero, c))[1]))
patch_n = np.array(patch_n); print("patch subunits:", len(patch_n), "unique matches:", len(set(patch_match)))
# keep density only within 6 A of the patch model (+ hero), so the map outlines the fitted subunits
allX = np.vstack([app(Hn(n), np.array([a.pos.tolist() for r in ch7 for a in r])) for n in patch_n])
zz, yy, xx = np.meshgrid(*[np.arange(s_) for s_ in dat.shape], indexing="ij")
vox = np.stack([xx, yy, zz], -1).reshape(-1, 3)*vs + lo
dmin, _ = cKDTree(allX).query(vox, distance_upper_bound=6.0)
dat = np.where(np.isfinite(dmin).reshape(dat.shape), dat, np.minimum(dat, 0)).astype(np.float32)
# write map with explicit origin so every reader places it identically
with mrcfile.new(OUT + "emd_14238.mrc", overwrite=True) as o:
    o.set_data(dat); o.voxel_size = vs
    o.header.origin.x, o.header.origin.y, o.header.origin.z = lo
    o.header.nxstart = o.header.nystart = o.header.nzstart = 0
# threshold suggestion: density at atoms of the patch model
Xall = np.vstack([app(Hn(n), ref_CA) for n in patch_n]); ijk = np.round(Xall/vs - ns).astype(int)
ok = np.all((ijk >= 0) & (ijk < 128), 1); dv = dat[ijk[ok, 2], ijk[ok, 1], ijk[ok, 0]]
print("density at CA: 10th pct", np.percentile(dv, 10).round(2), "median", np.median(dv).round(2), "map std", dat.std().round(2))

# ---------- GvpC (8GBS) -----------------------------------------------------------
s8 = gemmi.read_structure(RAW + "8GBS.cif")
A = np.vstack([[r["CA"][0].pos.tolist() for r in s8[0][c]] for c in ["A1", "A2", "A3", "A4"]])
chC = s8[0]["C"]; Cxyz = np.array([a.pos.tolist() for r in chC for a in r]); cC = Cxyz.mean(0)
gc = gemmi.read_structure(RAW + "8GBS.cif")
for c in ["A1", "A2", "A3", "A4"]: gc[0].remove_chain(c)
for r in gc[0]["C"]:
    r.name = "ALA"                                   # poly-UNK backbone+CB model -> Ala so it is treated as protein
    for a in r: a.pos = gemmi.Position(*(np.array(a.pos.tolist()) - cC)); a.b_iso = 0
gc.setup_entities(); gc.remove_empty_chains()
gc.make_mmcif_document().write_file(OUT + "gvpc_repeat.cif")
gvpc_M, gvpc_mol, gvpc_rms, gvpc_fil, gvpc_pos = [], [], [], [], []
mol_id = 0
for fi, ch in enumerate(chains):
    blocks = [ch[b:b+4] for b in range(0, len(ch)-3, 4)]
    # molecules of 5 repeats, 1 empty block between molecules (~83 % coverage, cf. GvpC:GvpA stoichiometry)
    for b0 in range(0, len(blocks)-4, 6):
        mol_blocks = blocks[b0:b0+5]; fits = []
        for ks in mol_blocks:
            Q = np.vstack([CA[k] for k in ks]); G = kabsch(A, Q)
            fits.append((G, np.sqrt(((app(G, A)-Q)**2).sum(1).mean()), ks))
        if max(f[1] for f in fits) > 4.0: continue                       # skip badly fitting (tip) regions
        if np.min([np.hypot(*cents[f[2][0], :2]) for f in fits]) < 90: continue   # no rods on the sharp tips
        for G, r, ks in fits:
            gvpc_M.append(G @ T(cC)); gvpc_mol.append(mol_id); gvpc_rms.append(r); gvpc_fil.append(fi); gvpc_pos.append(fpos[ks[0]])
        mol_id += 1
gvpc_M = np.array(gvpc_M); print("GvpC repeats:", len(gvpc_M), "molecules:", mol_id, "max fit rmsd", np.max(gvpc_rms).round(2))

np.savez_compressed(OUT + "instances.npz",
    shell_M=shell_M, shell_fil=fil, shell_pos=fpos, shell_cent=cents, hero=hero, hero_M=shell_M[hero],
    patch_M=np.array(patch_M), patch_n=patch_n, patch_match=np.array(patch_match), map_M=M_hero,
    gvpc_M=gvpc_M, gvpc_mol=np.array(gvpc_mol), gvpc_fil=np.array(gvpc_fil), gvpc_pos=np.array(gvpc_pos),
    fil_len=np.array([len(c) for c in chains]))
print("wrote", OUT)
