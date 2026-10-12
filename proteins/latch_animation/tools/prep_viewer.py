"""Build the data for the interactive Mol* viewer of the latch (biosensor) mechanism.

Usage:
  python prep_viewer.py --trim <film_dir>   once: copy the inputs from the Drylab film folder
                                            (latch_anim/rank016/film/ in the Second Brain Vault) into ../data/
  python prep_viewer.py [out_dir]           build the viewer data from ../data/ and the gas vesicle data
                                            in ../../gv_animation/data/; out_dir defaults to
                                            static/latch-viewer/ in the wiki repository
Writes:
  latch-assembly.cif  every part of the construct (binder, latch, linker + target-binding domain, target, GvpA
                      pentamer) plus one GvpA monomer with an mmCIF assembly operator per copy in the wall.
                      Each part is its own assembly, so latch-viewer.js can move and fade it on its own. The
                      polypeptide also comes as unbroken chains (binder + latch, binder + latch + linker +
                      domain, and latch + linker + domain from the hinge), so that it can be drawn as ONE cartoon
                      wherever its parts sit together; the B-factor column says which part a residue belongs to.
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
LINKER_LENGTHS = (9, 10, 11)  # residues tried for the linker loop (its ends are 25 A apart: fewer would have to run stretched)
LINKER_CLEAR = 3.5     # A it keeps from every atom it can meet (binder through the whole swing, pentamer, target, its own ends)
LINKER_TRIES = 30      # random starts per length (~2.5 s each); seeded, so a rebuild gives the same loop


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

# ---- the linker loop: a real peptide backbone, closed onto the domain ---------------------------------------------
# ideal backbone geometry (Engh & Huber)
B_NCA, B_CAC, B_CN, B_CO = 1.458, 1.525, 1.329, 1.231
A_NCAC, A_CACN, A_CNCA = 111.2, 116.2, 121.7
# where loop glycines sit in the Ramachandran plot (phi, psi, weight): polyproline II, beta, right- and left-handed
# helix, and the extended left-handed region only glycine reaches
GLY_BASINS = [(-75, 145, .30), (-120, 130, .20), (-65, -40, .20), (80, 10, .15), (90, -170, .10), (-90, 0, .05)]


def nerf(a, b, c, bond, angle, torsion):
    """the atom after a-b-c at `bond` from c, angle b-c-d = `angle`, dihedral a-b-c-d = `torsion` (degrees)"""
    bc = (c - b) / np.linalg.norm(c - b)
    n = np.cross(b - a, bc); n /= np.linalg.norm(n)
    m = np.cross(n, bc)
    th, ph = np.radians(angle), np.radians(torsion)
    d2 = np.array([-bond * np.cos(th), bond * np.sin(th) * np.cos(ph), bond * np.sin(th) * np.sin(ph)])
    return c + d2[0] * bc + d2[1] * m + d2[2] * n


def dihedral(a, b, c, d):
    b0, b1, b2 = a - b, c - b, d - c
    b1 /= np.linalg.norm(b1)
    v, w = b0 - (b0 @ b1) * b1, b2 - (b2 @ b1) * b1
    return np.degrees(np.arctan2(np.cross(b1, v) @ w, v @ w))


def gly_allowed(phi, psi):
    """coarse glycine Ramachandran: everything except the squeezed bands round phi = 0"""
    return abs(phi) > 45 and not (abs(phi) < 70 and abs(psi) < 35)


def linker(name, before, after):
    """latch C terminus -> first residue of the domain: a Gly loop built as a peptide, the way loop-modelling tools do.
    Residues are added to the latch's last residue with ideal bond lengths, angles and planar peptide bonds, their phi/psi
    drawn from where loop glycines sit (GLY_BASINS); a least-squares fit over all phi/psi at once then closes the chain onto the
    domain's first residue (N, CA, C within 0.15 A), each residue held softly in its basin. Of LINKER_TRIES starts per length, the builder keeps
    loops that closed with every residue still in an allowed region of the plot, that do not touch themselves, and that
    keep LINKER_CLEAR from everything they move past (the binder at every angle of the block's swing, the pentamer once
    docked, the target on its domain, latch and domain); of those, the one with the most room and the most typical angles.
    It moves with the block, so it is checked in the block's frame. Nothing is known about how the real linker folds (the
    real construct has a flexible one; the AlphaFold pilot used GGSGGS): this is one plausible conformation, drawn."""
    from scipy.spatial import cKDTree
    rng = np.random.default_rng(7)
    at_ = lambda r, n: np.array(r[n][0].pos.tolist())
    a, b = before[len(before) - 1], after[0]
    Na, CAa, Ca = at_(a, "N"), at_(a, "CA"), at_(a, "C")
    psi_a = dihedral(Na, CAa, Ca, at_(a, "O")) - 180.0           # the latch's own last psi: its O sits anti to the next N
    target = np.array([at_(b, "N"), at_(b, "CA"), at_(b, "C")])

    ends = lambda chain, h, t: np.array([q.pos.tolist() for r in list(chain)[h:len(chain) - t] for q in r])
    own = cKDTree(np.vstack([ends(before, 0, 2), ends(after, 2, 0), xyz(st[0]["T"])]))
    binder = cKDTree(xyz(st[0]["F"]))
    penta = cKDTree(np.vstack([xyz(st[0][c]) for c in "ABCDE"]))
    swings = [rot(hinge_ax, hinge_deg * f) for f in np.linspace(0, 1, 13)]

    def build(phis, psis, phi_b):
        """backbone N, CA, C of n new residues and of the virtual residue that should coincide with the domain's first"""
        X = [Na, CAa, Ca]
        prev_psi = psi_a
        for ph, ps in list(zip(phis, psis)) + [(phi_b, None)]:
            N = nerf(X[-3], X[-2], X[-1], B_CN, A_CACN, prev_psi)
            CA = nerf(X[-2], X[-1], N, B_NCA, A_CNCA, 180.0)
            C = nerf(X[-1], N, CA, B_CAC, A_NCAC, ph)
            X += [N, CA, C]
            prev_psi = ps
        return np.array(X[3:])                                    # 3 (n + 1) atoms: N, CA, C per residue

    CLOSE = 0.15         # A rmsd of the last N, CA, C onto the domain's first residue
    from scipy.optimize import least_squares
    wrap = lambda d: (d + 180.0) % 360.0 - 180.0

    def ccd(X, rounds):
        """cyclic coordinate descent: turn phi (N-CA) and psi (CA-C) bonds one at a time towards the target"""
        X = X.copy()
        k = len(X) // 3 - 1                                       # new residues
        for _ in range(rounds):
            for j in range(k + 1):
                for (p_, q) in ((3 * j, 3 * j + 1), (3 * j + 1, 3 * j + 2)):
                    if j == k and q != 3 * j + 1:                 # the domain's own psi is not ours to turn
                        continue
                    o, u = X[p_], X[q] - X[p_]
                    u /= np.linalg.norm(u)
                    mv = X[-3:] - o; tg = target - o
                    mv_p = mv - np.outer(mv @ u, u); tg_p = tg - np.outer(tg @ u, u)
                    th = np.arctan2((np.cross(mv_p, tg_p) @ u).sum(), (mv_p * tg_p).sum())
                    X[q + 1:] = (X[q + 1:] - o) @ rot(u, np.degrees(th)).T + o
        return X

    def close(start):
        """CCD from the drawn angles gets the chain's end to the domain (fast, ends ~1 A off); a short least-squares
        polish over all phi/psi then closes it, each angle allowed to move only a little from where CCD left it"""
        n = (len(start) - 1) // 2
        X = ccd(build(start[:n], start[n:2 * n], start[2 * n]), 40)
        if np.sqrt(((X[-3:] - target) ** 2).sum(1).mean()) > 2.5:
            return X, False
        tor = torsions(X)                                                   # (phi, psi) of the n new residues
        v0 = np.r_[[t[0] for t in tor], [t[1] for t in tor], dihedral(X[-4], X[-3], X[-2], X[-1])]   # + the closing phi
        def resid(v):
            Y = build(v[:n], v[n:2 * n], v[2 * n])
            return np.r_[(Y[-3:] - target).ravel() * 10.0, wrap(v - v0) / 90.0]
        v = least_squares(resid, v0, max_nfev=40 * len(v0)).x
        X = build(v[:n], v[n:2 * n], v[2 * n])
        return X, np.sqrt(((X[-3:] - target) ** 2).sum(1).mean()) < CLOSE

    def torsions(X):
        full = np.vstack([[Na, CAa, Ca], X])
        out = []
        for j in range(1, len(full) // 3 - 1):
            N, CA, C = full[3 * j], full[3 * j + 1], full[3 * j + 2]
            out.append((dihedral(full[3 * j - 1], N, CA, C), dihedral(N, CA, C, full[3 * j + 3])))
        return out

    def oxygens(X):
        """O of every new residue: in the peptide plane, opposite the bisector of C-CA and C-N(next)"""
        full = np.vstack([X, target[:1]])
        O = []
        for j in range(len(X) // 3 - 1):
            CA, C, Nn = full[3 * j + 1], full[3 * j + 2], full[3 * j + 3]
            v = (C - CA) / np.linalg.norm(C - CA) + (C - Nn) / np.linalg.norm(C - Nn)
            O.append(C + B_CO * v / np.linalg.norm(v))
        return np.array(O)

    def clearance(P, inner):
        c = own.query(inner)[0].min()          # latch and domain: not the residues at the two joins (bonded to them)
        for R in swings:
            c = min(c, binder.query((P - hinge_c) @ R.T + hinge_c)[0].min())
        return min(c, penta.query((P - hinge_c) @ swings[-1].T + hinge_c)[0].min())

    best, tried, closed_n = None, 0, 0
    for n in LINKER_LENGTHS:
        for _ in range(LINKER_TRIES):
            tried += 1
            pick = rng.choice(len(GLY_BASINS), size=n + 1, p=[w for *_, w in GLY_BASINS])
            pp = [(GLY_BASINS[i][0] + rng.normal(0, 15), GLY_BASINS[i][1] + rng.normal(0, 15)) for i in pick]
            X, ok = close(np.array([q[0] for q in pp[:n]] + [q[1] for q in pp[:n]] + [pp[n][0]]))
            if not ok:
                continue
            closed_n += 1
            tor = torsions(X)
            if not all(gly_allowed(ph, ps) for ph, ps in tor):
                continue
            Xn = X[:-3]                                           # the new residues (the last three are the domain's own)
            O = oxygens(X)
            P = np.vstack([Xn, O])
            res_of = np.r_[np.repeat(np.arange(n), 3), np.arange(n)]
            pairs = cKDTree(P).query_pairs(3.0, output_type="ndarray")
            if len(pairs) and (np.abs(res_of[pairs[:, 0]] - res_of[pairs[:, 1]]) >= 3).any():
                continue                                          # touches itself further along the chain
            inner = P[(res_of > 0) & (res_of < n - 1)]
            c = clearance(P, inner)
            if c < LINKER_CLEAR:
                continue
            # most room (up to 6 A counts), then the most typical angles: distance to the nearest basin centre
            typ = np.mean([min(np.hypot((ph - bp + 180) % 360 - 180, (ps - bs + 180) % 360 - 180) for bp, bs, _ in GLY_BASINS)
                           for ph, ps in tor])
            score = min(c, 6.0) - typ / 30.0
            if best is None or score > best[0]:
                best = (score, X, P, c, tor, n)
    if best is None:
        raise SystemExit("linker: no closed, allowed loop keeps %.1f A from everything (%d tries, %d closed)"
                         % (LINKER_CLEAR, tried, closed_n))
    _, X, P, c, tor, n = best
    O = oxygens(X)

    ch = gemmi.Chain(name)
    ch.add_residue(a)
    for j in range(n):
        res = gemmi.Residue(); res.name = "GLY"
        for an, el, pos in (("N", "N", X[3 * j]), ("CA", "C", X[3 * j + 1]), ("C", "C", X[3 * j + 2]), ("O", "O", O[j])):
            q = gemmi.Atom(); q.name = an; q.element = gemmi.Element(el); q.pos = gemmi.Position(*pos)
            res.add_atom(q)
        ch.add_residue(res)
    ch.add_residue(b)
    for i, r in enumerate(ch, 1):
        r.seqid = gemmi.SeqId(i, " ")
    st[0].add_chain(ch)
    return dict(residues=n, end_to_end=round(float(np.linalg.norm(at_(b, "CA") - CAa)), 1), clearance=round(float(c), 2),
                closure_rmsd=round(float(np.sqrt(((X[-3:] - target) ** 2).sum(1).mean())), 3),
                phi_psi=[[round(float(ph)), round(float(ps))] for ph, ps in tor], tried=tried, closed=closed_n)


link = linker("K", st[0]["L"], st[0]["P"])                     # latch C terminus -> target-binding domain
print("linker: %(residues)d Gly, ends %(end_to_end).1f A apart, closed to %(closure_rmsd).3f A, %(clearance).2f A from everything "
      "(%(closed)d of %(tried)d starts closed)" % link)
print("linker phi/psi:", link["phi_psi"])

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
for name in ("A", "B", "C", "D", "E", "G"):
    for r in st[0][name]:
        for a in r:
            a.b_iso = 0.0

# One unbroken cartoon. A cartoon stops at the end of its chain, so parts that are chains of their own show a seam where they
# join (a cap, a kink, a step in thickness). The viewer therefore draws the polypeptide as ONE chain wherever its parts sit
# together, and as separate pieces only while one part moves against another. Every chain below runs binder -> latch ->
# linker -> domain without a break (the joins are the 3.8 A CA-CA bonds of the design), with the part that each residue
# belongs to in the B-factor column: 0 binder, 1 GvpA-binding site, 2 clash with PD-L1 (all three binder residues, as above),
# 3 latch, 4 linker and target-binding domain. The viewer paints by it.
F_, L_, K_, P_ = (list(st[0][c]) for c in "FLKP")
hinge_res, gly = F_[-1:], K_[1:-1]                      # K runs latch end, 11 Gly, domain start: the Gly are its own
PARTS = {  # chain, what it is drawn for, residues (code: B-factor to give them; None keeps the binder's own)
    "R": ("arm", [(hinge_res, None), (L_, 3), (gly, 4), (P_, 4)]),   # latch + linker + domain, from the hinge: swings as one block
    "S": ("latch piece", [(hinge_res, None), (L_, 3)]),              # the latch alone: steps 2 and the clash view
    "U": ("domain piece", [(L_[-1:], 3), (gly, 4), (P_, 4)]),        # linker + domain: the clash view
    "V": ("binder + latch", [(F_, None), (L_, 3)]),                  # at rest after step 2
    "W": ("binder + latch + linker + domain", [(F_, None), (L_, 3), (gly, 4), (P_, 4)]),  # at rest in step 3
}
for cname, (_, parts) in PARTS.items():
    ch = gemmi.Chain(cname)
    for res, code in parts:
        for r in res:
            nr = ch.add_residue(r)
            if code is not None:
                for a in nr:
                    a.b_iso = float(code)
    for i, r in enumerate(ch, 1):
        r.seqid = gemmi.SeqId(i, " ")
    st[0].add_chain(ch)
P_latch, P_tbd = xyz(st[0]["L"]), xyz(st[0]["P"])        # for the camera frames below
for c in "JKLP":                                         # replaced by the chains above (J: not any more, K: the linker is inside them)
    if any(x.name == c for x in st[0]):
        st[0].remove_chain(c)
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
for c, (what, _) in PARTS.items():                    # the linker (11 Gly) is schematic: it is drawn, its fold is not known
    names[c] = "Unbroken cartoon, " + what + " (the linker between latch and domain is drawn, not predicted)"
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
assembly("latch", names["S"], one, "S")
assembly("tbd", names["U"], one, "U")
assembly("arm", names["R"], one, "R")
assembly("cont2", names["V"], one, "V")
assembly("cont3", names["W"], one, "W")
assembly("target", names["T"], one, "T")
assembly("gvpa", "GvpA pentamer the binder docks onto", one, "ABCDE")
# the wall in blocks of WALL_BLOCK neighbouring subunits (along each helical filament, tips last): Mol* draws a block of a
# few dozen copies faster than one object of 1725 (whose per-copy highlight data no longer fits the GPU's cache), and a
# block out of the picture is skipped as a whole. Same copies, same picture.
WALL_BLOCK = 32
order = sorted(range(len(keep)), key=lambda j: (fil[keep[j]], pos[keep[j]] < 0, pos[keep[j]]))
wall_ids = add_ops([wall_ops[j] for j in order])
wall_blocks = []
for b0 in range(0, len(wall_ids), WALL_BLOCK):
    bid = "wall-%02d" % (b0 // WALL_BLOCK)
    assembly(bid, "wall of the gas vesicle: GvpA copies %d-%d" % (b0 + 1, min(b0 + WALL_BLOCK, len(wall_ids))),
             wall_ids[b0:b0 + WALL_BLOCK], "G")
    wall_blocks.append(bid)

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
    # the PD-L1 binder and PD-L1 most likely come from AlphaFold Server predictions (file names fold_pdl1_*): confirm, see README.md
    f.write("# The target-binding domain (PD-L1 binder, in chains R, U, V, W) and PD-L1 (chain T) are probably derived from AlphaFold Server predictions. If so, this\n"
            "# information is subject to the AlphaFold Server Output Terms of Use found at\n"
            "# alphafoldserver.com/output-terms (non-commercial use).\n" + body)

# ------------------------------------------------------------------ camera frames
P_binder, P_target = (xyz(st[0][c]) for c in "FT")
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
    "wall": wall_blocks,
    "counts": {"wall copies": len(wall_ids), **{a: n for a, _, n in assemblies if not a.startswith("wall-")}},
    "geometry": {**{k: bh[k] for k in ("iface_covered_closed", "iface_covered_open", "n_iface")},
                 "placement": place["metrics"] if place else None, "linker": link},
}
with open(os.path.join(OUT, "latch-viewer.json"), "w") as f:
    json.dump(manifest, f, indent=1)

for fn in ("latch-assembly.cif", "latch-viewer.json"):
    print(f"{fn:20s} {os.path.getsize(os.path.join(OUT, fn)) / 1e3:8.1f} kB")
print("assemblies:", manifest["counts"])
