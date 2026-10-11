"""Place the target-binding domain + target (one rigid body) so that the latch has a job.

Usage:  python place_domain.py            writes ../data/placement.json (read by prep_viewer.py)
Requires: numpy, scipy, gemmi

Why: in the Drylab film the PD-L1 binder hangs off the latch at a pose that was only checked for "the swing frees the
site". The viewer wants more: with the latch CLOSED, PD-L1 sitting on its binder domain must clash clearly with the GvpA
binder (that is the mechanism: PD-L1 can only bind when the latch is open); with the latch OPEN, nothing may clash.
The domain is joined to the latch by a drawn, schematic linker (prep_viewer.py), so its pose relative to the latch is free:
we look for a rigid transform of (domain, PD-L1) that fulfils, in the film frame:

  closed  PD-L1 vs binder: >= CLASH_MIN atoms closer than 2.2 A (and not absurdly deep: <= CLASH_MAX closer than 3.0 A)
          domain, PD-L1, linker vs binder, latch: >= 3 A (the clash is PD-L1 vs binder only)
  swing   latch + linker + domain (PD-L1 is not there yet) swing about block_hinge: >= 2.5 A to the binder all along
  open    PD-L1 vs binder, GvpA pentamer: >= 3.5 A; domain vs pentamer >= 3 A; everything outside the vesicle wall
  fly     PD-L1 can come in to its bound pose from outside along a straight line (open frame), clear of everything

The transform with the largest smallest clearance wins. Seeded: the result is reproducible.
"""
import json, os
import numpy as np, gemmi
from scipy.spatial import cKDTree

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "..", "data")
TRIM = 26                  # keep in sync with prep_viewer.py
CLASH_MIN, CLASH_MAX = 40, 450
LINK_MIN, LINK_MAX = 14.0, 60.0   # A between the latch end and the domain's first residue (a straight drawn linker)
FLY_DIST = 90.0
MARGIN = 0.7               # the search asks for this much more than the limits above, so the result is not a razor's edge
rng = np.random.default_rng(20261007)


def read(fn):
    return gemmi.read_structure(os.path.join(DATA, fn))[0]


xyz = lambda ch, lo=0: np.array([a.pos.tolist() for r in list(ch)[lo:] for a in r])
ca = lambda ch, lo=0: np.array([r["CA"][0].pos.tolist() for r in list(ch)[lo:]])
B, L, P, T, G = read("binder.pdb")[0], read("latch.pdb")[0], read("tbd.pdb")[0], read("target.pdb")[0], read("gvpa.pdb")
pb, pl, pd_, pt = xyz(B), xyz(L), xyz(P, TRIM), xyz(T)
pg = np.vstack([xyz(c) for c in G])
lend = ca(L)[-1]                                    # latch C terminus (fixed: the latch does not move)
first0 = ca(P, TRIM)[0]                             # first residue of the domain before the move
geo = json.load(open(os.path.join(DATA, "geometry.json")))["block_hinge"]
hc, ax = np.array(geo["center"]), np.array(geo["axis"]) / np.linalg.norm(geo["axis"])


def rot_axis(deg):
    k = np.radians(deg)
    K = np.array([[0, -ax[2], ax[1]], [ax[2], 0, -ax[0]], [-ax[1], ax[0], 0]])
    return np.eye(3) + np.sin(k) * K + (1 - np.cos(k)) * K @ K


R_open = rot_axis(geo["deg"])
swing = lambda X, d: (X - hc) @ rot_axis(d).T + hc
tree = {k: cKDTree(v) for k, v in dict(b=pb, l=pl, g=pg).items()}
l_inner = cKDTree(xyz(L)[:-12])                     # the latch without its last residues (the linker starts there)


def dmin(tr, X):
    return float(tr.query(X)[0].min())


def rand_rot():
    q = rng.normal(size=4); q /= np.linalg.norm(q)
    w, x, y, z = q
    return np.array([[1 - 2 * (y * y + z * z), 2 * (x * y - z * w), 2 * (x * z + y * w)],
                     [2 * (x * y + z * w), 1 - 2 * (x * x + z * z), 2 * (y * z - x * w)],
                     [2 * (x * z - y * w), 2 * (y * z + x * w), 1 - 2 * (x * x + y * y)]])


body = np.vstack([pd_, pt])
cb = body.mean(0)
lo, hi = pb.min(0) - 45, pb.max(0) + 45


def linker_pts(f):
    return lend + (f - lend) * np.linspace(0.12, 0.88, 7)[:, None]


def evaluate(R, p, fine=False):
    """(penalty, details): penalty 0 = every constraint met. `fine` uses all atoms and a 10-degree swing; the search
    uses every second atom and a 25-degree swing"""
    st_ = 1 if fine else 2
    f = (first0 - cb) @ R.T + p
    dom, tar = ((pd_ - cb) @ R.T + p)[::st_], ((pt - cb) @ R.T + p)[::st_]
    d_link = float(np.linalg.norm(f - lend))
    pen = max(0.0, LINK_MIN - d_link) ** 2 + max(0.0, d_link - LINK_MAX) ** 2
    n22 = sum(1 for x in tree["b"].query_ball_point(tar, 2.2) if x) * st_
    n30 = sum(1 for x in tree["b"].query_ball_point(tar, 3.0) if x) * st_
    pen += (max(0, CLASH_MIN - n22) / 8.0) ** 2 + (max(0, n30 - CLASH_MAX) / 40.0) ** 2
    lk = linker_pts(f)
    M = MARGIN
    cl = dict(dom_b=dmin(tree["b"], dom) - 3.0 - M, dom_l=dmin(l_inner, dom) - 3.0 - M, tar_l=dmin(l_inner, tar) - 3.0 - M,
              lk_b=dmin(tree["b"], lk) - 3.0 - M, lk_l=dmin(l_inner, lk) - 2.5 - M)
    step = 10 if fine else 25
    both = np.vstack([dom, lk])
    sw = min(dmin(tree["b"], swing(both, a)) for a in range(step, 121, step)) - 2.5 - M
    tar_o, dom_o = swing(tar, geo["deg"]), swing(dom, geo["deg"])
    op = dict(tar_b=dmin(tree["b"], tar_o) - 3.5 - M, tar_g=dmin(tree["g"], tar_o) - 3.5 - M, dom_g=dmin(tree["g"], dom_o) - 3.0 - M,
              dom_b=dmin(tree["b"], dom_o) - 3.0 - M)
    margins = list(cl.values()) + [sw] + list(op.values())
    pen += sum(min(0.0, m) ** 2 for m in margins)
    return pen, dict(n22=n22, n30=n30, link=round(d_link, 1), clear=round(min(margins), 2), swing=round(sw, 2),
                     open=round(min(op.values()), 2), closed=round(min(cl.values()), 2))


FIB = []
for _i in range(72):
    _z = 1 - 2 * (_i + .5) / 72; _r = np.sqrt(1 - _z * _z); _a = _i * np.pi * (3 - np.sqrt(5))
    FIB.append((_r * np.cos(_a), _r * np.sin(_a), _z))
FIB = np.array(FIB)
FLY_MIN = 2.5               # A of clearance along the way in
DOCK = 12                   # A: the last stretch of the way in is the docking onto the domain and is not checked against it


def fly_clear(R, p, coarse=True):
    """(clearance, direction in the open frame): the straight way in with the most room for PD-L1"""
    dom, tar = (pd_ - cb) @ R.T + p, (pt - cb) @ R.T + p
    tar_o, dom_o = swing(tar, geo["deg"]), swing(dom, geo["deg"])
    obst = [(cKDTree(pb), 6), (cKDTree(pg), 6), (cKDTree(swing(xyz(L), geo["deg"])), 6), (cKDTree(dom_o), DOCK)]
    dirs = FIB
    if not coarse:
        n = 240
        dirs = np.array([(np.sqrt(1 - (1 - 2 * (i + .5) / n) ** 2) * np.cos(i * np.pi * (3 - np.sqrt(5))),
                          np.sqrt(1 - (1 - 2 * (i + .5) / n) ** 2) * np.sin(i * np.pi * (3 - np.sqrt(5))),
                          1 - 2 * (i + .5) / n) for i in range(n)])
    sub = tar_o[::4 if coarse else 2]
    best = (-1e9, None)
    for d in dirs:
        # the domain is the thing PD-L1 docks onto: its last DOCK A are the docking itself and are not checked
        m = min(min(float(t.query(sub + d * s_)[0].min()) for s_ in np.arange(s0, FLY_DIST + 1, 6)) for t, s0 in obst)
        if m > best[0]:
            best = (m, d)
    return best


def small_rot(k):
    U, _, Vt = np.linalg.svd(np.eye(3) * (1 - k) + rand_rot() * k)
    return U @ Vt


def restart(seed):
    global rng
    rng = np.random.default_rng(seed)
    R, p = rand_rot(), rng.uniform(pb.min(0) - 10, pb.max(0) + 10)
    cur, _ = evaluate(R, p)
    step_r, step_t = 0.35, 6.0
    for it in range(300):
        R2 = small_rot(step_r * rng.random()) @ R
        p2 = p + rng.normal(size=3) * step_t * rng.random()
        v, _ = evaluate(R2, p2)
        if v < cur:
            R, p, cur = R2, p2, v
        if it % 70 == 69:
            step_r *= 0.6; step_t *= 0.6
        if cur == 0 or (it == 130 and cur > 40):
            break
    if cur > 0:
        return seed, cur, R, p, None, None
    # phase 2: everything else holds; now also make room for PD-L1's way in
    def total(R_, p_):
        v_, det_ = evaluate(R_, p_)
        fc_ = fly_clear(R_, p_)[0]
        return v_ + max(0.0, FLY_MIN + MARGIN - fc_) ** 2, det_, fc_
    cur, det, fc = total(R, p)
    step_r, step_t = 0.12, 2.5
    for it in range(90):
        R2 = small_rot(step_r * rng.random()) @ R
        p2 = p + rng.normal(size=3) * step_t * rng.random()
        v, d2, f2 = total(R2, p2)
        if v < cur:
            R, p, cur, det, fc = R2, p2, v, d2, f2
        if cur == 0:
            break
    v, det = evaluate(R, p, fine=True)
    fl = fly_clear(R, p, coarse=False)
    ok = v == 0 and fl[0] >= FLY_MIN
    return seed, (0 if ok else 1), R, p, det, fl


if __name__ == "__main__":
    import multiprocessing as mp
    N = 640
    feas = []
    with mp.Pool() as pool:
        for seed, v, R, p, det, fl in pool.imap_unordered(restart, [20261007 + i for i in range(N)]):
            if v == 0 and det is not None:
                feas.append((min(det["clear"], fl[0] - FLY_MIN), R, p, det, fl))
                print("seed %d: %s fly clearance %.1f A" % (seed, det, fl[0]), flush=True)
    print("%d of %d restarts ended feasible (fine check)" % (len(feas), N))
    if not feas:
        raise SystemExit("no pose found: relax the constraints")
    feas.sort(key=lambda c: -c[0])
    best = feas[0]
    best = (best[0], best[1], best[2], best[3], best[4])
    _, R, p, det, fl = best
    bestf = fl
    fly_open = bestf[1]
    fly_closed = R_open.T @ fly_open                    # the viewer flies in the closed frame, then swings the block
    print("chosen: %s; fly-in direction (open frame) %s, clearance %.1f A" % (det, fly_open.round(3).tolist(), bestf[0]))
    M_ = R
    t = p - cb @ R.T
    out = dict(R=M_.tolist(), t=t.tolist(), fly_open=fly_open.tolist(), fly_closed=fly_closed.tolist(),
               fly_clearance=round(bestf[0], 2), metrics=det, trim=TRIM)
    json.dump(out, open(os.path.join(DATA, "placement.json"), "w"), indent=1)
    print("wrote placement.json")
