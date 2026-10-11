#!/usr/bin/env python3
"""Build the data of the five-design latch animation (static/latch-anim5/) from the PC search results.

Each design is the GvpC repeat 3-4 (the fixed motif, Serratia GvpC residues 92-164) + a latch + a PD-L1 binder + PD-L1;
the search (latch_search.py + finalize.py) gives a closed and an open state; the block latch + domain + PD-L1 swings
rigidly between them (one hinge). The animation is the biosensor viewer's (proteins/latch_animation/), so this script does
what prep_viewer.py does, per design, with these changes:
  - the binder is placed on the vesicle wall by superposing it on GvpC 92-164 of the molecule that the Model page's
    viewer puts on the wall (gvpc_af.cif + the 8GBS-based placement): every design docks at the same site;
  - the wall is the real GvpA shell (instances.npz); the GvpA copies in contact with the binder are their own
    assembly ("gvpa"), the rest is "wall";
  - hinge, fly-in path and camera frames are computed per design.
Everything is stored in the map frame of the Model page's viewer (the vesicle's long axis is z).

  python tools/prep_anim5.py "<vault>/iGEM/Wiki/Biosensing-Animation/final" [--ranks 1 7 8 9 10]
Writes latch-assembly.cif and latch-anim.json. Requires numpy, scipy, gemmi.
"""
import argparse, csv, json, os
import numpy as np, gemmi
from scipy.spatial import cKDTree
from scipy.spatial.transform import Rotation as Rot

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.abspath(os.path.join(HERE, '..', '..', '..'))
GV = os.path.join(REPO, 'proteins', 'gv_animation', 'data')
OUT = os.path.join(REPO, 'static', 'latch-anim5')

FLY = 90.0           # A the target starts from
AWAY = 65.0          # A the construct stands off the wall until the last step
GAP = 25.0           # A of room around the model for each camera frame
CLASH_A = 2.8        # A, heavy atoms: the clash view
SITE_A = 6.0         # A, binder atoms this close to a GvpA atom: GvpA-binding site
LINKER_STEP = 3.8
CYL_Z = 240.0
STRETCH = [(0, 32), (33, 65), (66, 98), (99, 131), (131, 163)]


def kabsch(P, Q):
    pc, qc = P.mean(0), Q.mean(0)
    U, _, Vt = np.linalg.svd((P - pc).T @ (Q - qc))
    R = Vt.T @ np.diag([1, 1, np.sign(np.linalg.det(Vt.T @ U.T))]) @ U.T
    return R, qc - R @ pc


def xyz(chain):
    return np.array([a.pos.tolist() for r in chain for a in r])


def mat(R, t):
    M = np.eye(4); M[:3, :3] = R; M[:3, 3] = t
    return M


def move(chain, M):
    for r in chain:
        for a in r:
            a.pos = gemmi.Position(*(M[:3, :3] @ np.array(a.pos.tolist()) + M[:3, 3]))


def fib(n):
    k = np.arange(n) + .5
    phi = np.arccos(1 - 2 * k / n)
    th = np.pi * (1 + 5 ** .5) * k
    return np.stack([np.cos(th) * np.sin(phi), np.sin(th) * np.sin(phi), np.cos(phi)], 1)


# ---------------------------------------------------------------- the wall and the docking site (Model page frame)
def wall_and_site():
    d = np.load(os.path.join(GV, 'instances.npz'))
    F = np.linalg.inv(d['map_M'])
    shell = np.array([F @ M for M in d['shell_M']])
    gvpc = np.array([F @ M for M in d['gvpc_M']])
    af = gemmi.read_structure(os.path.join(GV, 'gvpc_af.cif'))
    af_ca = np.array([r['CA'][0].pos.tolist() for r in af[0][0]])
    rep_ca = np.array([r['CA'][0].pos.tolist() for r in gemmi.read_structure(os.path.join(GV, 'gvpc_repeat.cif'))[0][0]])
    axis_pt, axis_dir = F[:3, 3], F[:3, 2] / np.linalg.norm(F[:3, 2])
    mols = []
    for m in sorted(set(d['gvpc_mol'].tolist())):          # as gv_animation/tools/prep_viewer.py: one rigid fit per 5-repeat molecule
        ks = np.where(d['gvpc_mol'] == m)[0]
        if np.abs(d['gvpc_M'][ks][:, 2, 3]).max() > CYL_Z:
            continue
        tgt = [(gvpc[k][:3, :3] @ rep_ca.T).T + gvpc[k][:3, 3] for k in ks]
        R, t = kabsch(np.vstack([af_ca[a:b][:32] for a, b in STRETCH]), np.vstack([t_[:32] for t_ in tgt]))
        rod = af_ca @ R.T + t
        v = np.vstack(tgt) - axis_pt; rad_t = np.linalg.norm(v - np.outer(v @ axis_dir, axis_dir), axis=1).min()
        v = rod - axis_pt; rad_r = np.linalg.norm(v - np.outer(v @ axis_dir, axis_dir), axis=1).min()
        o = rod.mean(0) - axis_pt; o -= (o @ axis_dir) * axis_dir
        t = t + max(0.0, rad_t - rad_r) * o / np.linalg.norm(o)
        mols.append(mat(R, t))
    # the molecule whose repeats 3-4 lie nearest the equator: flat wall, no cone
    zs = [abs(((af_ca[91:164] @ M[:3, :3].T + M[:3, 3]) - axis_pt) @ axis_dir).mean() for M in mols]
    return dict(shell=shell, af_ca=af_ca, M=mols[int(np.argmin(zs))], axis_pt=axis_pt, axis_dir=axis_dir)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('final')
    ap.add_argument('--ranks', type=int, nargs='+', default=[1, 7, 8, 9, 10])
    a = ap.parse_args()
    src = os.path.join(a.final, 'best')
    rows = {int(r['rank']): r for r in csv.DictReader(open(os.path.join(src, 'best_candidates.csv')))}
    os.makedirs(OUT, exist_ok=True)

    W = wall_and_site()
    mono = gemmi.read_structure(os.path.join(GV, 'gvpa_monomer.cif'))
    mono_xyz = xyz(mono[0][0])
    shell = W['shell']
    cen = shell[:, :3, 3]

    st = None
    designs, report = [], []
    docked_B = None
    for n, rk in enumerate(a.ranks, 1):
        r = rows[rk]
        tag = f'd{n}'
        closed = gemmi.read_structure(os.path.join(src, r['name'] + '_closed.pdb'))
        opened = gemmi.read_structure(os.path.join(src, r['name'] + '_open.pdb'))
        C = {ch.name: xyz(ch) for ch in closed[0]}
        O = {ch.name: xyz(ch) for ch in opened[0]}
        # 1. dock: binder -> GvpC 92-164 of the chosen wall molecule
        Bca = np.array([res['CA'][0].pos.tolist() for res in closed[0]['B']])
        R0, t0 = kabsch(Bca, W['af_ca'][91:164])
        rms = float(np.sqrt(((Bca @ R0.T + t0 - W['af_ca'][91:164]) ** 2).sum(1).mean()))
        D = W['M'] @ mat(R0, t0)
        tr = lambda P: P @ D[:3, :3].T + D[:3, 3]
        C = {k: tr(v) for k, v in C.items()}
        O = {k: tr(v) for k, v in O.items()}
        for ch in closed[0]:
            move(ch, D)
        # 2. the swinging block (latch + domain + PD-L1) is rigid: closed -> open is one rotation about a hinge axis
        Pc = np.vstack([C[k] for k in 'MPT']); Po = np.vstack([O[k] for k in 'MPT'])
        Rb, tb = kabsch(Pc, Po)
        res = float(np.sqrt(((Pc @ Rb.T + tb - Po) ** 2).sum(1).mean()))
        assert res < 0.05, (r['name'], 'block is not rigid', res)
        rv = Rot.from_matrix(Rb).as_rotvec(); theta = float(np.linalg.norm(rv)); ax = rv / theta
        c0 = np.linalg.lstsq(np.eye(3) - Rb, tb, rcond=None)[0]
        # junction: the binder's terminal residue where the latch is fused
        jB = C['B'][-4:].mean(0) if r['terminus'] == 'c' else C['B'][:4].mean(0)
        c = c0 + ax * ((jB - c0) @ ax)                    # the point of the axis nearest the junction
        treeB = cKDTree(C['B'])
        def sweep(ang):                                   # smallest distance of latch (without the 10 residues at the fused end) and domain to the binder on the way; PD-L1 flies in later
            lat = C['M'][40:] if r['terminus'] == 'c' else C['M'][:-40]
            sub = np.vstack([lat[::3], C['P'][::3]])
            best = 9e9
            for k in np.linspace(0.02, 1, 50):
                Rk = Rot.from_rotvec(ax * ang * k).as_matrix()
                best = min(best, treeB.query((sub - c) @ Rk.T + c)[0].min())
            return best
        # closed -> open is fixed by the search; which way round the hinge it turns is not: take the way with more room
        # for the latch and the domain (the animation shows them swing without PD-L1)
        clear, ang = max((sweep(theta), theta), (sweep(theta - 2 * np.pi), theta - 2 * np.pi))
        # 3. the way in of PD-L1: a straight line from far away to the bound pose in the open state, clear of binder, latch, domain
        obs = np.vstack([O['B'], O['M'], O['P']])
        tobs = cKDTree(obs)
        Tsub = O['T'][::4]
        best = (-1, None)
        for u in fib(900):
            m = 9e9
            for s in np.arange(6, FLY + 1, 3.0):
                m = min(m, tobs.query(Tsub + u * s)[0].min())
                if m < best[0]:
                    break
            if m > best[0]:
                best = (m, u)
        fly_clear, u_open = best
        fly_closed = Rb.T @ u_open                        # the viewer flies in the closed frame and swings the result
        # 4. clash view: binder / PD-L1 overlap with the latch shut
        tT = cKDTree(C['T'])
        fa = [(i, at) for i, rr in enumerate(closed[0]['B']) for at in rr]
        hits = tT.query_ball_point(np.array([at.pos.tolist() for _, at in fa]), CLASH_A)
        clash_b = {i for (i, _), h in zip(fa, hits) if h}
        tB = cKDTree(C['B'])
        ta = [(i, at) for i, rr in enumerate(closed[0]['T']) for at in rr]
        thits = tB.query_ball_point(np.array([at.pos.tolist() for _, at in ta]), CLASH_A)
        clash_t = {i for (i, _), h in zip(ta, thits) if h}
        n_pairs = int(sum(len(h) for h in thits))
        # 5. against the real wall: GvpA copies near the docked binder
        near = np.where(np.linalg.norm(cen - C['B'].mean(0), axis=1) < 70)[0]
        wall_xyz = np.vstack([mono_xyz @ shell[i][:3, :3].T + shell[i][:3, 3] for i in near])
        wtree = cKDTree(wall_xyz)
        wall = {}
        for state, X in (('closed', C), ('open', O)):
            for k, name in (('M', 'latch'), ('P', 'domain'), ('T', 'PD-L1')):
                dist = wtree.query(X[k])[0]
                wall[f'{state}_{name}'] = dict(min=round(float(dist.min()), 2), n_lt3=int((dist < 3.0).sum()))
        if docked_B is None:
            docked_B = C['B'].copy()
            bd = wtree.query(C['B'])[0]
            site = {i for i, rr in enumerate(closed[0]['B']) if any(wtree.query(at.pos.tolist())[0] < SITE_A for at in rr)}
            contact = np.array([i for i in near if cKDTree(mono_xyz @ shell[i][:3, :3].T + shell[i][:3, 3]).query(docked_B)[0].min() < SITE_A])
        # closed latch over the site: how many binding-site residues has it within 4 A, closed and open?
        site_xyz = np.array([at.pos.tolist() for i, rr in enumerate(closed[0]['B']) if i in site for at in rr])
        cover = {s: int((cKDTree(X['M']).query(site_xyz)[0] < 4.0).sum()) for s, X in (('closed', C), ('open', O))}
        designs.append(dict(
            tag=tag, rank=int(r['rank']), id=r['id'], moves=r['moves'], k=int(r['k_linker']), terminus=r['terminus'],
            n_binder=int(r['n_binder']), n_latch=int(r['n_latch']), quality=float(r['quality']), clearance=float(r['clearance_open_A']),
            hinge=dict(center=c.round(3).tolist(), axis=ax.round(4).tolist(), deg=round(float(np.degrees(ang)), 2)),
            fly=dict(dir=fly_closed.round(4).tolist(), dist=FLY),
            clash=dict(binder_residues=len(clash_b), target_residues=len(clash_t), atom_pairs=n_pairs, max_A=CLASH_A),
            geometry=dict(binder_vs_AF_rms=round(rms, 2), block_rms=round(res, 4), sweep_clear=round(float(clear), 2),
                          fly_clear=round(float(fly_clear), 2), wall=wall, site_cover=cover, n_site_atoms=len(site_xyz)),
            clash_b=sorted(clash_b), clash_t=sorted(clash_t), C=C, O=O, st=closed, ang=ang, c=c, ax=ax))
        report.append((tag, r['name'], round(rms, 2), round(res, 4), round(float(np.degrees(ang)), 1), round(float(clear), 2),
                       round(float(fly_clear), 2), wall, cover))
        print(tag, r['name'], 'dockRMS %.2f blockRMS %.4f swing %.1f deg (clear %.1f A) fly clear %.1f A' %
              (rms, res, np.degrees(ang), clear, fly_clear))
        print('   wall  closed: latch<3A %d domain<3A %d PD-L1<3A %d | open: latch %d domain %d PD-L1 %d | site atoms near latch closed/open %d/%d of %d' %
              (wall['closed_latch']['n_lt3'], wall['closed_domain']['n_lt3'], wall['closed_PD-L1']['n_lt3'],
               wall['open_latch']['n_lt3'], wall['open_domain']['n_lt3'], wall['open_PD-L1']['n_lt3'],
               cover['closed'], cover['open'], len(site_xyz)))

    # ------------------------------------------------------------ one model: per design F (binder), L, P, T, J, K; plus GvpA (G)
    st = gemmi.Structure()
    st.spacegroup_hm = 'P 1'
    st.add_model(gemmi.Model('1'))
    names = {}

    def chain_of(name, residues):
        ch = gemmi.Chain(name)
        for res in residues:
            ch.add_residue(res)
        for i, rr in enumerate(ch, 1):
            rr.seqid = gemmi.SeqId(i, ' ')
        return ch

    def straight(name, a_res, b_res):
        """latch end -> domain end: a straight, schematic Gly linker (CA every ~3.8 A), as latch_animation's viewer draws it"""
        p0, p1 = np.array(a_res['CA'][0].pos.tolist()), np.array(b_res['CA'][0].pos.tolist())
        m = max(1, int(round(np.linalg.norm(p1 - p0) / LINKER_STEP)) - 1)
        u = (p1 - p0) / np.linalg.norm(p1 - p0)
        w = np.cross(u, [0, 0, 1.0]); w /= np.linalg.norm(w)
        res_ = [a_res]
        for i in range(1, m + 1):
            q = p0 + (p1 - p0) * i / (m + 1)
            rs = gemmi.Residue(); rs.name = 'GLY'
            for an, el, pos in (('N', 'N', q - 1.2 * u), ('CA', 'C', q), ('C', 'C', q + 1.2 * u), ('O', 'O', q + 1.2 * u + 1.2 * w)):
                at = gemmi.Atom(); at.name = an; at.element = gemmi.Element(el); at.pos = gemmi.Position(*pos)
                rs.add_atom(at)
            res_.append(rs)
        res_.append(b_res)
        return chain_of(name, res_)

    desc = {'F': 'GvpC repeat 3-4 (Serratia GvpC residues 92-164, the fixed motif of the design)',
            'L': 'Latch (RFdiffusion backbone, poly-Gly: no sequence yet)',
            'P': 'Target-binding domain: PD-L1 binder (177 aa)', 'T': 'Target: PD-L1 ectodomain (211 aa)',
            'J': 'Bridge: binder to latch', 'K': 'Linker: latch to target-binding domain (schematic, straight)',
            'G': 'GvpA (B. megaterium, PDB 7R1C), one copy per wall subunit'}
    for d in designs:
        s_ = d['st'][0]
        n = d['tag'][1:]
        B, M_, P_, T_ = (s_[k] for k in 'BMPT')
        bj = (B[len(B) - 1], M_[0]) if d['terminus'] == 'c' else (M_[len(M_) - 1], B[0])
        free_end = M_[len(M_) - 1] if d['terminus'] == 'c' else M_[0]            # the latch end that carries the domain
        pa = np.array(free_end['CA'][0].pos.tolist())
        ends = [P_[0], P_[len(P_) - 1]]
        dist = [np.linalg.norm(np.array(e['CA'][0].pos.tolist()) - pa) for e in ends]
        join = ends[int(np.argmin(dist))]
        d['link_dist'] = round(float(min(dist)), 2)
        for ch_name, src_chain in (('F', B), ('L', M_), ('P', P_), ('T', T_)):
            st[0].add_chain(chain_of(ch_name + n, [src_chain[i] for i in range(len(src_chain))]))
        st[0].add_chain(chain_of('J' + n, list(bj)))
        st[0].add_chain(straight('K' + n, free_end, join))
        for ch in st[0]:                                                          # B-factor: roles the viewer colours by (README)
            if not ch.name.endswith(n) or ch.name[1:] != n:
                continue
        for i, rr in enumerate(st[0]['F' + n]):
            for at in rr:
                at.b_iso = 2.0 if i in d['clash_b'] else (1.0 if i in site else 0.0)
        for i, rr in enumerate(st[0]['T' + n]):
            for at in rr:
                at.b_iso = 1.0 if i in d['clash_t'] else 0.0
        for ch_name in 'LPJK':
            for rr in st[0][ch_name + n]:
                for at in rr:
                    at.b_iso = 0.0
        for k, v in desc.items():
            if k != 'G':
                names[k + n] = v + ' (design %d, %s)' % (d['rank'], d['id'])
    mono.rename_chain(mono[0][0].name, 'G')
    for at in (a for rr in mono[0]['G'] for a in rr):
        at.b_iso = 0.0
    st[0].add_chain(mono[0]['G'])
    names['G'] = desc['G']
    for ch in st[0]:
        for rr in ch:
            rr.subchain = ch.name
            rr.entity_id = ''
    st.setup_entities()
    for e in st.entities:
        ch = next(c for c in st[0] if c[0].subchain in e.subchains)
        e.full_sequence = [rr.name for rr in ch]
    st.assign_label_seq_id()
    sub_of = {ch.name: ch[0].subchain for ch in st[0]}
    edesc = {e.name: names[ch] for e in st.entities for ch, sub in sub_of.items() if sub in e.subchains}
    doc = st.make_mmcif_document()
    blk = doc.sole_block()
    ent = [(i, t) for i, t in blk.find('_entity.', ['id', 'type'])]
    blk.find_mmcif_category('_entity.').erase()
    loop = blk.init_mmcif_loop('_entity.', ['id', 'type', 'pdbx_description'])
    for i, t in ent:
        loop.add_row([i, t, gemmi.cif.quote(edesc.get(i, '?'))])
    for cat in ('_pdbx_struct_assembly.', '_pdbx_struct_assembly_gen.', '_pdbx_struct_oper_list.', '_cell.', '_symmetry.'):
        blk.find_mmcif_category(cat).erase()
    lab = {}
    for row in blk.find('_atom_site.', ['auth_asym_id', 'label_asym_id']):
        lab.setdefault(row[0], row[1])

    # ------------------------------------------------------------ assemblies
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
            out.append(str(ids[i]) if i == j else f'{ids[i]}-{ids[j]}')
            i = j + 1
        return ','.join(out)

    def assembly(aid, details, op_ids, chains):
        gens.append((aid, '(' + expr(op_ids) + ')', ','.join(lab[c] for c in chains)))
        assemblies.append((aid, details, len(op_ids)))

    one = add_ops([np.eye(4)])
    for d in designs:
        n, t = d['tag'][1:], d['tag']
        assembly(t + '_binder', desc['F'], one, ['F' + n])
        assembly(t + '_latch', desc['L'], one, ['L' + n, 'J' + n])
        assembly(t + '_tbd', desc['P'], one, ['P' + n, 'K' + n])
        assembly(t + '_target', desc['T'], one, ['T' + n])
    contact_set = set(contact.tolist())
    rest = [i for i in range(len(shell)) if i not in contact_set]
    assembly('gvpa', 'GvpA copies in contact with the binder', add_ops([shell[i] for i in sorted(contact_set)]), ['G'])
    assembly('wall', 'wall of the gas vesicle: %d GvpA copies' % len(rest), add_ops([shell[i] for i in rest]), ['G'])
    loop = blk.init_mmcif_loop('_pdbx_struct_assembly.', ['id', 'details', 'oligomeric_details', 'oligomeric_count'])
    for aid, det, k in assemblies:
        loop.add_row([aid, gemmi.cif.quote(det), '?', str(k)])
    loop = blk.init_mmcif_loop('_pdbx_struct_assembly_gen.', ['assembly_id', 'oper_expression', 'asym_id_list'])
    for aid, e, asym in gens:
        loop.add_row([aid, gemmi.cif.quote(e), asym])
    cols = ['id', 'type'] + [f'matrix[{i}][{j}]' for i in (1, 2, 3) for j in (1, 2, 3)] + [f'vector[{i}]' for i in (1, 2, 3)]
    loop = blk.init_mmcif_loop('_pdbx_struct_oper_list.', cols)
    for i, M in enumerate(ops, 1):
        loop.add_row([str(i), gemmi.cif.quote('general operation')]
                     + [f'{M[r_, s_]:.6f}' for r_ in range(3) for s_ in range(3)] + [f'{M[r_, 3]:.4f}' for r_ in range(3)])
    path = os.path.join(OUT, 'latch-assembly.cif')
    doc.write_file(path, gemmi.cif.Style.Pdbx)
    with open(path) as f:
        body = f.read()
    with open(path, 'w') as f:
        f.write('# Chains F (GvpC) and probably P and T (PD-L1 binder, PD-L1) are derived from AlphaFold Server predictions. If so,\n'
                '# this information is subject to the AlphaFold Server Output Terms of Use found at\n'
                '# alphafoldserver.com/output-terms (non-commercial use).\n' + body)

    # ------------------------------------------------------------ camera frames (map frame; the vesicle axis is z)
    ax_pt, ax_dir = W['axis_pt'], W['axis_dir']
    site_pt = docked_B.mean(0)
    out_ = site_pt - ax_pt
    out_ -= (out_ @ ax_dir) * ax_dir
    out_ /= np.linalg.norm(out_)
    away = out_ * AWAY
    near_c = cen[np.linalg.norm(cen - site_pt, axis=1) < 60]

    def sphere(P):
        c = (P.min(0) + P.max(0)) / 2
        return c, float(np.linalg.norm(P - c, axis=1).max() + GAP)

    def frame(P, dir_, up):
        c, rad = sphere(P)
        return dict(target=c.round(2).tolist(), radius=round(rad, 1), dir=dir_.round(4).tolist(), up=up.round(4).tolist())

    out_designs = []
    for d in designs:
        C, ang, c, ax = d['C'], d['ang'], d['c'], d['ax']
        Rs = Rot.from_rotvec(ax * ang).as_matrix()
        swing = lambda P: (P - c) @ Rs.T + c
        n = d['tag'][1:]
        P_b = C['B']
        P_l = np.vstack([C['M']])
        P_t = C['P']
        P_T = C['T']
        fly = np.array(d['fly']['dir'])
        everything = np.vstack([P_b, P_l, swing(P_l), P_t, swing(P_t), P_T, swing(P_T)])
        Vt = np.linalg.svd(everything - everything.mean(0), full_matrices=False)[2]
        vdir, vup = Vt[2], Vt[0]
        fr = {
            'binder': frame(P_b + away, vdir, vup),
            'latch': frame(np.vstack([P_b, P_l]) + away, vdir, vup),
            'tbd': frame(np.vstack([P_b, P_l, P_t]) + away, vdir, vup),
            'mechanism': frame(np.vstack([P_b, P_l, swing(P_l), P_t, swing(P_t), swing(P_T), swing(P_T + fly * FLY * 0.5)]) + away, vdir, vup),
        }
        v = P_T.mean(0) - P_b.mean(0); v /= np.linalg.norm(v)
        both = np.vstack([P_b, P_T])
        flat = both - np.outer(both @ v, v)
        w = np.linalg.svd(flat - flat.mean(0), full_matrices=False)[2]
        side = w[1] - (w[1] @ v) * v; side /= np.linalg.norm(side)
        fr['clash'] = frame(both + away, side, v)
        docked = np.vstack([P_b, swing(P_l), swing(P_t), swing(P_T)])
        fr['dock'] = frame(np.vstack([docked, docked + away, near_c]), -out_, ax_dir)
        o = dict(tag=d['tag'], rank=d['rank'], id=d['id'], moves=d['moves'], k=d['k'], terminus=d['terminus'],
                 n_binder=d['n_binder'], n_latch=d['n_latch'], quality=d['quality'], clearance=d['clearance'],
                 hinge=d['hinge'], fly=d['fly'], clash=d['clash'], geometry=d['geometry'], link_dist=d['link_dist'], frames=fr)
        out_designs.append(o)
    manifest = dict(frame='map frame of the Model page viewer (EMD-14238 / 7R1C); the construct stands off by `away` until the last step',
                    away=away.round(2).tolist(), n_site_residues=len(site), designs=out_designs,
                    counts={a: k for a, _, k in assemblies})
    with open(os.path.join(OUT, 'latch-anim.json'), 'w') as f:
        json.dump(manifest, f, indent=1)
    for fn in ('latch-assembly.cif', 'latch-anim.json'):
        print(f'{fn:20s} {os.path.getsize(os.path.join(OUT, fn)) / 1e3:8.1f} kB')
    print('assemblies:', manifest['counts'])
    print('away points %+.2f along the wall normal (should be > 0); site residues %d of %d; contact GvpA copies %d' %
          (1.0, len(site), len(docked_B), len(contact)))
    return manifest


if __name__ == '__main__':
    main()
