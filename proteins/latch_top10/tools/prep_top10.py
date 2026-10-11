#!/usr/bin/env python3
"""Build the data of the latch top-10 viewer (static/latch-top10/) from the PC search results.

Input:  <results>/final/best/best_candidates.csv + best/<name>_closed.pdb / _open.pdb
        (output of latch_search.py + finalize.py: iGEM/Wiki/Biosensing-Animation/final/ in the vault)
Output: d01_closed.pdb, d01_open.pdb ... d10_*.pdb  (ATOM lines only, B-factor column = role: 0 binder,
        1 latch, 2 target-binding domain, 3 target) and designs.json (card text + one camera frame per design).

  python tools/prep_top10.py "<vault>/iGEM/Wiki/Biosensing-Animation/final" [--n 10] [--out latch-top10]
"""
import argparse, csv, json, pathlib
import numpy as np

ROLE = {'B': 0, 'L': 1, 'M': 1, 'P': 2, 'T': 3}          # chain -> role (latch_core.py: B binder, L/M latch (fixed/moving part), P domain, T PD-L1)
STATIC = pathlib.Path(__file__).resolve().parents[3] / 'static'


def read(path):
    lines, xyz = [], []
    for l in open(path):
        if l.startswith('ATOM'):
            r = ROLE[l[21]]
            lines.append(l[:60] + f'{float(r):6.2f}' + l[66:].rstrip('\n') + '\n')
            xyz.append((float(l[30:38]), float(l[38:46]), float(l[46:54])))
    return lines, np.array(xyz)


def frames(Xc, Xo):
    """Camera frames: look along the thinnest direction of both states together, longest direction across the
    screen; per view (closed, open, both) the centre and the half-extents across (h) and up (v) the screen."""
    allX = np.vstack([Xc, Xo])
    _, _, vt = np.linalg.svd(allX - allX.mean(0), full_matrices=False)
    right, up, d = vt
    if d[2] < 0:                                   # a stable sign, so repeated builds give the same file
        d = -d
    out = dict(dir=[round(float(v), 4) for v in d], up=[round(float(v), 4) for v in up])
    for key, X in (('closed', Xc), ('open', Xo), ('both', allX)):
        P = X @ np.stack([right, up, d]).T         # coordinates along the three axes
        lo, hi = P.min(0), P.max(0)
        mid = (lo + hi) / 2
        out[key] = dict(target=[round(float(v), 2) for v in mid @ np.stack([right, up, d])],
                        h=round(float(hi[0] - lo[0]) / 2, 1), v=round(float(hi[1] - lo[1]) / 2, 1),
                        z=round(float(hi[2] - lo[2]) / 2, 1))
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('final')
    ap.add_argument('--n', type=int, default=10)
    ap.add_argument('--out', default='latch-top10', help='folder under static/ that receives the data')
    a = ap.parse_args()
    OUT = STATIC / a.out
    src = pathlib.Path(a.final) / 'best'
    rows = list(csv.DictReader(open(src / 'best_candidates.csv')))[:a.n]
    OUT.mkdir(parents=True, exist_ok=True)
    designs = []
    for i, r in enumerate(rows, 1):
        tag = f'd{i:02d}'
        closed, Xc = read(src / f"{r['name']}_closed.pdb")
        opened, Xo = read(src / f"{r['name']}_open.pdb")
        assert len(closed) == len(opened), r['name']
        open(OUT / f'{tag}_closed.pdb', 'w').writelines(closed + ['END\n'])
        open(OUT / f'{tag}_open.pdb', 'w').writelines(opened + ['END\n'])
        still = np.abs(Xc - Xo).max(1) < 1e-2      # atoms that do not move between the states
        designs.append(dict(
            tag=tag, rank=int(r['rank']), id=r['id'], moves=r['moves'], k=int(r['k_linker']), terminus=r['terminus'],
            n_binder=int(r['n_binder']), n_latch=int(r['n_latch']), quality=float(r['quality']),
            clearance=float(r['clearance_open_A']), still=int(still.sum()), atoms=len(Xc),
            frame=frames(Xc, Xo)))
        print(tag, r['name'], 'atoms', len(Xc), 'fixed', int(still.sum()), designs[-1]['frame']['both'])
    json.dump(dict(roles=['binder', 'latch', 'tbd', 'target'], designs=designs), open(OUT / 'designs.json', 'w'), indent=1)


if __name__ == '__main__':
    main()
