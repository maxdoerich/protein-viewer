"""GvpA x4 + GvpC (PDB 8GBS) with Coulombic surface potential -> static/charge-viewer/.

Runs inside ChimeraX (no GUI needed):
    /Applications/ChimeraX-1.12.app/Contents/bin/ChimeraX --nogui --exit --script proteins/gvpc_charge/tools/prep_viewer.py

The same figure as the vault's ChimeraX image (assets/scripts/chimerax_gvpa_gvpc_charge_2026-10-13.cxc): the four GvpA
chains A1-A4 and the GvpC chain C of PDB 8GBS, surfaces coloured by ChimeraX's Coulombic potential (AMBER ff14SB
charges, distance-dependent dielectric 4r, 1.4 A off the surface, kcal/(mol*e)).

GvpC in 8GBS is deposited as 33 x UNK (N, CA, C, O, CB only), so it has no charge of its own. Its side chains are taken
from the AlphaFold DB model AF-P09413-F1 (proteins/gv_animation/data/ana_gvpc_afdb.cif), residue by residue: the
deposited numbering 85-117 is UniProt P09413 repeat 3, each AlphaFold residue is superposed on its N, CA, C, CB in
8GBS and only the atoms beyond CB are copied. The 8GBS backbone stays as deposited.

Charges: GvpA as in the image (the deposited chains, so residues 2 and 66 are charged termini). GvpC is a piece cut
out of a 193-residue protein, so for the charge assignment it is extended by AlphaFold residues 84 and 118 (placed by
the overall superposition), which take the terminal charges; they are deleted again before anything is computed.
GvpA surfaces carry the potential of the four GvpA (as `coulombic #4/A1-A4` did), the GvpC surface that of GvpC
alone: its own charge, the same whether or not GvpC is shown. Where an AlphaFold rotamer runs into the wall (heavy atom < 2.5 A from GvpA),
ChimeraX's swapaa puts in the Dunbrack rotamer with the fewest clashes (ties: the most probable) instead.

Mol* colours a surface by atom, so each atom's B-factor is the mean potential of the surface points ChimeraX assigns
to it (atoms without surface points: their residue's mean, else 0). Writes gvpc-charge.cif and charge-viewer.json.
"""
import json
import os

import numpy as np

from chimerax.atomic import AtomicStructure
from chimerax.core.commands import run

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.normpath(os.path.join(HERE, '..', '..', '..'))
PDB = os.path.join(ROOT, 'proteins', 'gvpc_charge', 'data', '8gbs.cif')
AF = os.path.join(ROOT, 'proteins', 'gv_animation', 'data', 'ana_gvpc_afdb.cif')
OUT = os.path.join(ROOT, 'static', 'charge-viewer')
GVPA = ['A1', 'A2', 'A3', 'A4']
FIRST, LAST = 85, 117          # GvpC residues in 8GBS = UniProt P09413 repeat 3
OFFSET, DIEL = 1.4, 4.0        # ChimeraX coulombic defaults (as in the image)
ONE = {'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D', 'CYS': 'C', 'GLN': 'Q', 'GLU': 'E', 'GLY': 'G', 'HIS': 'H',
       'ILE': 'I', 'LEU': 'L', 'LYS': 'K', 'MET': 'M', 'PHE': 'F', 'PRO': 'P', 'SER': 'S', 'THR': 'T', 'TRP': 'W',
       'TYR': 'Y', 'VAL': 'V'}


def kabsch(p, q):
    """rotation R and translation t with R @ p + t ~ q (rows are points)"""
    pc, qc = p.mean(0), q.mean(0)
    u, _, vt = np.linalg.svd((p - pc).T @ (q - qc))
    d = np.sign(np.linalg.det(vt.T @ u.T))
    r = vt.T @ np.diag([1, 1, d]) @ u.T
    return r, qc - r @ pc


def opened(path):
    return [m for m in run(session, 'open "%s"' % path) if isinstance(m, AtomicStructure)][0]


def atoms_of(res):
    return {a.name: np.array(a.coord) for a in res.atoms}


# ---- read 8GBS and the AlphaFold model --------------------------------------------------------------------------
pdb = opened(PDB)
af = opened(AF)
pdb_res = {(r.chain_id, r.number): r for r in pdb.residues}
af_res = {r.number: r for r in af.residues}

gvpc = [pdb_res[('C', n)] for n in range(FIRST, LAST + 1)]
ca_pdb = np.array([r.find_atom('CA').coord for r in gvpc])
ca_af = np.array([af_res[n].find_atom('CA').coord for n in range(FIRST, LAST + 1)])
R, t = kabsch(ca_af, ca_pdb)
rmsd_ca = float(np.sqrt((((ca_af @ R.T + t) - ca_pdb) ** 2).sum(1).mean()))

# ---- one PDB file: GvpA as deposited (chains A-D), GvpC grafted (chain C -> E) ------------------------------------
lines, serial = [], 0
def atom_line(chain, resname, num, name, element, xyz):
    global serial
    serial += 1
    nm = name if len(name) == 4 else ' ' + name.ljust(3)
    return 'ATOM  %5d %-4s %3s %1s%4d    %8.3f%8.3f%8.3f  1.00  0.00          %2s' % (
        serial, nm, resname, chain, num, xyz[0], xyz[1], xyz[2], element)

PDB_CHAIN = dict(zip(GVPA, 'ABCD'))
for cid in GVPA:
    for r in sorted((r for r in pdb.residues if r.chain_id == cid), key=lambda r: r.number):
        for a in r.atoms:
            lines.append(atom_line(PDB_CHAIN[cid], r.name, r.number, a.name, a.element.name, a.coord))
    lines.append('TER')

seq, worst = '', []
for n in range(FIRST - 1, LAST + 2):
    src = af_res[n]
    if FIRST <= n <= LAST:
        dep = atoms_of(pdb_res[('C', n)])
        own = atoms_of(src)
        keys = [k for k in ('N', 'CA', 'C', 'CB') if k in dep and k in own]
        r_loc, t_loc = kabsch(np.array([own[k] for k in keys]), np.array([dep[k] for k in keys]))
        for a in src.atoms:
            xyz = dep[a.name] if a.name in dep else r_loc @ np.array(a.coord) + t_loc
            lines.append(atom_line('E', src.name, n, a.name, a.element.name, xyz))
        seq += ONE[src.name]
    else:                                   # flank: placed by the overall fit, carries the terminal charge only
        for a in src.atoms:
            lines.append(atom_line('E', src.name, n, a.name, a.element.name, R @ np.array(a.coord) + t))
lines += ['TER', 'END']
tmp = os.path.join(OUT, '_graft.pdb')
with open(tmp, 'w') as f:
    f.write('\n'.join(lines) + '\n')
run(session, 'close')

# ---- charges, surfaces, potentials ------------------------------------------------------------------------------
s = opened(tmp)
os.remove(tmp)
run(session, 'changechains #%s A,B,C,D,E A1,A2,A3,A4,C' % s.id_string)

# an AlphaFold rotamer comes from GvpC alone; where it runs into the wall (heavy atom < 2.5 A from GvpA), ChimeraX
# puts in the Dunbrack rotamer with the fewest clashes (ties: the most probable) instead (backbone untouched)
from scipy.spatial import cKDTree
BACKBONE = ['N', 'CA', 'C', 'O', 'CB']
def clashing():
    at = s.atoms
    wall = at[np.isin(at.chain_ids, GVPA)]
    side = at[(at.chain_ids == 'C') & ~np.isin(at.names, BACKBONE)]
    d, _ = cKDTree(wall.coords).query(side.coords)
    return sorted({a.residue for a, x in zip(side, d) if x < 2.5 and FIRST <= a.residue.number <= LAST},
                  key=lambda r: r.number)
replaced = clashing()
for r in replaced:
    run(session, 'swapaa #%s/C:%d %s rotLib Dunbrack criteria cp' % (s.id_string, r.number, r.name))
replaced = ['%s%d' % (ONE[r.name], r.number) for r in replaced]

# the united-atom charges `coulombic` itself assigns (hydrogens added to a copy, their charge folded into the heavy atoms)
from chimerax.add_charge import ChargeMethodArg
from chimerax.atomic import Residues
from chimerax.coulombic.coulombic import assign_charges
assign_charges(session, Residues(s.residues), None, ChargeMethodArg.default_value)
run(session, 'delete #%s/C:%d,%d' % (s.id_string, FIRST - 1, LAST + 1))     # the flanks have done their job

chains = {c: s.atoms[s.atoms.chain_ids == c] for c in GVPA + ['C']}
for c, at in chains.items():
    assert len(at.unique_residues) == (65 if c != 'C' else LAST - FIRST + 1), c
gvpa_atoms = chains['A1'] | chains['A2'] | chains['A3'] | chains['A4']

import chimerax.arrays  # noqa: F401  (lets _esp link)
from chimerax.coulombic._esp import potential_at_points

bfac = {}
stats = {}
# one surface for the four GvpA (the one Mol* draws: no faces inside the contacts between them) and one for GvpC
for c, at, spec in (('GvpA', gvpa_atoms, '/A1,A2,A3,A4'), ('GvpC', chains['C'], '/C')):
    spec = '#%s%s' % (s.id_string, spec)
    srf = run(session, 'surface %s enclose %s' % (spec, spec))[0]
    pts = srf.vertices + OFFSET * srf.normals
    val = potential_at_points(srf.scene_position.transform_points(pts), at.scene_coords,
                              np.array([a.charge for a in at], dtype=np.double), True, DIEL, os.cpu_count() or 1)
    v2a = srf.vertex_to_atom_map()
    sat = srf.atoms
    acc = np.zeros(len(sat)); cnt = np.zeros(len(sat))
    np.add.at(acc, v2a, val); np.add.at(cnt, v2a, 1)
    for a, sm, n in zip(sat, acc, cnt):
        if n:
            bfac[a] = sm / n
    for r in at.unique_residues:                          # atoms with no surface point of their own
        own = [bfac[a] for a in r.atoms if a in bfac]
        for a in r.atoms:
            if a not in bfac:
                bfac[a] = float(np.mean(own)) if own else 0.0
    stats[c] = {'min': round(float(val.min()), 2), 'mean': round(float(val.mean()), 2), 'max': round(float(val.max()), 2)}

for a in s.atoms:
    a.bfactor = float(np.clip(bfac.get(a, 0.0), -99.99, 99.99))

# ---- what is left after the swap: heavy-atom contacts of GvpC side chains with GvpA under 2.5 A ------------------
clashes = ['%s%d' % (ONE[r.name], r.number) for r in clashing()]

# ---- camera: as in the ChimeraX image, square onto the wall's outer face with GvpC upright across the four GvpA ----
ca_c = chains['C'].filter(chains['C'].names == 'CA').coords
axis = np.linalg.svd(ca_c - ca_c.mean(0))[2][0]                 # GvpC's helix axis
wall = gvpa_atoms.coords
normal = np.linalg.svd(wall - wall.mean(0))[2][2]               # the wall is one layer thick: its thinnest direction
normal -= axis * (normal @ axis); normal /= np.linalg.norm(normal)
if normal @ (ca_c.mean(0) - wall.mean(0)) < 0:                  # pointing out, to the side GvpC sits on
    normal = -normal
out, up = normal, axis
allc = s.atoms.coords
centre = allc.mean(0)

def net(at):
    return int(round(float(sum(a.charge for a in at))))

per_res = []
for r in chains['C'].unique_residues:
    q = float(sum(a.charge for a in r.atoms))
    if abs(q) > 0.5:
        per_res.append({'res': '%s%d' % (ONE[r.name], r.number), 'q': int(round(q))})

run(session, 'rename #%s gvpc_charge' % s.id_string)
run(session, 'save "%s" models #%s' % (os.path.join(OUT, 'gvpc-charge.cif'), s.id_string))
manifest = {
    'source': 'PDB 8GBS (GvpA x4, GvpC 85-117 backbone) + AF-P09413-F1 side chains; ChimeraX coulombic',
    'range': [-10, 10],
    'gvpc': {'residues': [FIRST, LAST], 'sequence': seq, 'net': net(chains['C']), 'charged': per_res,
             'af_rmsd_ca': round(rmsd_ca, 2), 'rotamers_replaced': replaced, 'clashes_left': clashes},
    'potential': stats,
    'centre': [round(float(x), 3) for x in centre],
    'radius': round(float(np.linalg.norm(allc - centre, axis=1).max()), 2),
    'view': {'dir': [round(float(x), 4) for x in -out], 'up': [round(float(x), 4) for x in up]}
}
with open(os.path.join(OUT, 'charge-viewer.json'), 'w') as f:
    json.dump(manifest, f, indent=1)
print('CHARGE-VIEWER', json.dumps(manifest))
