"""GvpC-repeat viewer data: gvpc_repeat_4gvpa.cif (4 GvpA, chains A-D, + one 33-residue GvpC repeat, chain F)
-> static/gvpc-viewer/gvpc-assembly.cif + gvpc-viewer.json.

gvpc-assembly.cif   the same atoms, plus two mmCIF assemblies ("gvpa" = A-D, "gvpc" = F) so gvpc-viewer.js can fade
                    them apart, and in the B-factor column of chain F the residue's job: 2 anchor, 1 other interface
                    residue, 0 the rest (the viewer paints by it). GvpA keeps B = 0.
gvpc-viewer.json    the two camera frames: the whole model, and GvpC alone seen from outside, its helix across the view.

Residue numbers below count inside the repeat (1-33); chain F is numbered 92-124 in the file (OFFSET = 91).
Needs numpy. Run from anywhere.
"""
import json, os
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
SRC = os.path.join(HERE, "..", "data", "gvpc_repeat_4gvpa.cif")
OUT = os.path.join(HERE, "..", "..", "..", "static", "gvpc-viewer")
ANCHOR = [4, 11, 12, 15, 19, 22]
INTERFACE = [1, 5, 16, 33]
OFFSET = 91
NOTICE = ("# Chain F (GvpC) is derived from an AlphaFold Server prediction. This information is subject to the\n"
          "# AlphaFold Server Output Terms of Use found at alphafoldserver.com/output-terms (non-commercial use).\n")
TAIL = """#
loop_
_pdbx_struct_assembly.id
_pdbx_struct_assembly.details
_pdbx_struct_assembly.oligomeric_details
_pdbx_struct_assembly.oligomeric_count
gvpa 'Four GvpA' ? 4
gvpc 'One GvpC repeat (33 residues)' ? 1
#
loop_
_pdbx_struct_assembly_gen.assembly_id
_pdbx_struct_assembly_gen.oper_expression
_pdbx_struct_assembly_gen.asym_id_list
gvpa (1) A,B,C,D
gvpc (1) F
#
loop_
_pdbx_struct_oper_list.id
_pdbx_struct_oper_list.type
_pdbx_struct_oper_list.matrix[1][1]
_pdbx_struct_oper_list.matrix[1][2]
_pdbx_struct_oper_list.matrix[1][3]
_pdbx_struct_oper_list.matrix[2][1]
_pdbx_struct_oper_list.matrix[2][2]
_pdbx_struct_oper_list.matrix[2][3]
_pdbx_struct_oper_list.matrix[3][1]
_pdbx_struct_oper_list.matrix[3][2]
_pdbx_struct_oper_list.matrix[3][3]
_pdbx_struct_oper_list.vector[1]
_pdbx_struct_oper_list.vector[2]
_pdbx_struct_oper_list.vector[3]
1 'identity operation' 1 0 0 0 1 0 0 0 1 0 0 0
"""

job = {OFFSET + p: 2 for p in ANCHOR}
job.update({OFFSET + p: 1 for p in INTERFACE})
lines, ca = [], {}
for l in open(SRC):
    t = l.split()
    if t and t[0] == "ATOM":
        if t[6] == "F":
            t[14] = str(job.get(int(t[8]), 0))
            l = " ".join(t) + "\n"
        if t[3] == "CA":
            ca.setdefault(t[6], []).append([float(x) for x in t[10:13]])
    lines.append(l)
os.makedirs(OUT, exist_ok=True)
open(os.path.join(OUT, "gvpc-assembly.cif"), "w").write(NOTICE + "".join(lines).rstrip("\n") + "\n" + TAIL)

# camera: from outside (away from GvpA) onto GvpC, its helix axis across the view
F = np.array(ca["F"]); A = np.vstack([np.array(ca[c]) for c in "ABCD"]); both = np.vstack([F, A])
unit = lambda v: v / np.linalg.norm(v)
axis = unit(np.linalg.svd(F - F.mean(0))[2][0])
out = F.mean(0) - A.mean(0); out = unit(out - axis * (out @ axis))           # perpendicular to the helix, away from GvpA
def frame(pts, dir, right):
    up = np.cross(right, dir); lo, hi = [], []
    for ax in (right, up):
        p = pts @ ax; lo.append(p.min()); hi.append(p.max())
    mid = right * (lo[0] + hi[0]) / 2 + up * (lo[1] + hi[1]) / 2 + dir * (pts @ dir).mean()
    r = lambda x: round(float(x), 4)
    return dict(target=[round(float(x), 2) for x in mid], dir=[r(x) for x in dir], up=[r(x) for x in up],
                halfWidth=round((hi[0] - lo[0]) / 2 + 6, 1), halfHeight=round((hi[1] - lo[1]) / 2 + 6, 1),
                radius=round(float(np.linalg.norm(pts - mid, axis=1).max()) + 6, 1))   # the camera fits the extents on screen
tilt = lambda a: unit(-out * np.cos(a) + np.cross(axis, -out) * np.sin(a))   # -out turned about the helix axis by a
json.dump({"offset": OFFSET, "anchor": ANCHOR, "interface": INTERFACE,
           "frames": {"all": frame(both, tilt(np.radians(50)), axis), "gvpc": frame(F, tilt(0), axis)}},
          open(os.path.join(OUT, "gvpc-viewer.json"), "w"), indent=1)
print("wrote", OUT)
