# Latch animation, five designs: the biosensor viewer's five steps for five designs of the PC search

Page `content/Project/Latch Animations Test.md`. Five buttons above the model (ranks 1, 7, 8, 9, 10 of the search's top 10:
all single-hinge designs), under it the five steps of `proteins/latch_animation/`: GvpC repeat → + latch → + target-binding
domain → + target (latch opens, PD-L1 flies in; chip *Bind with latch closed*) → + gas vesicle (docks onto the wall).
Another design starts again at step 1.

## Files
| what | where |
|---|---|
| viewer script, styles | `static/latch-anim5/latch-anim.js` (a patched copy of `latch-viewer.js`), `latch-anim.css` (load `gv-viewer.css`, `latch-viewer.css` first) |
| built data | `static/latch-anim5/latch-assembly.cif`, `latch-anim.json` (generated, don't edit) |
| builder | `tools/prep_anim5.py` |

## Rebuild
```sh
python tools/prep_anim5.py "<vault>/iGEM/Wiki/Biosensing-Animation/final" [--ranks 1 7 8 9 10]
```
Needs `proteins/gv_animation/data/` (instances.npz, gvpc_af.cif, gvpc_repeat.cif, gvpa_monomer.cif) and numpy, scipy, gemmi.
Only single-hinge designs work: the swing is one rigid rotation of latch + domain + PD-L1 (the script asserts it).

## What it does differently from `latch_animation`
- **The binder is GvpC.** These designs are the GvpC cage–latch designs: chain B is the fixed motif, Serratia GvpC residues
  92–164 (repeat 3–4; fits `gvpc_af.cif` with 0.2–0.4 Å rmsd). It is docked by superposing it on those residues of the
  GvpC molecule the Model page's viewer puts on the wall (equator, flat part), so all designs sit at the same site.
- **The wall is the real GvpA shell.** The six GvpA copies within 6 Å of the repeat are their own assembly (`gvpa`).
- **Hinge** = rotation closed → open of the block M+P+T (Kabsch, rms 0.001 Å), centre = the axis point nearest the junction.
  The search fixes only the two end states; the script takes the way round the axis with more room for latch and domain.
- **Way in of PD-L1**: straight line in the open state, ≥ 2.5 Å from binder, latch and domain (900 directions tried).
- Per design the card says (from `geometry` in the json) how many latch / domain / PD-L1 atoms lie closer than 3 Å to GvpA
  with the latch shut and open: the search never saw the wall.

## Honest limits
- Poly-Gly RFdiffusion backbones, no sequence; the k-residue GS linker latch → domain is drawn as a straight line.
- Swing clearance (latch beyond the first 10 residues + domain vs. binder, on the way): 0.7 Å for S0139 and 1.4 Å for
  S0161, so there the latch or domain grazes the repeat in the animation; 4.8–5.8 Å for the others (see `geometry.sweep_clear`).
  The swing is **not** checked against the wall.
- The GvpC site is the Model page's illustrative placement (8GBS layout fitted with an AlphaFold model), not a measured pose.
- Chains F, P, T are probably AlphaFold Server output (see the header of the cif and `NOTICE.md`); confirm, as for `latch_animation`.
