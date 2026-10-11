# Latch top-10 viewer: the ten best designs of the PC search, in 3D

One page (`content/Project/Latch Top10 Test.md`), same look and tools as the biosensor viewer (`proteins/latch_animation/`,
`static/latch-viewer/`), but instead of five steps a row of ten buttons above the stage, one per design.
Chips under the stage: **Closed**, **Open**, **Both** (open drawn faint over closed). Camera, rotation and chosen state
carry over from design to design; arrow keys step through the designs.

## Files
| what | where |
|---|---|
| viewer script, styles | `static/latch-top10/latch-top10.js`, `latch-top10.css` (load `gv-viewer.css`, then `latch-viewer.css`, first) |
| built data | `static/latch-top10/d01_closed.pdb` … `d10_open.pdb`, `designs.json` (generated, don't edit) |
| builder | `tools/prep_top10.py` |

## Rebuild
```sh
python tools/prep_top10.py "<vault>/iGEM/Wiki/Biosensing-Animation/final"   # [--n 10] [--out latch-top10]
```
A second run (blocking-peptide latches on the GvpC binder, `content/Project/Latch Top10 Blocking Peptide Test.md`) uses the
same viewer: `--out latch-top10-bp` writes its data to `static/latch-top10-bp/`, and the page's mount element sets
`data-base` (that folder) and `data-binder="GvpC binder"` (legend label).
Reads `final/best/best_candidates.csv` (ranked) and the `_closed.pdb` / `_open.pdb` pairs of the PC search
(`latch_search.py` + `finalize.py`, C:\jobs\latch_pdl1). Chains B binder, M latch, P target-binding domain, T PD-L1; the
builder writes the part into the B-factor column (0..3), which the viewer colours by.

## Honest limits
- Poly-Gly RFdiffusion backbones, no sequence; the linker of `k` GS residues is **not drawn** (k = 0 is a direct fusion).
- Hinge motions are placed by a geometric search (closed: PD-L1 overlaps the binder; open: nothing touches), not predicted.
- The score (2.20–2.27 over the top 60) and the open clearance (capped at 8 Å) saturate, so the ranking among the ten
  means little. PD-L1 stands in for the wet-lab target.
- No GvpA binding-site marking: the PC output has no site residues.
