# Latch (biosensor) viewer: interactive Mol* figure

Five steps, one construct: **GvpA binder** → **+ latch** (blocks the binding site) → **+ target-binding domain** →
**+ target** (the latch opens first, then PD-L1 binds the open state) → **+ gas vesicle** (the binder docks onto GvpA).
Step 4 has a chip, **Bind with latch closed**: the block swings shut with PD-L1 on its domain and the residues of binder
and PD-L1 that overlap are painted reddish orange (the cartoon stays a cartoon), to look at the clash from every side.
Same look and tools as the gas-vesicle viewer (`proteins/gv_animation/`, `static/gv-viewer/`).

## Files
| what | where |
|---|---|
| viewer script, colours | `static/latch-viewer/latch-viewer.js`, `latch-viewer.css` (layout comes from `gv-viewer.css`, load it first) |
| built data | `static/latch-viewer/latch-assembly.cif`, `latch-viewer.json` (generated, don't edit) |
| inputs | `data/*.pdb`, `data/geometry.json` (trimmed copies of the Drylab film folder) |
| builder | `tools/prep_viewer.py` |
| domain placement (search) | `tools/place_domain.py` → `data/placement.json` (read by the builder) |

## Rebuild
```sh
python tools/prep_viewer.py --trim "<vault>/iGEM/Drylab-Binder-Design/Blocking-Peptide-Target-Domain/Scripts/latch_anim/rank016/film"
python tools/place_domain.py         # only if the constraints or the inputs change: ~15 min on 8 cores, writes data/placement.json
python tools/prep_viewer.py          # ~4 min (linker search); needs numpy + scipy + gemmi; also reads ../gv_animation/data/ for the wall
```
`--trim` is only needed when the Drylab pipeline produced new geometry (new binder, other target). It keeps the closed latch
(last trajectory frame) and the first model of the others; everything else in the film folder is not used.

## Put it on a page
Code stays in `static/`; the page loads it with a few lines of `::js`, like the Model page's vesicle viewer:

```
::custom id="latch-viewer"
<figure class="figure-block">
  <div class="figure-frame">
    <span class="figure-ring" aria-hidden="true"></span><span class="figure-title">3D model</span>
    <div data-latch-viewer role="group" aria-label="Interactive 3D model: how the latch construct senses a target">
      <p>An interactive 3D model of the biosensor appears here. It needs JavaScript and WebGL.</p>
    </div>
  </div>
  <figcaption class="figure-caption"><strong>Figure N: …</strong></figcaption>
</figure>
::js
(function () {
  var t = document.querySelector('link[href*="tokens.css"]');
  if (!t) return;
  var base = t.href.replace(/tokens\.css.*$/, '');
  ['gv-viewer/gv-viewer.css', 'latch-viewer/latch-viewer.css'].forEach(function (f) {
    var l = document.createElement('link');
    l.rel = 'stylesheet'; l.href = base + f;
    document.head.appendChild(l);
  });
  var s = document.createElement('script');
  s.src = base + 'latch-viewer/latch-viewer.js';
  document.head.appendChild(s);
})();
::end
```

## How it is built
- Frame: the film frame of `latch_animation.py` (the rank016 complex). Chains: F binder (112 aa), L latch (120 aa, last
  RFdiffusion trajectory frame = closed design), P target-binding domain (PD-L1 binder, residues 27–177), T target
  (PD-L1, 18–228), A–E GvpA pentamer, G one GvpA monomer for the wall. These are the chains the script works on; the
  file has the polypeptide as **unbroken chains** instead of F, L and P (a cartoon stops at the end of its chain, so parts
  that are chains of their own show a seam where they join): R latch + linker + domain from the hinge (binder's last
  residue first), S binder's last residue + latch, U latch's last residue + linker + domain, V binder + latch, W binder +
  latch + linker + domain. The 11-residue Gly linker (K while building) is inside them. The B-factor column says which part
  a residue belongs to (0 binder, 1 binding site, 2 clash, 3 latch, 4 linker + domain); the viewer paints by it.
- **One cartoon (`latch-viewer.js`).** Where the parts sit together (steps 2 and 3 at rest) the viewer draws V or W, one
  chain from the binder through latch and linker into the domain. Where one part moves against another (the latch swinging
  in, the block swinging open, the clash view's dimmed latch) it draws pieces (F, S, U, R), and the unbroken chain fades in
  over them when the move is done. Across the hinge, once the block is open, the cartoon is really in two pieces.
- Every part is its own mmCIF assembly; every move is one matrix per part (no mesh sequences, no trajectory).
- Latch, domain and target are **one rigid block**. It swings by `block_hinge` of `film_geometry.json` (120° about the
  binder's C terminus; interface residues covered 12 → 0 of 19). Order in step 4: latch + domain swing open first, PD-L1
  is not there yet; then PD-L1 flies in along a straight line to the open domain.
- **Where the domain sits (`place_domain.py`).** The domain is joined to the latch by a drawn linker, so its pose relative to
  the latch is free. The search moves domain + PD-L1 (one rigid body, their relative pose is the film's) until: with the
  latch closed PD-L1 clashes with the binder (≥ 40 PD-L1 atoms closer than 2.2 Å to the binder, at most 450 closer than
  3 Å); domain, linker and PD-L1 stay ≥ 3 Å from binder and latch; latch + domain + linker swing without touching the binder
  (≥ 2.5 Å); open, PD-L1 and domain stay ≥ 3.5 / 3 Å from binder and GvpA pentamer; PD-L1 has a straight way in with ≥ 2.5 Å
  to binder, pentamer, latch and domain (the last 12 Å are the docking itself and are not checked against the domain).
  The search asks for 0.7 Å more than each limit. The pose found: **13 PD-L1 residues × 17 binder residues, 206 atom pairs
  closer than 2.8 Å** (closed), linker ends 25 Å apart (drawn as below). 1 of 640 restarts met everything, so the room is small.
- Clash colouring: the prep script stores 2 in the B-factor column of binder residues that clash and 1 on PD-L1 residues
  that clash (heavy atoms < 2.8 Å); the viewer's clash view paints them (`--lv-clash` in `latch-viewer.css`).
- Wall: the 1730-copy shell of the Huber 2023 model (Zenodo 10.5281/zenodo.6458345). The pentamer is fitted onto five
  consecutive wall subunits, those five are left out of the wall and replaced by the pentamer. The Drylab film takes the
  fit with 1.08 Å rmsd, which lands 245 Å from the vesicle's equator where the wall curves into the cone. The viewer
  accepts 1.7 Å (1.6 Å found) and gets a flat piece of the cylinder. The 1725 remaining copies are written as 54
  assemblies of 32 neighbouring subunits (`wall-00` …, listed in `latch-viewer.json` `wall`): Mol* draws blocks of a few
  dozen copies much faster than one object of 1725, and skips blocks that are out of the picture. Same copies, same picture.
- The construct stands off the wall by 65 Å (`away`) until step 5, as in the film.
- Binding-site residues are marked by B-factor = 1 in chain F; the viewer colours by it.

## Honest limits (what the figure caption has to say)
- Latch: RFdiffusion backbone, poly-Gly, **no sequence** (ProteinMPNN would be needed for a real construct).
- Target-binding domain: in the Drylab film it hangs off the latch end by a 3.8 Å CA–CA bridge, and its first 26
  residues are a helix that runs 30 Å away from the join, so the folded domain sits far from the latch. The viewer
  leaves those 26 residues out and builds an **11-residue Gly linker loop** (`linker()` in `prep_viewer.py`, `TRIM`)
  from the latch end to the first residue of the domain (25.3 Å apart). It is built the way loop-modelling tools build
  loops: residues are added to the latch's last residue with ideal bond lengths and angles and planar peptide bonds,
  their φ/ψ drawn from where loop glycines sit (β, polyproline II, both helical regions); cyclic coordinate descent and
  a least-squares polish then close the chain onto the domain's first residue (0.015 Å). Of 90 seeded starts (9–11
  residues) the builder keeps loops with every residue in an allowed region of the Ramachandran plot, no contact with
  themselves, and room to everything they move past (the binder at 13 angles of the swing, the pentamer once docked,
  PD-L1 on its domain, latch and domain); of those the one with the most room and the most typical angles: **≥ 4.9 Å
  from every atom**, φ/ψ in `latch-viewer.json` `geometry.linker`. The build takes ~4 minutes because of this search.
  It is **one plausible conformation, not a prediction**: nothing is known about how the real linker folds (the real
  construct has a flexible one; the AlphaFold pilot used GGSGGS). No energy, no linker design.
  Moving the domain itself next to the latch end was tried and does not fit (no clash-free position within 5 Å of the latch end).
- The swing is a **placed** hinge motion, checked only for "site freed, nothing penetrates, also on the way". Not a
  prediction.
- The clash is **placed** too: the domain pose was searched so that PD-L1 clashes with the binder when the latch is
  shut. It shows what the design intends (PD-L1 can only bind the open state), not a measured or predicted binding. The binder
  model has no side chains, so overlaps are underestimated. That the latch really opens before PD-L1 binds is a
  conformational-selection picture (the latch breathes, PD-L1 holds it open), not something we have measured. Whether PD-L1 binds past the latch is not known (design intent); a real test would be AlphaFold or docking
  scores of the open and closed state.
- Target is PD-L1 as a stand-in: the wet-lab target in the Drylab notes is GFP.

## Still to do
- `NOTICE.md`: chains P and T are probably AlphaFold Server output (file names `fold_pdl1_*`). The CIF header already
  carries the AlphaFold Server notice; **confirm the source**, then add the matching entry to `NOTICE.md` (and drop the
  header notice if they are not AlphaFold). Not edited here: `NOTICE.md` has uncommitted changes from the vesicle viewer.
- Optional: the target-binding domain could grow chain-wise from the latch (split chain P into chunks, one assembly each).
