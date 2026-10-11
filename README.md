# Protein viewer

Interactive Mol* (5.13) viewers for protein structures, static or animated, extracted from the iGEM Heidelberg 2026
wiki so they can be developed on their own. Nothing here depends on the wiki's generator, templates or content.

## Run it

```sh
./serve.sh            # http://localhost:8000  (any static server works; file:// does not, the viewers fetch data)
```

`index.html?v=<key>` shows one viewer; the buttons on top switch between them:

| key | viewer | what it shows |
|---|---|---|
| `gv` | `static/gv-viewer/` | four steps: one GvpA → a patch in its EMD-14238 density → the 1730-copy vesicle → + GvpC |
| `hydro4` | same, `data-preview="hydro4"` | four GvpA coloured by hydropathy, inside / outside of the wall |
| `latch` | `static/latch-viewer/` | five steps of one biosensor construct: binder → + latch → + domain → + target (latch opens) → on the wall |
| `anim5` | `static/latch-anim5/` | the same five steps for five latch designs, picked by buttons |
| `top10` | `static/latch-top10/` | ten latch designs, each closed / open / both |
| `top10bp` | same script, data from `static/latch-top10-bp/` | the ten blocking-peptide designs |

## Layout

```
index.html          demo page: wires one viewer per ?v=
serve.sh
static/             what a browser loads; every viewer finds Mol* at ../molstar/ and its data in its own folder
  tokens.css        28 design tokens, extracted from the wiki's tokens.css (values unchanged)
  molstar/          vendored Mol* 5.13 (MIT)
  gv-viewer/        gv-viewer.js/.css  (the base: layout, step buttons, card, legend, tools)
                    gv-assembly.cif, gv-shell.json, emd-14238.mrc, gv-viewer.json  (generated)
  latch-viewer/     latch-viewer.js/.css  + latch-assembly.cif, latch-viewer.json  (generated)
  latch-anim5/      latch-anim.js/.css    + latch-assembly.cif, latch-anim.json    (generated)
  latch-top10/      latch-top10.js/.css   + dNN_closed/open.pdb, designs.json      (generated)
  latch-top10-bp/   data only, shown by latch-top10.js
  dev/              dev-mode.js/.css: the developer mode of index.html (below)
proteins/           how the data is made (one README per folder, with the honest limits of each figure)
  gv_animation/     prep_data.py (raw PDB/EMDB/Zenodo → data/), prep_viewer.py (→ static/gv-viewer/),
                    build_gv_animation.py (the 60 s Blender 5.1 + Molecular Nodes film), captions.txt
  latch_animation/  prep_viewer.py (→ static/latch-viewer/), place_domain.py, data/
  latch_anim5/      prep_anim5.py (→ static/latch-anim5/)
  latch_top10/      prep_top10.py (→ static/latch-top10/ or --out latch-top10-bp)
```

The other CSS files layer on `gv-viewer.css`, so it always loads first, then `latch-viewer.css`, then the viewer's own.
CSS is token-only (`var(--…)`); the scene colours are `--gv-*` / `--lv-*` properties at the top of those two files,
read back by the scripts at start-up, so a re-skin never touches JS.

## Developer mode

The **Developer mode** button on top of `index.html` edits the viewer that is shown, in place. Everything edited is one
JSON object per viewer in `localStorage` (`pv-config:<key>`), applied on every load: switch the mode off and the
page shows the edited viewer and nothing else. **Export** / **Import** move that object to and from a file.

- **Text box, legend**: drag the tab to move, the corner to resize; click the tab and write the text of the current
  step in the panel. They lie behind the scene, so the protein is never covered.
- **Notes**: select a residue, *Add note to selection*. The box sits at a point in 3D and a line runs to its residue;
  it shows on the step it was made on, or on every step.
- **Markup** of all three: `# line` label with the blue bar, `## line` orange title, `**bold**`, `*italic*`,
  `[pine] label` a legend key (cayenne, tangerine, peach, frost, pine, cyan, grey, ink, white).
- **Buttons**: tick which of the viewer's buttons stay.
- **Selection**: click residues in the scene (Shift for a range) or drag over the **sequence view**, then pick one of
  the site's nine colours, a scheme (secondary structure, pLDDT, hydropathy) or a display mode (cartoon, surface,
  molecular surface, spheres). With nothing selected it applies to everything shown. A drag over the sequence means
  every copy of that protein; a click in the scene means that one copy.
- **Sequence view**: one row per protein on screen, with a name field.

pLDDT is read from the B-factor column. **None of the bundled files has it there**: the prep scripts use that column
for hydropathy (`gv-viewer`) or part / site flags (the latch viewers), so the scheme is only right for a model whose
B-factors are pLDDT.

A viewer takes part by setting `host.__pv = { plugin, lib, built }` and firing `pv:ready` (see the end of each
`boot()`); `index.html` loads `dev/dev-mode.css` and `dev/dev-mode.js` after it and marks the mount with `data-pv-key`.

## Put a viewer on another page

Same wiring as `index.html`: `tokens.css`, the CSS files, the script, and a mount element in the body:

```html
<link rel="stylesheet" href="static/tokens.css">
<link rel="stylesheet" href="static/gv-viewer/gv-viewer.css">
<div data-gv-viewer></div>
<script src="static/gv-viewer/gv-viewer.js"></script>
```

Mount attributes: `data-latch-viewer`, `data-latch-anim`, `data-latch-top10`; `data-base="…/"` loads the data from
another folder (a CDN), `data-binder="…"` renames the binder in the top-10 legend. Mol* is only fetched when the figure
nears the viewport.

## Rebuild the data

Needs `numpy scipy gemmi` (and `mrcfile` for `prep_data.py`). Run from anywhere; paths resolve from the script.

```sh
python proteins/gv_animation/tools/prep_viewer.py       # → static/gv-viewer/        (from proteins/gv_animation/data/)
python proteins/latch_animation/tools/prep_viewer.py    # → static/latch-viewer/
python proteins/latch_anim5/tools/prep_anim5.py  <final-folder>   # → static/latch-anim5/
python proteins/latch_top10/tools/prep_top10.py  <final-folder>   # → static/latch-top10/
```

`<final-folder>` is the Drylab PC-search output (`…/Biosensing-Animation/final`), **not in this repo**.
`prep_data.py` needs the original PDB 7R1C / 8GBS, EMD-14238 and the Zenodo model, also not in this repo. The built
data in `static/` and the intermediates in `proteins/*/data/` are enough to run everything else.

## Things to know

- **Copied from the working tree** of branch `feat/gv-viewer`, uncommitted changes included. Not the bundles
  (`gv-viewer.bundle`, `biosensing-animations.bundle`) and none of the wiki.
- **Font.** `tokens.css` loads Nunito Sans from the iGEM static server (as the wiki does); offline it falls back to
  `sans-serif`. Self-host the `.ttf` here if the tool leaves iGEM.
- **Provenance and licences**: see `NOTICE.md`. Two items are unresolved (AlphaFold Server terms for GvpC and the
  chains marked in `proteins/latch_*/README.md`; Zenodo licence). Settle them before publishing the data elsewhere.
- The latch figures are *placed* geometries, not predictions; each `proteins/*/README.md` lists what the figure can
  and cannot claim.
