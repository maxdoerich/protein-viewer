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
| `charge` | `static/charge-viewer/` | four GvpA + one GvpC repeat of PDB 8GBS, surfaces coloured by Coulombic potential (the vault's ChimeraX image); a chip removes / adds GvpC; no text box (the wiki caption describes it), the scale stays in front of the protein. GvpC's side chains are AlphaFold's on the 8GBS backbone (`proteins/gvpc_charge/README.md`) |
| `open` | `static/open-viewer/` | any structure: fetch a PDB ID (PDBe) or upload a `.pdb` / `.cif` / `.mmcif` / `.bcif` file; drawn as a cartoon in the family's colours, with the same text box, legend, buttons and developer mode |

## Layout

```
index.html          demo page: wires one viewer per ?v=
serve.sh
static/             what a browser loads; every viewer finds Mol* at ../molstar/ and its data in its own folder
  tokens.css        28 design tokens, extracted from the wiki's tokens.css (values unchanged)
  molstar/          vendored Mol* 5.13 (MIT)
  gv-viewer/        gv-viewer.js/.css  (the base: layout, step buttons, card, legend, tools)
                    gv-assembly.cif, gv-shell.json, emd-14238.mrc, gv-viewer.json  (generated)
                    gv-snapshot.js  (the "Save image" button, for every viewer: a 3840 px PNG of the current frame)
  latch-viewer/     latch-viewer.js/.css  + latch-assembly.cif, latch-viewer.json  (generated)
  charge-viewer/    charge-viewer.js/.css  + gvpc-charge.cif, charge-viewer.json  (generated)
  open-viewer/      open-viewer.js/.css  (fetch by PDB ID / upload a file; no data of its own; the bar above it is its only CSS)
  dev/              dev-mode.js/.css: the developer mode of index.html (below)
proteins/           how the data is made (one README per folder, with the honest limits of each figure)
  gv_animation/     prep_data.py (raw PDB/EMDB/Zenodo → data/), prep_viewer.py (→ static/gv-viewer/),
                    build_gv_animation.py (the 60 s Blender 5.1 + Molecular Nodes film), captions.txt
  latch_animation/  prep_viewer.py (→ static/latch-viewer/), place_domain.py, data/
  gvpc_charge/      prep_viewer.py (a ChimeraX script → static/charge-viewer/), data/8gbs.cif
```

The other CSS files layer on `gv-viewer.css`, so it always loads first, then `latch-viewer.css`, then the viewer's own.
CSS is token-only (`var(--…)`); the scene colours are `--gv-*` / `--lv-*` properties at the top of those two files,
read back by the scripts at start-up, so a re-skin never touches JS.

## Developer mode

The **Developer mode** button on top of `index.html` edits the viewer that is shown, in place. Everything edited is one
JSON object per viewer in `localStorage` (`pv-config:<key>`), applied on every load: switch the mode off and the
page shows the edited viewer and nothing else. **Export scene** / **Import scene** move it to and from a file
([scene files](#scene-files)).

- **Text box, legend**: drag the tab to move, the corner to resize; click the tab and write the text of the current
  step in the panel. They lie behind the scene, so the protein is never covered.
- **Notes**: select a residue, *Add note to selection*. The box sits at a point in 3D and a line runs to its residue;
  it shows on the step it was made on, or on every step.
- **Markup** of all three: `# line` label with the blue bar, `## line` orange title, `**bold**`, `*italic*`,
  `[pine] label` a legend key (cayenne, tangerine, peach, frost, pine, cyan, grey, ink, white).
- **Buttons**: tick which of the viewer's buttons stay.
- **Selection**: click residues in the scene (Shift for a range) or drag over the **sequence view**, then pick one of
  the site's nine colours, a scheme (secondary structure, pLDDT, hydropathy) or a display mode (cartoon, surface,
  molecular surface, spheres, or cartoon + surface: the cartoon with a see-through surface over it). With nothing selected it applies to everything shown. A drag over the sequence means
  every copy of that protein; a click in the scene means that one copy.
- **Sequence view**: one row per protein on screen, with a name field.
- **Structures** (open viewer): the list of what is in the scene. Per structure: *Select* (the whole of it), *Align to…*
  (another structure, or the current selection: residues are paired by sequence, then the structure is moved by the
  least-squares superposition of their CA atoms), a name field, *Delete*. *New structure from selection* copies the
  selected residues into a structure of their own, in place, so it can be aligned, painted or deleted separately.

pLDDT is read from the B-factor column. **None of the bundled files has it there**: the prep scripts use that column
for hydropathy (`gv-viewer`) or part / site flags (the latch viewers), so the scheme is only right for a model whose
B-factors are pLDDT.

A viewer takes part by setting `host.__pv = { plugin, lib, built }` (the open viewer adds `skey`, the unique name of each
structure, and `host.__pvObjects` for the list above) and firing `pv:ready` (see the end of each
`boot()`); `index.html` loads `dev/dev-mode.css` and `dev/dev-mode.js` after it and marks the mount with `data-pv-key`.

## Scene files

Developer mode's **Export scene** writes one JSON file, **Import scene** rebuilds the viewer from it, in place, on any
device (no server, no path):

```json
{ "format": "protein-viewer-scene", "version": 1, "viewer": "open",
  "structures": [                                       // open viewer: every structure of the scene, with its name and its pose (4×4)
    { "name": "1CRN", "pose": [1,0,0,0, 0,1,0,0, 0,0,1,0, 0,0,0,1], "pdb": "1CRN" }
  ],                                                    // or "file": "x.pdb", "text": "…" / "file": "x.bcif", "base64": "…" 
  "camera":    { "position": [...], "target": [...], … },
  "config":    { "box": {…}, "text": {…}, "hide": [], "paint": [], "mode": {}, "notes": [], "names": {}, "seq": false } }
```

- `config` is everything edited in developer mode: text and legend (text, position, size), notes, colours, display
  modes, hidden buttons, protein names. Every step's edits are in it; which step is on screen is not.
- `structures` and `camera` exist for the **open** viewer only. An uploaded file (or a structure made from a
  selection) travels *inside* the JSON (text, or base64 for `.bcif`), so the file is as big as the structures.
  A file with a single `"structure": {…}`, as written before there were several, still imports. The other viewers keep their structure in their folder,
  so their scene is the config alone.
- A scene imports into the viewer it was made in (`viewer`); anywhere else it says which `?v=` to open.
- A config exported before scene files existed (no `format`) still imports.

## Open a structure

`?v=open`: type a PDB ID and **Fetch**, or **Upload file** (several at once are fine). Each one is *added* to the scene
(delete one from developer mode's list of structures), and drawn like the others: polymer only (protein and nucleic acid; water and ligands are left out), cartoon, secondary-structure colours,
the same text box / legend / buttons, all editable in developer mode. Edits are kept in `pv-config:open`: the text of
a box belongs to the structure it was written for (the one added last), a colour to the sequence it was applied to
(or to the residues clicked, in that structure), a display mode to the structure it was applied to. The PDB IDs
fetched (with their names and poses) are remembered and fetched again on the next visit; an uploaded file or a
structure made from a selection is not (export a scene to keep them: the files go inside).

## Save image

The camera button among the round buttons saves the current frame as a PNG, 3840 px on its long side whatever the
size on screen (same proportions as the figure, so the same framing). It holds the protein, the background, the text
box, legend and notes; no buttons. Mol* draws the protein at that size; the boxes are rasterised from their computed
styles, with the page's web font inlined. The text box is left out on a phone, where it sits under the figure, not on it.
Add `static/gv-viewer/gv-snapshot.js` after the viewer script to get the button on another page.

## Put a viewer on another page

Same wiring as `index.html`: `tokens.css`, the CSS files, the script, and a mount element in the body:

```html
<link rel="stylesheet" href="static/tokens.css">
<link rel="stylesheet" href="static/gv-viewer/gv-viewer.css">
<div data-gv-viewer></div>
<script src="static/gv-viewer/gv-viewer.js"></script>
```

Mount attributes: `data-latch-viewer`, `data-charge-viewer`, `data-open-viewer`; `data-base="…/"` loads the data from
another folder (a CDN). Mol* is only fetched when the figure nears the viewport.

## Rebuild the data

Needs `numpy scipy gemmi` (and `mrcfile` for `prep_data.py`). Run from anywhere; paths resolve from the script.

```sh
python proteins/gv_animation/tools/prep_viewer.py       # → static/gv-viewer/        (from proteins/gv_animation/data/)
python proteins/latch_animation/tools/prep_viewer.py    # → static/latch-viewer/
/Applications/ChimeraX-1.12.app/Contents/bin/ChimeraX --nogui --exit --script proteins/gvpc_charge/tools/prep_viewer.py   # → static/charge-viewer/
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
  chains marked in `proteins/latch_animation/README.md`; Zenodo licence). Settle them before publishing the data elsewhere.
- The latch figure is a *placed* geometry, not a prediction; each `proteins/*/README.md` lists what the figure can
  and cannot claim.
