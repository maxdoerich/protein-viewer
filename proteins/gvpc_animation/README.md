# GvpC repeat viewer: interactive Mol* figure

One structure (four GvpA, chains A-D, and one 33-residue GvpC repeat, chain F). The camera sits on GvpC, GvpA is see-through, and the GvpC
residues are coloured by their job: anchor residues dark blue (pine), the other interface residues light blue (frost), surface residues orange.

Residue numbers count inside the repeat (1-33); chain F is numbered 92-124 in the file (`OFFSET = 91`).

| job | repeat residues | file numbers |
|---|---|---|
| anchor | 4, 11, 12, 15, 19, 22 | 95, 102, 103, 106, 110, 113 |
| other interface | 1, 5, 16, 33 | 92, 96, 107, 124 |

## Files
| what | where |
|---|---|
| viewer script, colours | `static/gvpc-viewer/gvpc-viewer.js`, `gvpc-viewer.css` (layout comes from `gv-viewer.css`, load it first) |
| built data | `static/gvpc-viewer/gvpc-assembly.cif`, `gvpc-viewer.json` (generated, don't edit) |
| input | `data/gvpc_repeat_4gvpa.cif` (`SerGvpCR3_4GvpA_trimmed_repacked`) |
| builder | `tools/prep_viewer.py` (needs numpy) |

## Change it
- **Which residues are anchor / interface:** `ANCHOR` / `INTERFACE` at the top of `tools/prep_viewer.py`, then
  `python tools/prep_viewer.py`. The viewer paints by the B-factor column of chain F (2 anchor, 1 interface, 0 rest).
- **Colours:** `--gc-*` at the top of `gvpc-viewer.css`. **How see-through GvpA gets:**
  `GVPA_SEE_THROUGH` in `gvpc-viewer.js`. **Camera:** the `frames` in the prep script (the `gvpc` frame).
- **Colour GvpA per chain:** colour the four chains separately in dev mode
  (Selection → colour by chain), or build one assembly per chain in `prep_viewer.py`.

## Put it on a page
Same wiring as the other viewers: `tokens.css`, `gv-viewer/gv-viewer.css`, `gvpc-viewer/gvpc-viewer.css`, the script, and
`<div data-gvpc-viewer></div>`.

## Honest limits
- The GvpC repeat is a **model** (AlphaFold Server prediction, repacked), not an experimental structure; the CIF header
  carries the AlphaFold Server notice. See `NOTICE.md`.
- Anchor / interface residues are as given by the team; this viewer does not compute contacts.
