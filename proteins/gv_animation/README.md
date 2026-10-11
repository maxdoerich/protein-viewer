# Gas vesicle animation (Blender 5.1 + Molecular Nodes)

A 60-second animation that tells the story in six shots:

1. The whole gas vesicle.
2. Zoom into the wall.
3. All subunits dissolve except one. This single GvpA is shown first as a cartoon, then as a surface coloured by hydropathy.
4. Neighbouring copies assemble into ribs around it, and the EMD-14238 density fades in around the fitted model.
5. Zoom out while the whole shell grows. The two helical halves light up and meet at the seam.
6. GvpC binds and forms the cage.

## What you need
- **Blender 5.1**
- **Molecular Nodes**: Edit → Preferences → Get Extensions → search "Molecular Nodes" → Install. On 5.1 this installs version 4.5.13.

## Run it
1. Unzip. Keep `build_gv_animation.py` and the `data` folder together.
2. Open Blender and save the empty file into the unzipped folder, e.g. `gas_vesicle.blend`. The render will then be written to `render/` next to it.
3. Switch to the **Scripting** tab. In the text editor, click **Open**, choose `build_gv_animation.py`, then click **Run Script** (▶).
   The build takes about a minute. It clears the current scene first.
4. Press **Space** in the viewport to play it back, or scrub the timeline.
   The viewport preview is lighter than the final render. That's normal.
5. **Render → Render Animation** (Ctrl+F12). The output is `render/gas_vesicle.mp4`.

Rendering time with EEVEE on a recent GPU is roughly 1–4 s per frame at 1080p (1800 frames), so plan 30–120 min.
For a quick test first, set `PREVIEW = True` at the top of the script, run it again and render.

## Settings (top of the script)
| Setting | What it does |
|---|---|
| `RESOLUTION`, `FPS` | Output format |
| `ENGINE` | `"EEVEE"` is fast. `"CYCLES"` gives nicer light but is much slower (GPU recommended). |
| `OUTPUT_FORMAT` | `"MP4"` writes one file. `"PNG"` writes an image sequence, which is safer for long renders because you can resume them. |
| `COLORS` | All colours as hex codes |
| `DENSITY_THRESHOLD` | Contour level of the EMD-14238 map |
| `T` | The timeline in frames. Move numbers here to re-time shots. |
| `KEYS` (camera section) | Camera poses: frame, target, azimuth, elevation, distance |

After changing anything, click **Run Script** again. It rebuilds everything.

**Rendering in parts:** set the Start/End frame in the Output properties, use `OUTPUT_FORMAT = "PNG"`, and render.
Blender skips frames that already exist if you untick *Overwrite*.

## How the scene is built (scientific notes)
- **GvpA**: PDB 7R1C (B. megaterium GvpA2, residues 2–66), centred on its Cα centroid.
  Hydropathy is the Kyte–Doolittle scale per residue, stored in the B-factor column of `gvpa_monomer.cif`.
- **Whole shell (1730 subunits)**: the pseudo-atomic model of a complete B. megaterium gas vesicle from Zenodo (10.5281/zenodo.6458345).
  The 7R1C monomer was rigid-body superposed onto every subunit.
  The fit is Cα RMSD median 1.2 Å, max 5.6 Å, mostly at the seam and the tips. That's invisible at this scale.
  The two 852-subunit filaments (one per half) were traced with the helical operator of 7R1C (twist −3.874°, rise 0.525 Å). They drive the half-shell highlight.
- **Local patch + density (shot 4)**: 45 subunits generated with the exact helical lattice of 7R1C, so they sit precisely in EMD-14238.
  The map was masked to 6 Å around this model.
  During the zoom-out the patch morphs onto the full-shell positions (≤5 Å shift) before the rest of the shell grows in.
- **GvpC**: PDB 8GBS (Anabaena flos-aquae). Its four GvpA chains were superposed onto four consecutive GvpA of the B. megaterium lattice (Cα RMSD 2–3.8 Å), which carries the 33-residue GvpC helix into place along the rib.
  GvpC is shown as molecules of 5 repeats with a small gap between molecules. This coverage is illustrative. There is no GvpC on the sharp tips.
  The deposited GvpC is a poly-UNK backbone, so it was converted to poly-Ala.

`tools/prep_data.py` regenerates everything in `data/` from the original files.
It is not needed for rendering. It needs gemmi, mrcfile, scipy and numpy.

## Troubleshooting
- **"Molecular Nodes is not installed/enabled"**: install it via Get Extensions (see above).
- **"Could not find the 'data' folder"**: open the script from disk (Text → Open), not by pasting it. Also keep `data/` next to it.
- **The density is too faint or too bulky**: change `DENSITY_THRESHOLD` (0.6–1.2) or the density alpha keys in the script.
