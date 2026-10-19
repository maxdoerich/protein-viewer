# GvpC charge viewer (`?v=charge`)

Four GvpA and one GvpC repeat of PDB 8GBS (*Anabaena flos-aquae*; Dutka et al., Structure 31, 518 (2023),
doi 10.1016/j.str.2023.03.011), each surface coloured by its Coulombic potential. It is the scene of the vault's ChimeraX
image (`SecondBrain/assets/scripts/chimerax_gvpa_gvpc_charge_2026-10-13.cxc`), with GvpC's charge added and one
chip, **Remove GvpC / Add GvpC** (GvpC fades off the wall and back). The text box gives GvpC's net charge and its
charged residues, read from `charge-viewer.json`.

## Rebuild

```sh
/Applications/ChimeraX-1.12.app/Contents/bin/ChimeraX --nogui --exit --script proteins/gvpc_charge/tools/prep_viewer.py
```

Reads `data/8gbs.cif` (RCSB download, 2026-10-19) and `../gv_animation/data/ana_gvpc_afdb.cif` (AF-P09413-F1 v6),
writes `static/charge-viewer/gvpc-charge.cif` and `charge-viewer.json`. About 30 s, no GUI needed.

## How it is made, and its honest limits

- **GvpA**: the four deposited chains A1–A4 (residues 2–66) as they are. Some side chains are incomplete in 8GBS
  (ChimeraX warns, e.g. A1 Val30), as in the ChimeraX image; their charges are ChimeraX's for the atoms present.
  Residues 2 and 66 are treated as charged termini, as `coulombic` does by default.
- **GvpC has no side chains in 8GBS** (33 × UNK: N, CA, C, O, CB), so the deposited structure carries no GvpC charge.
  Its residue numbers 85–117 are UniProt P09413 repeat 3 (UniProt "Repeat" features, queried 2026-10-18). The side
  chains are taken from the AlphaFold DB model AF-P09413-F1, residue by residue: each AlphaFold residue is superposed
  on its N, CA, C, CB in 8GBS and only the atoms beyond CB are copied, so the 8GBS backbone stays as deposited.
  The AlphaFold repeat as a whole sits on that backbone with a Cα RMSD of 0.99 Å.
- **Three rotamers replaced.** Where an AlphaFold side chain ran into the wall (heavy atom < 2.5 Å from GvpA: L88,
  R103, Q106), ChimeraX `swapaa … rotLib Dunbrack criteria cp` put in the library rotamer with the fewest clashes
  (ties: the most probable). AlphaFold predicted GvpC alone, so its rotamers know nothing of GvpA. None is left
  under 2.5 Å. **GvpC's charge distribution is therefore a model**, not a measurement; the sequence (which residues
  are charged, the net charge) is not affected by it, the positions of the charges are.
- **Charges** are those ChimeraX `coulombic` assigns: AMBER ff14SB, hydrogens added to a copy and folded into the
  heavy atoms, standard protonation (His neutral). GvpC is a piece cut out of a 193-residue protein, so for the
  assignment it is extended by AlphaFold residues 84 and 118, which take the terminal charges and are then deleted:
  residues 85 and 117 carry internal-residue charges. Net charge of repeat 3: **−1** (K86, R103, K109 positive; E87,
  E90, E108, E113 negative).
- **Potential**: ChimeraX's own `potential_at_points` (distance-dependent dielectric ε = 4r, 1.4 Å off the surface,
  kcal/(mol·e)), the defaults of the ChimeraX image. The GvpA surface carries the potential of the four GvpA (as
  `coulombic #4/A1-A4` did), the GvpC surface that of GvpC alone: its own charge, the same whether or not GvpC is shown.
  GvpA is one surface around all four chains (the image had one per chain; the outside is the same, the faces between
  the chains are gone).
- **Per atom, not per surface point.** Mol* colours a surface by atom, so each atom's B-factor is the mean potential
  of the surface points ChimeraX assigns to it (atoms without any: their residue's mean). Colours therefore change
  from atom to atom rather than smoothly, and Mol*'s surface is not ChimeraX's triangle for triangle. Same range and
  colours as the image: −10 Cayenne Red, 0 white, +10 Pine Blue.
- The scene's light (cel shading, the family's ambient occlusion) makes pale oranges look a little yellower than in
  the ChimeraX image; the data are the same.
