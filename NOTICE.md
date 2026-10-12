# Third-party attribution (protein viewer)

Carried over from the wiki's `NOTICE.md`; only what this folder ships.

| asset | version | source | license |
|---|---|---|---|
| `static/molstar/` | 5.13.0 | [Mol*](https://molstar.org/) (`build/viewer` of the npm package, source maps removed) — © 2017–2026 Mol* contributors | [MIT](https://github.com/molstar/molstar/blob/master/LICENSE) (copy in that folder) |

## Structural data

`static/gv-viewer/` is derived by `proteins/gv_animation/tools/prep_viewer.py` from PDB 7R1C, PDB 8GBS and EMD-14238
(PDB and EMDB data are [CC0](https://www.wwpdb.org/about/usage-policies)), from the full-vesicle pseudo-atomic model of
Huber et al. (Zenodo 10.5281/zenodo.6458345), and from a GvpC structure predicted with AlphaFold Server (job
`Ser_GvpC_single`, kept as `proteins/gv_animation/data/gvpc_af.cif`). The latch viewer (`latch-viewer/`) reuses that wall and GvpC, plus the team's RFdiffusion designs from the Drylab PC search.

**AlphaFold Server output.** The GvpC chain of `gv-assembly.cif` is a Derivative of AlphaFold Server output, provided
under the [AlphaFold Server Output Terms of Use](https://alphafoldserver.com/output-terms): non-commercial use only, with a
conspicuous notice (the first lines of that file; repeat it in any figure caption).

## Still to confirm

- Zenodo 10.5281/zenodo.6458345: check the licence on the record and add it above.
- The latch chains P and T (PD-L1 domain / target) are probably AlphaFold Server output too (file names `fold_pdl1_*`);
  confirm the source and extend the notice above.
- The wiki is CC BY 4.0, which allows commercial use; the AlphaFold terms do not. Decide, with an advisor, what the
  notice on any published copy of this data must say.
