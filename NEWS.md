# HPCell NEWS

## HPCell 0.6.0

### New features

* Added a `species` argument to `initialise_hpc()` (`"human"` by default) and
  threaded it through empty-droplet filtering, live-cell QC, cell-type
  annotation, and CellChat.
* Mitochondrial genes are now resolved from the species Ensembl database
  (`EnsDb.Hsapiens.v86` or `EnsDb.Mmusculus.v79`) by chromosome `MT`, for both
  gene symbols and Ensembl IDs. Ribosome genes are matched with `^(RPL|RPS)`
  for human and `^(Rpl|Rps)` for mouse, mapping Ensembl IDs to symbols first.
* Human cells are annotated with Blueprint, Monaco, and optional Azimuth.
  Mouse cells are annotated with `celldex::MouseRNAseqData` (fine and coarse
  labels). CellChat uses `CellChatDB` and `PPI` for the selected species.

## HPCell 0.5.3

### Bug fixes

* Fixed `transform_utility()` to no longer assume the first assay of the input
  object is named `"X"`. The canonical output assay name `"X"` is now set
  explicitly via a single `assay_name` variable rather than being inferred from
  the input, so inputs whose first assay carries any other name are handled
  correctly.

### Warnings

* `transform_utility()` now emits a warning when the input SCE's assay is not
  named `"X"`, reporting the original name before renaming it to `"X"` for
  downstream consistency.