# HPCell NEWS

## HPCell 0.6.0

### Bug fixes

* Fixed Azimuth-based cell-type label transfer (`annotation_label_transfer()`) failing
  with SeuratObject >= 5.0.0. Azimuth 0.5.0 still calls the defunct
  `GetAssayData(slot = ...)` interface; HPCell now temporarily patches the
  Azimuth import environment so `slot` is mapped to `layer` during
  `RunAzimuth()`.
* Updated `GetAssayData()` calls from deprecated `slot = "counts"` to
  `layer = "counts"` in `empty_droplet_id()`, `empty_droplet_threshold()`, and
  `alive_identification()` for Seurat 5 compatibility.
* Fixed `transform_utility()` so the `identity` transform also applies
  `limit_max_to_scale()` before transformation. Previously only `expm1` was
  pre-scaled, which could leave high-count samples unscaled and cause downstream
  failures in the assay transformation pipeline.

### Improvements

* `non_batch_variation_removal()` now catches `SCTransform()` failures for edge
  cases with very few overdispersion genes, emits a warning, and returns
  `NULL` instead of stopping the pipeline.
* `initialise_hpc()` gains a new `default_controller` argument to set the
  default `crew` controller for targets that do not specify their own
  controller via `tar_resources()`.

### Documentation

* Expanded and corrected roxygen documentation for `initialise_hpc()`, including
  the new `default_controller` argument.

### Development

* Added CellNexus 2025 pipeline scripts under `dev/cellnexus-2025-scripts/`
  (steps 1–10 for census download, metadata preparation, HPCell execution,
  local-cache splitting, pseudobulk preparation, and metadata unification).
* Updated CellNexus 2024 pipeline scripts, including renamed
  `step9_unify_and_update_sce_metadata.R` and new steps for pseudobulk
  preparation and CellNexus querying.
* Added `.positai` and `.claude` to `.Rbuildignore` and `.gitignore`.
