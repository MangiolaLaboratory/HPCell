# HPCell NEWS

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