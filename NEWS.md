# hallfast 0.1.0

* First release.
* `hall_bw_annual()` (matrix interface) and `hall_bw_long()` (long-format
  interface): annual-resolution solver for the Hall et al. (2011) adult body
  weight model. The web appendix is the reference specification.
* Exact closed forms for adaptive thermogenesis, extracellular fluid and
  glycogen; one scalar ODE for lean mass, solved with adaptive
  Dormand-Prince 5(4). About 140-320 times faster than daily RK4 (package
  'bw'), with maximum error below 1e-4 kg in a realistic population.
* Optional proportional appetite feedback (Hall 2024) through `appetite_k`,
  with published fitted values in `hall_appetite_presets`.
* Default physical activity level of 1.6 (SACN 2011 median 1.63; Hall 2024
  used 1.65), instead of the appendix's 1.5.
* `model = "bw"` reproduces the 'bw' package v1.0.0, which differs from the
  appendix in its glycogen energy term and in how K is set when baseline
  energy intake is supplied.
* `hall_bw_reference()`: slow, independent high-precision reference solver
  (full four-compartment system, fixed-step RK4) for validation.
* Optional OpenMP multi-threading (`threads`, option `hallfast.threads`).
