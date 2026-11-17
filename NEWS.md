# hallfast 0.3.0

* New `hall_state()` and `hall_step()`: a carried state (one row per
  individual) advanced one interval at a time, for microsimulations that
  simulate their population year by year and cannot pass the whole history of
  intake changes at once. `Y` calls of `hall_step()` give exactly (bit for bit)
  the weights of one `hall_bw_annual()` call over the `Y` intervals; rows are
  independent, so a state can be subset or bound together as individuals enter
  and leave. States can start from body weight and height or, as in the BMI
  interfaces, from BMI with height optional (fractional changes only).
* The C++ header moved from `src/hall_model.h` to
  `inst/include/hallfast/hall_model.h`, so packages can use the model directly
  with `LinkingTo: hallfast` and `#include "hallfast/hall_model.h"`.

# hallfast 0.2.0

* New `EIchange_unit` argument in `hall_bw_annual()`, `hall_bw_long()` and
  `hall_bw_reference()`. With `"fraction"`, `EIchange` is a fraction of each
  person's baseline energy intake (for example `-0.05` for a 5% cut). The
  default, `"kcal"`, keeps the previous behaviour.
* New `hall_bmi_annual()` and `hall_bmi_long()`: take baseline BMI and return
  BMI (`BMI_0 * BW_t / BW_0`). With fractional changes, the default, height is
  optional; missing heights are replaced by `hall_ref_height` (Health Survey
  for England 2019 adult means). At a given BMI the proportional change in
  BMI then depends very little on height, whereas with a fixed kcal/day
  change it scales roughly with 1/height^2, so height is required in that
  case.
* New vignette section "BMI and proportional changes in intake".

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
