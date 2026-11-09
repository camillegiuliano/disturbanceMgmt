# disturbanceMgmt 0.0.1.9001

* new functions `processPartialDist()` and `plantCohorts()`, brought in from the
  `LandRCBM_partialDisturbance` module and documented with roxygen2. `processPartialDist()`
  removes tree cohorts from the pixels hit by a disturbance in a given year, splitting a
  `pixelGroup` where only part of it is affected; `plantCohorts()` replants them.
* `data.table`, `terra`, and `LandR` added to `Imports`; `LandR` is pinned to
  `PredictiveEcology/LandR@development` in `Remotes`.
* `plantCohorts()` now falls back to its own `initialB` default rather than looking up
  `LandR::plantNewCohorts()`, and its optional assertion is gated on
  `isTRUE(getOption("LandR.assertions"))`, so it no longer errors when that option is
  unset.

# disturbanceMgmt 0.0.1.9000

* initial package scaffold: `DESCRIPTION`, `R/`, `tests/testthat/`, and CI workflows
  calling the shared actions at `PredictiveEcology/actions@main`.
