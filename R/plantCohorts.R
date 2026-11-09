#' Plant new cohorts after a disturbance
#'
#' Adds one new, young cohort to `cohortData` for each unique combination of species,
#' ecoregion, and `pixelGroup` in `newPixelCohortData` — typically the cohorts that a
#' disturbance has just removed. New cohorts are planted at age 2 years with a fixed starting
#' biomass, and with zero mortality and zero annual net primary productivity.
#'
#' This is a variant of `LandR::plantNewCohorts()` that omits the `Provenance` column, for
#' workflows whose cohorts do not carry provenance information. Behaviour is otherwise the
#' same. If your cohorts do track provenance, call `LandR::plantNewCohorts()` instead.
#'
#' @param newPixelCohortData A `data.table` of the cohorts to be replaced, with columns
#'   `speciesCode`, `ecoregionGroup`, `age`, `B`, and either `pixelGroup` or `pixelIndex`.
#'   Where only `pixelIndex` is present, `pixelGroup` is looked up from `pixelGroupMap`.
#'   Rows that duplicate an existing combination of `speciesCode`, `ecoregionGroup`,
#'   `pixelGroup`, and `age` are dropped, so each species is planted at most once per group.
#' @param cohortData A `data.table` with columns `pixelGroup`, `ecoregionGroup`,
#'   `speciesCode`, `age` (years), and `B` (biomass, in g m-2), to which the new cohorts are
#'   appended.
#' @param pixelGroupMap A `SpatRaster` whose pixel values are `pixelGroup` numbers
#'   corresponding exactly to the `pixelGroup` column of `cohortData`. Used only to resolve
#'   `pixelIndex` to `pixelGroup`.
#' @param initialB The starting biomass of each new cohort, in g m-2. Defaults to 10. `NA` or
#'   `NULL` falls back to that same default.
#' @param currentTime The current simulation time, e.g. `time(sim)`.
#' @param successionTimestep The time between successive seed dispersal events; in LANDIS-II
#'   this is called the "Succession Timestep". Currently unused — accepted for compatibility
#'   with `LandR::plantNewCohorts()`.
#' @param trackPlanting Logical. If `TRUE`, the new cohorts are marked with a `planted`
#'   column set to `TRUE`, so that planted cohorts can be told apart from those that
#'   established some other way. Defaults to `FALSE`.
#'
#' @return The input `cohortData` with the new cohorts appended, as a `data.table`. New rows
#'   carry `age` 2, `B` equal to `initialB`, and `mortality` and `aNPPAct` of zero.
#'
#' @seealso [processPartialDist()], and `LandR::plantNewCohorts()` in the \pkg{LandR}
#'   package.
#'
#' @export
#' @importFrom data.table rbindlist set :=
plantCohorts <- function (newPixelCohortData, cohortData, pixelGroupMap, initialB = 10,
                          currentTime, successionTimestep, trackPlanting = FALSE)
{
  namesNCD <- names(newPixelCohortData)
  if (!isTRUE("pixelGroup" %in% namesNCD)) {
    if (isTRUE("pixelIndex" %in% namesNCD)) {
      newPixelCohortData[, `:=`(pixelGroup, as.vector(pixelGroupMap[])[pixelIndex])]
    }
    else {
      stop("newPixelCohortData must have either pixelIndex or pixelGroup")
    }
  }
  if (!is.null(newPixelCohortData[["pixelIndex"]])) {
    set(newPixelCohortData, NULL, "pixelIndex", NULL)
  }
  duplicates <- duplicated(newPixelCohortData[, .(speciesCode,
                                                  ecoregionGroup, pixelGroup, age)])
  newCohortData <- newPixelCohortData[!duplicates]
  newCohortData[, `:=`(age, 2)]
  if (isTRUE(is.na(initialB)) || is.null(initialB)) {
    initialB <- formals(plantCohorts)$initialB
  }
  newCohortData[, `:=`(B, initialB)]
  newCohortData <- newCohortData[, .(pixelGroup, ecoregionGroup,
                                     speciesCode, age, B, mortality = 0L, aNPPAct = 0L)]
  if (isTRUE(getOption("LandR.assertions"))) {
    if (isTRUE(NROW(unique(newCohortData, by = c("pixelGroup",
                                                 "age", "speciesCode"))) != NROW(newCohortData))) {
      stop("Duplicated new cohorts in a pixelGroup. Please debug plantCohorts().")
    }
  }
  if (trackPlanting) {
    newCohortData[, `:=`(planted, TRUE)]
  }
  cohortData <- rbindlist(list(cohortData, newCohortData),
                          fill = TRUE, use.names = TRUE)
  return(cohortData)
}
