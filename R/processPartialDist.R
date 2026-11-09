#' Apply a partial disturbance to a landscape
#'
#' Removes selected tree cohorts from the pixels affected by a disturbance in a given year,
#' and optionally replants them. A *partial* disturbance is one that need not affect every
#' pixel of a `pixelGroup`: where only some of a group's pixels are disturbed, the group is
#' split, so that the disturbed pixels get a new `pixelGroup` and the undisturbed pixels keep
#' the original one. Which cohorts are removed is controlled by species and by an age window,
#' so a disturbance can be made to take, say, only mature individuals of one species.
#'
#' The disturbance is applied for a single year: `partialDistTable` is first subset to the
#' rows whose `distYear` equals `currentTime`.
#'
#' @param cohortData A `data.table` with columns `pixelGroup`, `ecoregionGroup`,
#'   `speciesCode`, `age` (years), and `B` (biomass, in g m-2). Each row is one cohort of one
#'   species in one `pixelGroup`.
#' @param partialDistTable A `data.table` describing the disturbance, with columns
#'   `speciesCode` (which species the disturbance removes) and `distYear` (the simulation
#'   year it occurs in). It may also carry `ageMin` and `ageMax`, the inclusive age window,
#'   in years, within which a cohort of that species is eligible for removal; where those
#'   columns are absent they are filled with the `ageMin` and `ageMax` arguments.
#' @param pixelGroupMap A `SpatRaster` whose pixel values are `pixelGroup` numbers
#'   corresponding exactly to the `pixelGroup` column of `cohortData`.
#' @param replanting Logical. If `TRUE` (the default), one new cohort is planted per species
#'   per disturbed `pixelGroup`, via [plantCohorts()], even where several cohorts of that
#'   species were removed. If `FALSE`, removed cohorts are not replaced.
#' @param currentTime The current simulation time, e.g. `time(sim)`. Used to select the rows
#'   of `partialDistTable` to apply.
#' @param partialDistLoc Where on the landscape the disturbance happens. One of: `NULL` (the
#'   default), meaning the whole landscape is eligible; an `sf` polygon layer, which is
#'   rasterized onto `pixelGroupMap`; or a `SpatRaster` in which disturbed pixels have the
#'   value `1`. As a convenience, a raster that is entirely `NA` is treated as its own
#'   bounding box, so that every pixel of `pixelGroupMap` falling inside that extent is
#'   disturbed.
#' @param ageMin Lowest cohort age, in years, eligible for removal; used only for rows of
#'   `partialDistTable` that do not supply their own `ageMin`. Defaults to 50.
#' @param ageMax Highest cohort age, in years, eligible for removal; used only for rows of
#'   `partialDistTable` that do not supply their own `ageMax`. Defaults to 500.
#'
#' @return A `list` of length two:
#'   \describe{
#'     \item{`cohortData`}{the updated `data.table`, with disturbed cohorts removed, any
#'       replanted cohorts added, `totalBiomass` recalculated per `pixelGroup`, and
#'       `pixelGroup` stored as an integer.}
#'     \item{`pixelGroupMap`}{the updated `SpatRaster`, in which partially disturbed groups
#'       have been split and the disturbed pixels carry their new `pixelGroup` values.}
#'   }
#'
#' @seealso [plantCohorts()]
#'
#' @export
#' @importFrom data.table copy data.table rbindlist :=
#' @importFrom terra crds ext ncell rasterize setValues values vect
processPartialDist <- function(cohortData, partialDistTable, pixelGroupMap, replanting = TRUE, currentTime, partialDistLoc = NULL, ageMin = 50L, ageMax = 500L) {

  ## Subset partialDistTable to disturbances of current year and add age columns if not provided
  partialDistTable <- partialDistTable[distYear == currentTime]
  if (!"ageMin" %in% names(partialDistTable)){
    partialDistTable$ageMin <- ageMin
  }
  if (!"ageMax" %in% names(partialDistTable)){
    partialDistTable$ageMax <- ageMax
  }

  ## ID pixels that can be disturbed
  if (!is.null(partialDistLoc)) {
    if (inherits(partialDistLoc, "sf")) {
      partialDistLoc <- terra::vect(partialDistLoc)
      partialDistLoc <- terra::rasterize(
        partialDistLoc,
        pixelGroupMap,
        field = 1,     # all disturbed area gets value 1
        background = NA)}

    distValues <- terra::values(partialDistLoc, mat = FALSE)

    if (all(is.na(distValues))) { #if empty raster on which partialDistTable is to be applied
      ## no pixel is flagged, so fall back to the raster's extent: every pixel of
      ## pixelGroupMap inside that bounding box is treated as disturbed
      coords <- terra::crds(pixelGroupMap, df = TRUE)
      ext <- terra::ext(partialDistLoc)

      distArea <- coords$x >= ext$xmin &
        coords$x <= ext$xmax &
        coords$y >= ext$ymin &
        coords$y <= ext$ymax

      distPixels <- which(distArea)

    } else { #if raster has values 1=disturbed 0=undisturbed
      distPixels <- which(distValues == 1)
    }

  } else {#partialDistLoc is NUlL and partialDistTable is applied over entire landscape
    distPixels <- seq_len(terra::ncell(pixelGroupMap))
  }

  ## determine if all pixels in a pixelGroup is disturbed or now
  pixelGroups <- terra::values(pixelGroupMap, mat = FALSE)
  pixelGroupsClean <- pixelGroups[!is.na(pixelGroups) & pixelGroups > 0]
  distPixelGroupsClean <- pixelGroups[distPixels]
  distPixelGroupsClean <- distPixelGroupsClean[!is.na(distPixelGroupsClean) & distPixelGroupsClean > 0]
  pixelGroupIDs <- sort(unique(pixelGroupsClean))
  totalTab <- tabulate(pixelGroupsClean, nbins = max(pixelGroupIDs))
  distTab <- tabulate(distPixelGroupsClean, nbins = max(pixelGroupIDs))
  totalPixelCount <- totalTab[pixelGroupIDs]
  distPixelCount <- distTab[pixelGroupIDs]

  pixelGroupDist <- data.table(
    pixelGroup = pixelGroupIDs,
    disturbedPixels = distPixelCount,
    totalPixels = totalPixelCount)

  pixelGroupDist[, disturbed := disturbedPixels > 0]
  pixelGroupDist[, partial := disturbedPixels > 0 & disturbedPixels < totalPixels]

  partialPG <- pixelGroupDist[partial == TRUE, pixelGroup]
  fullPG    <- pixelGroupDist[disturbed == TRUE & partial == FALSE, pixelGroup]

  ## create new pixel groups from groups that are not entirely disturbed (new group = disturbed, old group = undisturbed)
  maxPixelGroup <- max(cohortData$pixelGroup)
  # newPixelGroupmMap <- pixelGroups
  newPixelGroupmMap <- terra::values(pixelGroupMap, mat = FALSE, na.rm = FALSE)
  newCohorts <- list()
  newPGs <- integer()

  for (pg in partialPG) {
    # create new pixel group
    maxPixelGroup <- maxPixelGroup + 1
    newPixelGroup <- maxPixelGroup

    newPGs <- c(newPGs, newPixelGroup)

    pgPixels <- which(pixelGroups == pg)
    disturbedSubset <- intersect(pgPixels, distPixels)

    # update pixelGroupMap with new pixelGroup in the disturbed area
    newPixelGroupmMap[disturbedSubset] <- newPixelGroup

    ## duplicate appropriate cohortData row for new disturbed group
    newRows <- copy(cohortData[pixelGroup == pg])
    newRows[, pixelGroup := newPixelGroup]

    newCohorts[[length(newCohorts) + 1]] <- newRows
  }

  ## update cohortData to add the new pixelGroup
  if (length(newCohorts) > 0)
    cohortData <- rbindlist(list(cohortData, rbindlist(newCohorts)))

  ## Update pixelGroupMap raster
  pixelGroupMap <- terra::setValues(pixelGroupMap, newPixelGroupmMap)

  ## create table of disturbed cohorts
  # disturbedPixelGroups <- unique(newPixelGroupmMap[distPixels])
  disturbedPixelGroups <- unique(c(fullPG, newPGs))
  distCohorts <- cohortData[pixelGroup %in% disturbedPixelGroups]
  distCohorts <- partialDistTable[
    distCohorts,
    on = .(speciesCode, ageMin <= age, ageMax >= age),
    nomatch = 0]

  ## clean up the columns to match cohortData
  ## NOTE: a non-equi join renames the result columns after the join's left-hand side, so
  ## the column now called `ageMin` holds each matched cohort's actual `age`. Assigning
  ## `age = ageMin` therefore restores the real age rather than overwriting it.
  distCohorts[, `:=`(
    age = ageMin,
    ageMin = NULL,
    ageMax = NULL)]
  colsToKeep <- c("speciesCode", "ecoregionGroup", "age", "B", "pixelGroup")
  distCohorts <- distCohorts[, ..colsToKeep]

  ## remove the distCohorts from sim$cohortData
  cohortData <- cohortData[!distCohorts, on = .(speciesCode, pixelGroup, age)]

  ## replanting
  if (replanting == TRUE){
  ## one new cohort per species per pixelGroup, even if several of its cohorts were cut
  rePlanting <- plantCohorts(
    newPixelCohortData = unique(distCohorts, by = c("pixelGroup", "speciesCode")),
    cohortData = cohortData,
    pixelGroupMap = pixelGroupMap,
    initialB = 10,
    currentTime = currentTime,
    trackPlanting = FALSE)

  ## recalc total biomass w/ replanting
  cohortData <- rePlanting[
    , totalBiomass := sum(B),
    by = pixelGroup]
  } else {
    #recalc total biomass if cohorts are removed and not replanted
    cohortData <- cohortData[
      , totalBiomass := sum(B),
      by = pixelGroup]
  }
  #ensure cohortData$pixelGroup is an integer
  cohortData[, pixelGroup := as.integer(pixelGroup)]
  ## return new cohortData and pixelGroupMap
  return(list(
    cohortData = cohortData,
    pixelGroupMap = pixelGroupMap
  ))
}
