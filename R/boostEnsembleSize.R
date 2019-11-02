# "`-''-/").___..--''"`-._
#  (`6_ 6  )   `-.  (     ).`-.__.`)   WE ARE ...
#  (_Y_.)'  ._   )  `._ `. ``-..-'    PENN STATE!
#    _ ..`--'_..-_/  /--'_.' ,'
#  (il),-''  (li),'  ((!.-'
#
# Author: Weiming Hu <weiming@psu.edu>
#         Geoinformatics and Earth Observation Laboratory (http://geolab.psu.edu)
#         Department of Geography and Institute for CyberScience
#         The Pennsylvania State University
#


#' EITrans::boostEnsembleSize
#'
#' EITrans::boostEnsembleSize boosts the number of ensemble memebers using
#' an empirical approach.
#'
#' @author Weiming Hu \email{weiming@@psu.edu}
#'
#' @details
#' An important assumption is that AnEn.LOO.test is generated using the same config
#' as AnEn except for the config$test_times_compare and config$num_members being changed.
#'
#' @md
#' @export
boostEnsembleSize <- function(
  AnEn, config, AnEn.LOO.test = NA,
  scale.size = 10, show.progress = T,
  silent = F, avoid.duplicates = T,
  analog.member = 'analogs', similarity.member = 'similarity') {

  # Sanity checks
  stopifnot(class(AnEn) == 'AnEn')
  stopifnot(config$mode == 'independentSearch')
  stopifnot(config$preserve_similarity)

  if (identical(AnEn.LOO.test, NA)) {
    config$test_times_compare <- config$search_times_compare
    config$max_num_sims <- scale.size
    config$num_members <- scale.size
    config$preserve_mapping <- F
    config$preserve_similarity <- F

    if (!silent) cat(
      'Carrying out AnEn LOO test with', config$num_members,
      'members and for', length(config$test_times_compare),
      'test times ...\n')
    AnEn.LOO.test <- generateAnalogs(config)

  } else {
    stopifnot(dim(AnEn[[analog.member]])[1] == dim(AnEn.LOO.test$analogs)[1])
    stopifnot(length(config$search_times_compare) == dim(AnEn.LOO.test$analogs)[2])
    stopifnot(dim(AnEn[[analog.member]])[3] == dim(AnEn.LOO.test$analogs)[3])
  }

  if (dim(AnEn.LOO.test$analogs)[4] < (scale.size - 1)) {
    stop("The precomputed Analogs do not have enough ensemble members to select.")
  }

  new.dims <- dim(AnEn[[analog.member]])
  new.dims[4] <- new.dims[4] * scale.size
  analogs.boost <- array(NA, dim = new.dims)

  # Keep the ensemble members from the original analogs
  analogs.boost[, , , 1:dim(AnEn[[analog.member]])[4], ] <- AnEn[[analog.member]]

  # Where should I insert the extra ensemble members
  insert.start <- dim(AnEn[[analog.member]])[4] + 1
  insert.end <- dim(AnEn[[analog.member]])[4] * scale.size
  insert.offset <- insert.end - insert.start

  if (!silent) cat('Boosting ensemble size from', dim(AnEn[[analog.member]])[4],
                   'to', new.dims[4], '...\n')
  rm(new.dims)

  if (show.progress) {
    pb <- txtProgressBar(max = prod(dim(AnEn[[analog.member]])[1:3]), style = 3)
    counter <- 0
  }

  na.exists <- F

  for (i.grid in 1:dim(AnEn[[analog.member]])[1]) {
    for (i.test in 1:dim(AnEn[[analog.member]])[2]) {
      for (i.flt in 1:dim(AnEn[[analog.member]])[3]) {

        # Get the past day index in forecast times for this partitular AnEn
        past.fcst.day <- AnEn[[similarity.member]][i.grid, i.test, i.flt, , 3]

        if (any(is.na(past.fcst.day))) {
          na.exists <- T

        } else {
          # Get the past day index in search search forecast times for this particular AnEn/
          past.fcst.day <- sapply(past.fcst.day, function(x) {which(
            config$forecast_times[x] == config$search_times_compare)})

          # If any of the time is not found, the returned variable should be a list
          if (is.list(past.fcst.day)) {
            stop('Some AnEn search times cannot be found in config.')
          }

          # Get the past time index in observation time for most similar past AnEn
          v <- AnEn.LOO.test$analogs[i.grid, past.fcst.day, i.flt, , 3]

          # Convert the matrix to vector using column-wise order
          v <- as.vector(v)

          if (avoid.duplicates) {
            # Remove duplicates. I didn't use the function unique() because I want
            # to keep the original order of these days.
            #
            v <- v[!duplicated(v)]
          }

          if (length(v) < 1 + insert.offset) {
            stop('There are not enough unique similar days to select from. Try adjusting scale.size and config$num_members.')
          }

          analogs.boost[i.grid, i.test, i.flt, insert.start:insert.end, 3] <- v[1:(1 + insert.offset)]
          analogs.boost[i.grid, i.test, i.flt, insert.start:insert.end, 2] <- i.grid
          analogs.boost[i.grid, i.test, i.flt, insert.start:insert.end, 1] <-
            config$search_observations[config$observation_id, i.grid, v[1:(1 + insert.offset)]]
        }

        if (show.progress) {
          counter <- counter + 1
          setTxtProgressBar(pb, counter)
        }

      }
    }
  }

  if (na.exists & !silent) {
    warning('NA values found in the original AnEn during boosting.')
  }

  if (show.progress) {
    close(pb)
  }

  # Garbage collection
  rm(AnEn, AnEn.LOO.test, config)
  garbage <- gc(reset = T)


  if (!silent) cat('Done (boostEnsembleSize)!\n')

  return(analogs.boost)
}
