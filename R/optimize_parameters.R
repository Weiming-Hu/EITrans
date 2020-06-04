#' EITrans::optimize_parameters
#'
#' EITrans::optimize_parameters carries out a search for the optimized parameter
#' combination for EITrans calibration with the input training ensembles forecasts
#' and observations.
#'
#' @author Weiming Hu <weiming@@psu.edu>
#'
#' @param grid_search A data frame with columns `left_deltas`, `right_deltas`,
#' `left_infinity`, `right_infinity`, and `multiplier`.
#' @param ens_train A 4-dimensional array as training ensemble forecasts.
#' Dimensions should be `[stations, times, lead times, members]`.
#' @param obs_train A 3-dimensional array for observations that correspond to the
#' ensemble forecasts. Dimensions should be `[stations, times, lead times]`.
#' @param verbose Whether to print progress information.
#'
#' @return A list with member offsets calculated with each of the parameter combination
#' and the rank histogram of the calibrated training forecast ensemble with the
#' corresponding parameter combination.
#'
#' @import progress
#' @import foreach
#'
#' @md
#' @export
optimize_parameters <- function(grid_search, ens_train, obs_train, verbose = F) {

  # Sanity check
  stopifnot(is.data.frame(grid_search))
  stopifnot(names(grid_search) %in% c('left_deltas', 'right_deltas',
                                      'left_infinity', 'right_infinity',
                                      'multiplier'))

  stopifnot(length(dim(ens_train)) == 4)
  stopifnot(length(dim(obs_train)) == 3)
  stopifnot(all.equal(dim(ens_train)[-4], dim(obs_train)))

  # Initialize a progress bar
  num_options <- nrow(grid_search)
  pb <- progress_bar$new(format = "[:bar] :percent eta: :eta",
                         total = num_options, clear = F)

  if (verbose) opts <- list()
  else opts <- list(progress = function(n) pb$tick())

  # Grid search for the best parameter combination
  if (verbose) cat('Grid search for best parameter combination out of',
                   num_options, 'possibilities ...\n')

  results <- foreach(index = 1:num_options, .options.snow = opts) %dopar% {

    # Calculate member offset
    offset <- member_offset(
      ensembles = ens_train,
      observations = obs_train,
      left_delta = grid_search$left_deltas[index],
      right_delta = grid_search$right_deltas[index],
      left_infinity_estimator = grid_search$left_infinity[index],
      right_infinity_estimator = grid_search$right_infinity[index],
      verbose = F, pre_sorted = T)

    # Apply the scaling factor
    offset <- offset * grid_search$multiplier[index]

    # Apply offset on the training dataset
    ens_train_calibrated <- apply_offset(
      ens = ens_train, offset = offset,
      pre_sorted = T, verbose = F)

    # Calculate rank histogram of the calibrated ensembles
    rh_dev_calibrated <- RAnEnExtra::verifyRankHist(
      anen.ver = ens_train_calibrated, obs.ver = obs_train,
      show.progress = F, pre.sort = T)

    list(offset = offset, rank = rh_dev_calibrated$rank)
  }

  # Check for any errors
  if (any(sapply(results, inherits, what = 'try-error'))) {
    stop('Grid search failed')
  }

  return(results)
}
