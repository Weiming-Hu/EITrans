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

#' EITrans::EITrans
#'
#' EITrans::EITrans is the main function for ensemble forecast calibration.
#'
#' @author Weiming Hu <weiming@@psu.edu>
#'
#' @param ens A 4-dimensional array for ensemble forecasts. Dimensions should be
#' `[stations, times, lead times, members]`.
#' @param ens_times A vector for ensemble forecast times.
#' @param ens_flts A vector for ensemble forecast lead times.
#' @param ens_times_train Training times from the ensemble forecast times.
#' @param ens_times_dev Development times from the ensemble forecast times.
#' @param ens_times_test Testing times from the ensemble forecast times.
#' @param obs A 3-dimensional array for observations that correspond to the
#' ensemble forecasts. Dimensions should be `[stations, times, lead times]`.
#' @param left_deltas A vector of left-edge deltas to experiment.
#' @param rigth_deltas A vector of right-edge deltas to experiment.
#' @param infinity_estimator A vector of values to experiment for estimating
#' the ensemble spread.
#' @param multiplier A vector of values to experiment for adjusting the
#' ensemble member offset.
#' @param circular_ens Whether the ensemble forecast variable is circular.
#' @param member_weights Weights for each ensemble members when finding similar
#' historical ensemble forecasts.
#'
#' @return A list with the calibrated ensemble forecasts and intermediate results.
#'
#' @import foreach
#'
#' @examples
#'
#' \dontrun{
#'
#' # If you are using MPI. Remember that you need to lauch this program with mpirun.
#' # I didn't have too much luck with the MPI spawn mechannism due to the hanging problem
#' # from Rmpi.
#' #
#' cl <- doMPI::startMPIcluster()
#' doMPI::registerDoMPI(cl)
#'
#' # If you are using doSNOW
#' cl <- snow::makeCluster()
#' doSNOW::registerDoSNOW(cl)
#'
#' eitrans_results <- EITrans(
#'   ens = ens$analogs,
#'   ens_times = ens$test_times,
#'   ens_flts = ens$flts,
#'
#'   ens_times_train = ens$test_times[1:(test_start - 366)],
#'   ens_times_dev = ens$test_times[(test_start - 365):(test_start - 1)],
#'   ens_times_test = ens$test_times[test_start:test_end],
#'
#'   obs = ens$obs_aligned,
#'
#'   left_deltas = seq(-0.046, 0.12, by = 0.002),
#'   right_deltas = seq(-0.02, 0.044, by = 0.002),
#'   infinity_estimator = seq(0.1, 0.5, by = 0.1),
#'   multiplier = seq(0.1, 1.1, by = 0.1))
#'
#' # If you are using MPI
#' doMPI::closeCluster(cl)
#' Rmpi::mpi.exit()
#'
#' # If you are using doSNOW
#' snow::stopCluster(cl)
#' }
#'
#' @md
#' @export
EITrans <- function(ens, ens_times, ens_flts,
                    ens_times_train, ens_times_dev, ens_times_test,
                    obs, left_deltas, right_deltas,
                    infinity_estimator,
                    multiplier,
                    circular_ens = F,
                    member_weights = NULL) {

  # Sanity checks
  cat('Start EITrans calibration ...\n')
  stopifnot(packageVersion('RAnEn') >= '4.0.8')
  stopifnot(all(ens_times_test %in% ens_times))
  stopifnot(length(intersect(ens_times_test, ens_times_dev)) == 0)
  stopifnot(length(intersect(ens_times_test, ens_times_train)) == 0)

  stopifnot(length(dim(obs)) == 3)
  stopifnot(all.equal(dim(obs), dim(ens)[-4]))

  stopifnot(infinity_estimator > 0)

  grid_search <- expand.grid(left_deltas = left_deltas,
                             right_deltas = right_deltas,
                             infinity_estimator = infinity_estimator,
                             multiplier = multiplier)

  lapply(1:nrow(grid_search), function(i) {
    check_delta(grid_search$left_deltas[i],
                grid_search$right_deltas[i],
                dim(ens)[4])})

  if (!getDoParRegistered()) stop(
    'Register your workers (doMPI for multinode and doSNOW for multicores) for parallel processing!')


  #################################################################################
  #                             Training Stage                                    #
  #################################################################################


  #
  # 1. Identify similar ensemble forecasts for ens_times_dev using ens_times_train
  #

  cat('Pre-sort ensemble forecast memebers ...\n')
  ens <- aperm(apply(ens, 1:3, sort, na.last = T), c(2, 3, 4, 1))

  cat('Identify similar ensemble forecasts for the dev period ...\n')
  dev_AnEn <- univariate_ensemble_analogs(
    ens = ens,
    ens_times = ens_times,
    ens_flts = ens_flts,
    ens_times_train = ens_times_train,
    ens_times_dev = ens_times_dev,
    circular_ens = circular_ens,
    member_weights = member_weights)


  #
  # 2. Prepare variables for grid search
  #

  cat('Prepare for grid search ...\n')

  ens_dev <- ens[, (ens_times %in% ens_times_dev), , , drop = F]
  obs_dev <- obs[, (ens_times %in% ens_times_dev), , drop = F]


  #
  # 3. Grid search for each forecast lead time
  #

  cat('Grid search with', nrow(grid_search), 'combinations for',
      dim(ens)[3],  'forecast lead times ...\n')

  best_combinations <- list()

  for (flt_index in 1:dim(ens)[3]) {

    cat('Process forecast lead time', flt_index, '/', dim(ens)[3], '...\n')

    rh_dev <- RAnEnExtra::verifyRankHist(
      anen.ver = ens_dev[, , flt_index, , drop = F],
      obs.ver = obs_dev[, , flt_index, drop = F],
      show.progress = F, pre.sort = T)

    results <- foreach(index = 1:nrow(grid_search)) %dopar% {

      # Initialization
      ens_similar <- array(NA, dim = dim(ens_dev[, , flt_index, , drop = F]))
      obs_similar <- array(NA, dim = dim(obs_dev[, , flt_index, drop = F]))

      # Use the most similar ensembles
      for (i in 1:dim(ens_similar)[1]) {
        for (j in 1:dim(ens_similar)[2]) {
          most_similar <- dev_AnEn$similarity_time_index[i, j, flt_index, 1]
          ens_similar[i, j, 1, ] <- ens[i, most_similar, 1, ]
          obs_similar[i, j, 1] <- obs[i, most_similar, 1]
        }
      }

      # Calculate member offset
      offset <- member_offset(
        ensembles = ens_similar,
        observations = obs_similar,
        left_delta = grid_search$left_deltas[index],
        right_delta = grid_search$right_deltas[index],
        infinity_estimator = grid_search$infinity_estimator[index],
        verbose = F, pre_sorted = T)

      # Apply the scaling factor
      offset <- offset * grid_search$multiplier[index]

      # Calculate quality of this offset
      ens_dev_calibrated <- aperm(
        apply(ens_dev[, , flt_index, , drop = F], 1:3,
              function(x) x + offset),
        c(2, 3, 4, 1))

      rh_dev_calibrated <- RAnEnExtra::verifyRankHist(
        anen.ver = ens_dev_calibrated,
        obs.ver = obs_dev[, , flt_index, drop = F],
        show.progress = F, pre.sort = T)

      list(offset = offset,
           rank = rh_dev_calibrated$rank,
           left_delta = grid_search$left_deltas[index],
           right_delta = grid_search$right_deltas[index],
           infinity_estimator = grid_search$infinity_estimator[index],
           multiplier = grid_search$multiplier[index])
    }

    if (any(sapply(results, inherits, what = 'try-error'))) {
      warning('Grid search failed. Error messages are returned')
      return(results)
    }

    # Idenfity best offset with the lowest sd
    best_index <- which.min(sapply(results, function(x) {sd(x$rank)}))

    # Store this best combination
    best_combinations[[flt_index]] <- list(
      offset = results[[best_index]]$offset,
      index = best_index,
      rh_dev = rh_dev,
      candidates = results)
  }


  #################################################################################
  #                             Testing Stage                                     #
  #################################################################################

  #
  # 4. Prepare variables for testing
  #

  ens_test <- ens[, (ens_times %in% ens_times_test), , , drop = F]

  #
  # 5. Calibrate test ensembles
  #

  cat('Calibrating test ensembles ...\n')

  for (flt_index in 1:dim(ens)[3]) {
    ens_test[, , flt_index, ] <- aperm(
      apply(ens_test[, , flt_index, , drop = F], 1:3,
            function(x) x + best_combinations[[flt_index]]$offset),
      c(2, 3, 4, 1))
  }

  calibration_results <- list(
    ens_test = ens_test,
    best_combinations_by_flt = best_combinations,
    grid_search = grid_search)
  class(calibration_results) <- 'AnEn'

  cat('EITrans Done!\n')
  return(calibration_results)
}
