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
#' EITrans::EITrans is the main function for ensemble forecast calibration with the empirical
#' inverse transform function method.
#'
#' @details
#' This function uses `foreach` parallel mechanism. Parallelization is handled
#' by users creating the parallel backend. Please see examples. It is recommended
#' to use `doSNOW` for handling the parallel backend.
#'
#' @author Weiming Hu <weiming@@psu.edu>
#'
#' @param ens A 4-dimensional array for ensemble forecasts. Dimensions should be
#' `[stations, times, lead times, members]`.
#' @param ens_times A vector for ensemble forecast times.
#' @param ens_flts A vector for ensemble forecast lead times.
#' @param ens_times_train Training times from the ensemble forecast times.
#' @param ens_times_test Testing times from the ensemble forecast times.
#' @param obs A 3-dimensional array for observations that correspond to the
#' ensemble forecasts. Dimensions should be `[stations, times, lead times]`.
#' @param deltas A vector of deltas to experiment. It can also be a list with
#' `left` and `right` to explicitly control deltas for edges.
#' @param infinity_estimators A vector of values to experiment for estimating
#' the ensemble spread. Or a list with named members `left` and `right` for
#' different left/right infinity estimators.
#' @param multipliers A vector of values to experiment for adjusting the
#' ensemble member offset.
#' @param circular_ens Whether the ensemble forecast variable is circular.
#' @param pre_sorted Whether the ensemble members are presorted.
#' @param save_intermediate Whether to save intermediate data.
#' @param optimize_lead_time Whether to calibrate each forecast lead time
#' respectively. Theoretically, this would lead to better calibration results,
#' but also significantly more computation because calibration needs to be
#' evaluated per lead time.
#'
#' @return A list with the calibrated ensemble forecasts and intermediate results.
#'
#' @import foreach
#' @import progress
#'
#' @examples
#'
#' \dontrun{
#' # Use doSNOW to launch parallel backend
#' library(doSNOW)
#' cl <- snow::makeCluster(4) # Use 4 cores
#' registerDoSNOW(cl)
#'
#' # If you know a host file, you can also launch on remote server
#' nodefile <- Sys.getenv("PBS_NODEFILE")
#' nodes <- readLines(nodefile)
#' cl <- makeCluster(nodes, type = "SOCK")
#' registerDoSNOW(cl)
#'
#' eitrans_results <- EITrans(
#'   ens = ens$analogs,
#'   ens_times = ens$test_times,
#'   ens_flts = ens$flts,
#'
#'   ens_times_train = ens$test_times[1:(test_start - 366)],
#'   ens_times_test = ens$test_times[test_start:test_end],
#'
#'   obs = ens$obs_aligned,
#'
#'   deltas = seq(-0.02, 0.044, by = 0.002),
#'   infinity_estimators = seq(0.1, 0.5, by = 0.1),
#'   multipliers = seq(0.1, 1.1, by = 0.1))
#'
#' # If you are using doSNOW
#' stopCluster(cl)
#' }
#'
#' @md
#' @export
EITrans <- function(ens, ens_times, ens_flts,
                    ens_times_train, ens_times_test, obs,
                    deltas, infinity_estimators, multipliers,
                    circular_ens = F,
                    pre_sorted = F, save_intermediate = F,
                    optimize_lead_time = F) {


	#################
	# Sanity checks #
	#################

	cat('Start EITrans calibration ...\n')
	stopifnot(packageVersion('RAnEn') >= '4.0.8')

	stopifnot(length(dim(ens)) == 4)
	stopifnot(length(dim(obs)) == 3)
	stopifnot(dim(ens)[2] == length(ens_times))
	stopifnot(dim(ens)[3] == length(ens_flts))
	stopifnot(all.equal(dim(obs), dim(ens)[-4]))

	stopifnot(all(ens_times_test %in% ens_times))
	stopifnot(all(ens_times_train %in% ens_times))
	stopifnot(length(intersect(ens_times_test, ens_times_train)) == 0)

	if (inherits(deltas, 'list')) {
		stopifnot(all(names(deltas) %in% c('left', 'right')))
		left_deltas <- deltas$left
		right_deltas <- deltas$right

	} else {
		left_deltas <- right_deltas <- deltas
	}

	if (inherits(infinity_estimators, 'list')) {
		stopifnot(all(names(infinity_estimators) %in% c('left', 'right')))
		left_infinity <- infinity_estimators$left
		right_infinity <- infinity_estimators$right

	} else {
		left_infinity <- right_infinity <- infinity_estimators
	}

	stopifnot(all(multipliers > 0))
	stopifnot(all(left_infinity > 0))
	stopifnot(all(right_infinity > 0))

	grid_search <- expand.grid(left_deltas = left_deltas,
														 right_deltas = right_deltas,
														 left_infinity = left_infinity,
														 right_infinity = right_infinity,
														 multiplier = multipliers)

	lapply(1:nrow(grid_search), function(i) {
		check_delta(grid_search$left_deltas[i],
								grid_search$right_deltas[i],
								dim(ens)[4])})

	rm(right_infinity, left_infinity, left_deltas, right_deltas,
		 multipliers, deltas, infinity_estimators)

	if (!getDoParRegistered()) warning(
		'Create and register your workers with doSNOW for parallel processing!')

	# Get dimensions of ensembles
	num_stations <- dim(ens)[1]
	num_times <- dim(ens)[2]
	num_flts <- dim(ens)[3]
	num_members <- dim(ens)[4]
	num_test_times <- length(ens_times_test)


	#############################
	# Generate ensemble analogs #
	#############################

	# Pre-sort the ensembles to speed up computation
	if (!pre_sorted) {
		cat('Pre-sort ensemble forecast memebers ...\n')
		ens <- aperm(apply(ens, 1:3, sort, na.last = T), c(2, 3, 4, 1))
	}

	# Identify most similar historical ensembles for each of the test ensembles
	cat('Identify analog ensemble index for the test period ...\n')
	AnEn <- univariate_ensemble_analogs(
		ens = ens,
		ens_times = ens_times,
		ens_flts = ens_flts,
		ens_times_train = ens_times_train,
		ens_times_dev = ens_times_test,
		circular_ens = circular_ens)


	# Convert indices to values
	l <- index_to_value(AnEn$similarity_time_index, ens, obs)
	ens_train <- l$ens
	obs_train <- l$obs


	#########################
	# Parameter grid search #
	#########################

	# Initialize a list for return
	eitrans_results <- list(optimize_lead_time = optimize_lead_time)
	class(eitrans_results) <- 'AnEn'

	if (optimize_lead_time) {

	  #
	  # If lead time optimization is enabled, we need to apply the search across
	  # each lead time and come up with an individual set of offset values for members
	  # at each lead time. This would lead to significantly more computation but
	  # theoretically better calibration performance.
	  #

	  results <- list()

	  for (flt_index in 1:num_flts) {
	    cat('Grid search for lead time', flt_index, '/', num_flts, '...\n')

	    results[[flt_index]] <- optimize_parameters(
	      grid_search = grid_search,
	      ens_train = ens_train[, , flt_index, , drop = F],
	      obs_train = obs_train[, , flt_index, drop = F],
	      verbose = T)
	  }

	} else {

	  #
	  # If lead time optimization is not enabled, we could simply apply the search
	  # across all lead times and come up with a single set of offset values for
	  # members at all lead times.
	  #

	  results <- optimize_parameters(
	    grid_search = grid_search,
	    ens_train = ens_train,
	    obs_train = obs_train,
	    verbose = T)
	}

	# Define the metric for selecting the best offset values
	metric <- function(result) sd(result$rank)

	# Identify the best parameter combination based on the lowest standard
	# deviation of the rank histogram of the calibrated ensemble.
	#
	if (optimize_lead_time) {
	  best_index <- sapply(results, function(x) which.min(sapply(x, metric)))
	  best_offset <- sapply(1:length(best_index),
	                        function(index) results[[index]][[best_index[index]]]$offset)
	} else {
	  best_index <- which.min(sapply(results, metric))
	  best_offset <- results[[best_index]]$offset
	}

	# Save results to the output list
	eitrans_results$best_offset <- best_offset
	eitrans_results$best_index <- best_index
	eitrans_results$grid_search <- grid_search

	# Save the rank histogram of the original train ensembles
	if (save_intermediate) {
	  cat('Calculate the rank histogram of the training ensembles ...\n')

	  if (optimize_lead_time) {
	    eitrans_results$train_rank_original <- sapply(1:num_flts, function(flt_index) {
	      RAnEnExtra::verifyRankHist(
	        anen.ver = ens_train[, , flt_index, , drop = F],
	        obs.ver = obs_train[, , flt_index, drop = F],
	        show.progress = F, pre.sort = T)$rank
	    })

	  } else {
	    eitrans_results$train_rank_original <-
	      RAnEnExtra::verifyRankHist(
	        anen.ver = ens_train, obs.ver = obs_train,
	        show.progress = F, pre.sort = T)$rank
	  }

		# Save the rank histograms of the calibrated train ensembles
		eitrans_results$train_rank_calibrated <- results
	}


	################################################################
	# Calibrate test ensembles with the best parameter combination #
	################################################################

	# Extract test ensembles
	ens_test <- ens[, (ens_times %in% ens_times_test), , , drop = F]

	# Calibrate test ensembles
	if (optimize_lead_time) {
	  cat('Calibrate ensembles by individual lead times ...\n')

	  for (flt_index in 1:num_flts) {
	    ens_test[, , flt_index, ] <- apply_offset(
	      ens = ens_test[, , flt_index, , drop = F],
	      offset = best_offset[, flt_index],
	      pre_sorted = T,
	      verbose = F)
	  }

	} else {
	  ens_test <- apply_offset(
	    ens = ens_test,
	    offset = best_offset,
	    pre_sorted = T,
	    verbose = T)
	}

	eitrans_results$analogs_calibrated <- ens_test

	if (save_intermediate) {
		eitrans_results$analogs_most_similar <- ens_train
	}

	cat('EITrans Done!\n')
	return(eitrans_results)
}
