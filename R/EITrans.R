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
#' @param member_weights Weights for each ensemble members when finding similar
#' historical ensemble forecasts.
#' @param pre_sorted Whether the ensemble members are presorted.
#' @param save_intermediate Whether to save intermediate data.
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
                    circular_ens = F, member_weights = NULL,
                    pre_sorted = F, save_intermediate = F) {


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
		circular_ens = circular_ens,
		member_weights = member_weights)

	# Initialize arrays for the most similar historical ensembles
	test_times_index <- ens_times %in% ens_times_test
	ens_train <- array(NA, dim = c(num_stations, num_test_times, num_flts, num_members))
	obs_train <- array(NA, dim = c(num_stations, num_test_times, num_flts))

	# Initilize a progress bar
	pb <- progress_bar$new(format = "[:bar] :percent eta: :eta",
												 total = num_stations * num_test_times, clear = F)

	# Convert test ensembles to the most similar historical ensembles
	cat('Index analog values ...\n')
	for (station_index in 1:num_stations) {
		for (time_index in 1:num_test_times) {
			for (flt_index in 1:num_flts) {
				most_similar <- AnEn$similarity_time_index[station_index, time_index, flt_index, 1]

				ens_train[station_index, time_index, flt_index, ] <- ens[station_index, most_similar, flt_index, ]
				obs_train[station_index, time_index, flt_index] <- obs[station_index, most_similar, flt_index]
			}

			pb$tick()
		}
	}


	#########################
	# Parameter grid search #
	#########################

	# Initialize a progress bar
	num_options <- nrow(grid_search)
	pb <- progress_bar$new(format = "[:bar] :percent eta: :eta", total = num_options, clear = F)
	opts <- list(progress = function(n) pb$tick())

	# Grid search for the best parameter combination
	cat('Grid search for best parameter combination out of', num_options, 'possibilities ...\n')
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

		# Calculate quality of this offset
		ens_train_calibrated <- apply(ens_train, 1:3, function(x) x + offset)

		# Permutate array dimensions to be consistent with the original ensemble
		ens_train_calibrated <- aperm(ens_train_calibrated, c(2, 3, 4, 1))

		# Calculate rank histogram of the calibrated ensembles
		rh_dev_calibrated <- RAnEnExtra::verifyRankHist(
			anen.ver = ens_train_calibrated, obs.ver = obs_train,
			show.progress = F, pre.sort = T)

		list(offset = offset, rank = rh_dev_calibrated$rank)
	}

	# Check for any errors
	if (any(sapply(results, inherits, what = 'try-error'))) {
		warning('Grid search failed. Error messages are retruned.')
		return(results)
	}

	# Initialize a list for return
	eitrans_results <- list()
	class(eitrans_results) <- 'AnEn'

	# Identify the best parameter combination index
	best_index <- which.min(sapply(results, function(x) {sd(x$rank)}))

	# Store this best combination
	best_offset <- results[[best_index]]$offset

	# Save results to the output list
	eitrans_results$best_offset <- best_offset
	eitrans_results$best_index <- best_index
	eitrans_results$grid_search <- grid_search

	if (save_intermediate) {
		# Save the rank histogram of the original train ensembles
		eitrans_results$train_rank_original <-
			RAnEnExtra::verifyRankHist(
				anen.ver = ens_train, obs.ver = obs_train,
				show.progress = F, pre.sort = T)$rank

		# Save the rank histograms of the calibrated train ensembles
		eitrans_results$train_rank_calibrated <- results
	}


	################################################################
	# Calibrate test ensembles with the best parameter combination #
	################################################################

	# Extract test ensembles
	ens_test <- ens[, (ens_times %in% ens_times_test), , , drop = F]

	# Calibrate test ensembles
	eitrans_results$analogs_calibrated <- apply_offset(
	  ens = ens_test,
	  offset = best_offset,
	  pre_sorted = T,
	  verbose = T)

	if (save_intermediate) {
		eitrans_results$analogs_most_similar <- ens_train
	}

	cat('EITrans Done!\n')
	return(eitrans_results)
}
