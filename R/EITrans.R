# This is the EITrans main function script.

EITrans <- function(ens, ens_times, ens_flts,
										ens_times_train, ens_times_dev, ens_times_test,
										obs, left_deltas, right_deltas,
										infinity_estimator,
										multiplier,
										circular_ens = F,
										member_weights = NULL,
										grid_search_cores = 1) {
	
	# Sanity checks
	stopifnot(length(dim(ens)) == 4)
	stopifnot(dim(ens)[2] == length(ens_times))
	stopifnot(dim(ens)[3] == length(ens_flts))
	
	stopifnot(all(ens_times_train %in% ens_times))
	stopifnot(all(ens_times_dev %in% ens_times))
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
	
	
	#################################################################################
	#                             Training Stage                                    #
	#################################################################################
	
	
	#
	# 1. Identify similar ensemble forecasts for ens_times_dev using ens_times_train
	#
	
	cat('Pre-sort ensemble forecast memebers ...\n')
	ens <- aperm(apply(ens, 1:3, sort, na.last = T), c(2, 3, 4, 1))
	
	cat('Identify similar ensemble forecasts for the dev period ...\n')
	
	# Convert ensembles to RAnEn::Forecasts
	fcsts <- RAnEn::generateForecastsTemplate()
	fcsts$Data <- aperm(ens, c(4, 1, 2, 3))
	fcsts$ParameterNames <- paste('member', 1:dim(fcsts$Data)[1], sep = '_')
	if (circular_ens) fcsts$ParameterCirculars <- fcsts$ParameterNames
	fcsts$Times <- ens_times
	fcsts$FLTs <- ens_flts
	
	# Generate observation placeholder
	placeholder_obs <- RAnEn::generateObservationsTemplate()
	placeholder_obs$ParameterNames <- 'Placeholder'
	placeholder_obs$Times <- unique(rep(ens_times, each = length(ens_flts)) + 
																		as.numeric(ens_flts))
	placeholder_obs$Data <- array(1, dim = c(1, dim(fcsts$Data)[2],
																					 length(placeholder_obs$Times)))
	
	# Set up config
	config <- new(RAnEn::Config)
	config$num_similarity <- 1
	config$quick <- F
	config$operation <- F
	config$save_analogs <- F
	config$save_analogs_time_index <- F
	config$save_similarity <- F
	config$save_similarity_time_index <- T
	if (!is.null(member_weights)) config$weights <- member_weights
	config$verbose <- 3
	
	# Find similar ensemble forecasts from the training period
	dev_AnEn <- RAnEn::generateAnalogs(fcsts, placeholder_obs,
																		 ens_times_dev, ens_times_train, config)
	
	# Housekeeping
	rm(config, placeholder_obs)
	
	
	#
	# 2. Prepare variables for grid search
	#
	
	cat('Prepare for grid search ...\n')
	
	ens_train <- ens[, (ens_times %in% ens_times_train), , , drop = F]
	obs_train <- obs[, (ens_times %in% ens_times_train), , drop = F]
	
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
		
		results <- pbmcapply::pbmclapply(1:nrow(grid_search), function(index) {
			
			ens_similar <- array(NA, dim = dim(ens_dev[, , flt_index, , drop = F]))
			obs_similar <- array(NA, dim = dim(obs_dev[, , flt_index, drop = F]))
			
			# Use the most similar ensembles
			for (i in 1:dim(ens_similar)[1]) {
				for (j in 1:dim(ens_similar)[2]) {
					most_similar <- dev_AnEn$similarity_time_index[i, j, flt_index, 1]
					ens_similar[i, j, 1, ] <- ens_train[i, most_similar, 1, ]
					obs_similar[i, j, 1] <- obs_train[i, most_similar, 1]
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
			
			list(offset = offset, rank = rh_dev_calibrated$rank,
					 left_delta = grid_search$left_deltas[index],
					 right_delta = grid_search$right_deltas[index],
					 infinity_estimator = grid_search$infinity_estimator[index],
					 multiplier = grid_search$multiplier[index])},
			
			mc.cores = grid_search_cores)
		
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
	
	cat('EITrans Done!\n')
	return(list(ens_test = ens_test,
							best_combinations_by_flt = best_combinations,
							grid_search = grid_search))
}

check_delta <- function(left_delta, right_delta, num_members) {
	sample <- seq(0 + left_delta, 1 + right_delta, length.out =  num_members + 2)
	sample <- sample[2:(num_members+1)]
	if (any(sample <= 0)) stop(paste('Left delta (', left_delta, ') too small!'))
	if (any(sample >= 1)) stop(paste('Right delta (', right_delta, ') too large!'))
	return(NULL)
}

member_offset <- function(ensembles, observations,
													left_delta, right_delta,
													infinity_estimator = 1,
													verbose = T, pre_sorted = F) {
	
	# Sanity check
	stopifnot(length(dim(ensembles)) == 4)
	stopifnot(length(dim(observations)) == 3)
	stopifnot(all.equal(dim(ensembles)[1:3], dim(observations)))
	if (!pre_sorted) stop('Must presort ensemble members!')
	
	# Define the problem dimensions
	num_stations <- dim(ensembles)[1]
	num_times <- dim(ensembles)[2]
	num_flts <- dim(ensembles)[3]
	num_members <- dim(ensembles)[4]
	
	# Sample from the cumulative probability space
	sample <- seq(0 + left_delta, 1 + right_delta, length.out =  num_members + 2)
	sample <- sample[2:(num_members+1)]
	if (any(sample <= 0)) return('Left delta too small!')
	if (any(sample >= 1)) return('Right delta too large!')
	
	# Generating rank histogram
	if (verbose) cat('Generating rank histogram ...\n')
	rh <- RAnEnExtra::verifyRankHist(anen.ver = ensembles,
																	 obs.ver = observations,
																	 show.progress = verbose,
																	 pre.sort = pre_sorted)
	cum_sum <- c(0, cumsum(rh$rank))
	
	# Calculate offset from the observation value
	for (member_i in 1:dim(ensembles)[4]) {
		ensembles[, , , member_i] <- abind::adrop(
			ensembles[, , , member_i, drop = F], drop = 4) - observations
	}
	
	# Calculate the average offset from the observation values for all ensembles
	if (verbose) cat('Calculating average offset of members ...\n')
	average_offset <- apply(ensembles, 4, mean, na.rm = T)
	average_offset <- c(
		
		# The left-most value is to represent the smallest possible value (hence -Inf)
		average_offset[1] - infinity_estimator * (mean(average_offset) - average_offset[1]),
		
		# The original values
		average_offset,
		
		# The right-most value is to represent the largest possible value (hence Inf)
		average_offset[num_members] + infinity_estimator * (
			average_offset[num_members] - mean(average_offset)))
	
	# Create the new offset
	resampled_offset <- sapply(sample, function(x) {
		vec <- x - cum_sum
		i <- which(vec < 0)[1]
		
		approx(x = cum_sum[c(i-1, i)], y = average_offset[c(i-1, i)], xout = x)$y
	})
	
	ensemble_offset <- resampled_offset - average_offset[-c(1, num_members + 2)]
	
	if (verbose) cat('Done!\n')
	return(ensemble_offset)
}
