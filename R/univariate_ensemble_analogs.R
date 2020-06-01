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

#' EITrans::univariate_ensemble_analogs
#'
#' EITrans::univariate_ensemble_analogs generates analogs using univariate ensemble
#' forecasts. It only returns the time indices where the most similar historical
#' ensemble forecasts are found.
#'
#' @author Weiming Hu <weiming@@psu.edu>
#'
#' @param ens A 4-dimensional array for ensemble forecasts. Dimensions should be
#' `[stations, times, lead times, members]`.
#' @param ens_times A vector for ensemble forecast times.
#' @param ens_flts A vector for ensemble forecast lead times.
#' @param ens_times_train Training times from the ensemble forecast times.
#' @param ens_times_dev Development times from the ensemble forecast times.
#' @param circular_ens Whether the ensemble forecast variable is circular.
#' @param member_weights Weights for each ensemble members when finding similar
#' historical ensemble forecasts.
#'
#' @import RAnEn
#' @import foreach
#' @import progress
#'
#' @md
#' @export
univariate_ensemble_analogs <- function(ens, ens_times, ens_flts,
                                        ens_times_train, ens_times_dev,
                                        circular_ens = F,
                                        member_weights = NULL) {

  # Sanity checks
  stopifnot(length(dim(ens)) == 4)
  stopifnot(is.logical(circular_ens))
  stopifnot(dim(ens)[2] == length(ens_times))
  stopifnot(dim(ens)[3] == length(ens_flts))
  stopifnot(all(ens_times_train %in% ens_times))
  stopifnot(all(ens_times_dev %in% ens_times))
  if (!is.null(member_weights)) stopifnot(
    is.numeric(member_weights) & length(member_weights) == dim(ens)[4])

  # Initialize a progress bar
  pb <- progress_bar$new(
    format = "[:bar] :percent eta: :eta",
    total = dim(ens)[1], clear = F)

  opts <- list(progress = function(n) pb$tick())

  # Generate analogs for each station separately
  results <- foreach(station_index = 1:dim(ens)[1],
                     .options.snow = opts) %dopar% {

    # Take care of threads
    num_threads <- RAnEn::getNumThreads()
    if (foreach::getDoParRegistered()) RAnEn::setNumThreads(1)

    # Convert ensembles to RAnEn::Forecasts
    fcsts <- RAnEn::generateForecastsTemplate()
    fcsts$Data <- aperm(ens[station_index, , , , drop = F], c(4, 1, 2, 3))
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
    config$verbose <- 1

    # Find similar ensemble forecasts from the training period
    AnEn <- RAnEn::generateAnalogs(fcsts, placeholder_obs,
                                   ens_times_dev, ens_times_train, config)

    # Revert thread control
    RAnEn::setNumThreads(num_threads)

    return(AnEn$similarity_time_index)
  }

  return(list(similarity_time_index = abind::abind(results, along = 1)))
}
