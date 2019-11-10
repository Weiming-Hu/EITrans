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

#' EITrans::heuristicFilter
#'
#' EITrans::heuristicFilter is designed to improve the reliability of
#' Analog Ensemble forecasts by using the leave-one-out method and the
#' inverse transformation function. The filter works better with ensembles
#' with large amount of members (usually at the level of a hunderd).
#'
#' @author Weiming Hu \email{weiming@@psu.edu}
#'
#' @param AnEn The AnEn results from \code{\link{generateAnalogs}}.
#' @param config The configuration used by \code{\link{generateAnalogs}}.
#' @param final.ensemble.size How many ensemble members to keep in the
#' final ensemble forecasts.
#' @param historical.similar.count The number of past similar forecasts
#' to include when considering calculating the historical rank histogram.
#' @param LOO.test.size The number of historical forecasts to
#' carry out leave-one-out test.
#' @param keep.table Whether to keep the summary table of similariy.
#' @param keep.LOO.rank Whether to keep the verification rank histogram
#' for the leave-one-out tests.
#' @param keep.LOO Whether to keep the analogs and similarity results
#' from the leave-one-out test if they are computed.
#' @param do.not.append If set to True, results will not be appended to
#' AnEn, and the filtered analog ensembles will be returned. At the
#' same time, all `keep*` parameters become invalid.
#' @param show.progress Whether to show the progress bar.
#' @param silent Whether to be silent.
#' @param i.station A vector of station indices of which analogs will be filtered.
#' @param i.test.day A vector of test day indices of which analogs will be filtered.
#' @param AnEn.LOO.test The precomputed AnEn results for leave-one-out tests.
#' If not provided, LOO tests will be carried out at the spot. See details for
#' instructions on precopmuting LOO tests.
#' @param  heuristic.ensemble.size The number of ensemble members to take from
#' the LOO test results to construct the heuristic rank histogram. By default,
#' it is the `num_members` member from input `config`.
#' @param member.name The name of the member in AnEn to process. By default, it
#' is `analogs`. But this parameter can be helpful if, for example, you want to
#' process another member called `analogs.alternative`.
#'
#' @details
#' To precompute AnEn.LOO.test, simply copy the values from `search_times_compare`
#' to `test_times_compare` in the `config` which is passed to this function, and
#' use the `generateAnalogs` to generate the `AnEn.LOO.test`. To save memory and space,
#' you can optionally omit similarity information in the results by setting
#' `config$keep_similarity = F`. **The only required member is `analogs`**.
#'
#' You might wonder when you need to generate an updated version of `AnEn.LOO.test`.
#' Usually, if you change the settings in `config` which would affect the search
#' data, for example, `search_times_compare`, and the calculation of `AnEn`, for
#' example, `weights`, you will need to regenerate the `AnEn.LOO.test`. Otherwise, you
#' do not need to regenerate it. For example, you don't need to regenerate `AnEn.LOO.test`
#' if you only change `test_times_compare`.
#'
#' @md
#' @export
heuristicFilter <- function (AnEn, config, final.ensemble.size,
                             historical.similar.count = 1,
                             LOO.test.size = 10,
                             keep.table = T, keep.LOO.rank = T,
                             keep.LOO = T, do.not.append = F,
                             show.progress = T, silent = F,
                             i.station = NA, i.test.day = NA,
                             AnEn.LOO.test = NA,
                             heuristic.ensemble.size = config$num_members,
                             member.name = 'analogs') {

  # Sanity checks
  stopifnot(!config$quick)
  stopifnot(config$preserve_similarity)
  stopifnot(!config$advanced)
  stopifnot(config$mode == 'independentSearch')

  if (identical(i.test.day, NA)) {
    i.test.day <- 1:dim(AnEn[[member.name]])[2]
  }

  if (identical(i.station, NA)) {
    i.station <- 1:dim(AnEn[[member.name]])[1]
    config$forecasts <- config$forecasts[, i.station, , , drop = F]
    config$search_observations <- config$search_observations[, i.station, , drop = F]
  }

  # Select the most similar past forecasts to build heuristic ranks
  if (!silent) cat('Building heuristic ranks based on similar past forecsats ...\n')

  # Counts the number of times each past forecast is selected
  summary <- tabulate(AnEn$similarity[
    i.station, i.test.day, , 1:historical.similar.count, 3])
  names(summary) <- 1:length(summary)
  summary <- sort(summary, decreasing = T, na.last = T)
  stopifnot(LOO.test.size <= length(summary))

  # Get the most similar past forecast days for the leave-one-out tests
  if (LOO.test.size > length(summary)) {
    warning('LOO.test.size exceeds the number of available days. Set to the maximum.')
    LOO.test.size <- length(summary)
  }

  i.LOO.test.days <- as.numeric(names(summary[1:LOO.test.size]))

  config$test_times_compare <- config$forecast_times[i.LOO.test.days]


  if (identical(AnEn.LOO.test, NA)) {
    # No precomputed LOO provided
    if (!silent) cat('Compute LOO tests ...\n')

    # Leave-one-out test
    config$num_members <- heuristic.ensemble.size
    config$max_num_sims <- heuristic.ensemble.size
    config$preserve_mapping <- F
    AnEn.LOO.test <- generateAnalogs(config)

  } else {
    # Precomputed LOO provided
    if (!silent) cat('Use precomputed LOO test results ...\n')

    stopifnot(dim(AnEn.LOO.test$analogs)[1] == dim(AnEn[[member.name]])[1])
    stopifnot(dim(AnEn.LOO.test$analogs)[2] == length(config$search_times_compare))
    stopifnot(dim(AnEn.LOO.test$analogs)[3] == dim(AnEn[[member.name]])[3])
    stopifnot(dim(AnEn.LOO.test$analogs)[4] >= heuristic.ensemble.size)

    AnEn.LOO.test$analogs <- AnEn.LOO.test$analogs[
      i.station, i.LOO.test.days, , 1:heuristic.ensemble.size, , drop = F]
  }

  # Generate verifications for leave-one-out tests
  if (!silent) cat('Querying observations ...\n')
  obs.all <- alignObservations(
    config$search_observations, config$observation_times,
    config$test_times_compare, config$flts,
    show.progress = show.progress, silent = silent)

  # Reduce the observation dimensions
  obs <- obs.all[config$observation_id, , , , drop = F]
  dim(obs) <- dim(obs)[-1]

  # Reduce the analog dimensions
  anen <- AnEn[[member.name]][i.station, i.test.day, , , 1, drop = F]
  dim(anen) <- dim(anen)[-5]

  anen.LOO <- AnEn.LOO.test$analogs[, , , , 1, drop = F]
  dim(anen.LOO) <- dim(anen.LOO)[-5]

  # Generate heuristic ranks for leave-one-out analogs
  if (!silent) cat('Generating heuristic rank hitograms ...\n')
  rh.LOO.test <- verifyRankHist(anen.LOO, obs, show.progress = show.progress)

  # Heuristic filter
  if (!silent) cat('Applying the inverse transform function ...\n')
  analogs.hf <- inverseTransform(rh.LOO.test$rank, anen, final.ensemble.size)

  if (is.null(analogs.hf)) {
    return(NULL)
  }

  if (do.not.append) {
    return(analogs.hf)
  } else {
    AnEn$analogs.hf <- analogs.hf

    if (keep.table) {
      if (!silent) {
        cat('Appending the similarity summary table ...\n')
      }
      AnEn$count.table <- summary
    }

    if (keep.LOO.rank) {
      if (!silent) {
        cat('Appending the rank for leave-one-out tests ...\n')
      }
      AnEn$heuristic.rank <- rh.LOO.test
    }

    if (keep.LOO) {
      if (!silent) {
        cat('Appending the leave-one-out test results ...\n')
      }
      AnEn$analogs.LOO <- AnEn.LOO.test$analogs

      if ('similarity' %in% names(AnEn)) {
        AnEn$similarity.LOO <- AnEn.LOO.test$similarity
      }

      AnEn$LOO.forecast.times <- config$test_times_compare
    }

    if (!silent) {
      cat('Done (heuristicFilter)!\n')
    }
    return(AnEn)
  }
}
