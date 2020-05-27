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
#' EITrans::member_offset calculates an value offset for each member of the input
#' ensembles so that the resulting rank histogram shall be flat.
#'
#' @author Weiming Hu <weiming@@psu.edu>
#'
#' @details
#' Input ensembles and observations should match up for the first three dimensions
#' which should be #stations, #times, and #lead times. The 4th dimension of ensembles
#' should be ensemble members.
#'
#' The length of the output of the function will be a vector of the length of the
#' ensemble members. When you actually use these offset values to correct new ensembles,
#' **make sure that you SORT your ensemble members from smallest to largest**.
#'
#' @param ensembles A 4-dimensional array with the dimensions
#' `[stations, times, forecast lead times, ensemble members]`
#' @param observations A 3-dimensional array with the dimensions
#' `[stations, times, forecast lead times]`
#' @param left_delta A scalar addition offset for the minimum during sampling
#' @param right_delta A scalar addition offset for the maximum during sampling
#' @param infinity_estimator A scalar to scaling the spread of the ensemble
#' @param verbose Whether to be verbose
#' @param pre_sorted Whether the input ensemble members are pre-sorted. Use pre-sorted
#' ensembles will save you about 80% of the execution time. To pre-sort your ensembles,
#' use `aperm(apply(ensembles, 1:3, sort, na.last = T), c(2, 3, 4, 1))`.
#'
#' @return A vector for ensemble member offset. **Make sure that you have sorted
#' your ensemble values before you add the offset to each member**.
#'
#' @import RAnEnExtra
#' @md
#' @export
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
