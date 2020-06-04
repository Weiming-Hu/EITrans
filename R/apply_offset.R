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

#' EITrans::apply_offset
#'
#' EITrans::apply_offset adds the offset values to a forecast ensemble. Offset values
#' are calculated from EITrans::member_offset. Or you could simply use EITrans::EITrans
#' which has built in a complete workflow for EITrans calibration.
#'
#' @param ens A 4-dimensional array for ensemble forecasts. Dimensions should be
#' `[stations, times, lead times, members]`.
#' @param offset Member offset values used to calibrate ensembles.
#' @param pre_sorted Whether the ensemble members are pre-sorted from the lowest to
#' the highest. Using pre-sorted ensembles can save runtime.
#' @param verbose Whether to print messages.
#'
#' @return A calibrated ensemble
#'
#' @md
#' @export
apply_offset <- function(ens, offset, pre_sorted = F, verbose = F) {

  # Sanity check
  stopifnot(length(dim(ens)) == 4)
  stopifnot(dim(ens)[4] == length(offset))

  if (!pre_sorted) {
    if (verbose) cat('Sort ensemble forecast memebers ...\n')
    ens <- aperm(apply(ens, 1:3, sort, na.last = T), c(2, 3, 4, 1))
  }

  if (verbose) cat('Calibrate ensembles ...\n')
  num_members <- dim(ens)[4]

  for (member_index in 1:num_members) {
    ens[, , , member_index] <-
      ens[, , , member_index, drop = F] +
      offset[member_index]
  }

  return(ens)
}
