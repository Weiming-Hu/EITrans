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

#' EITrans::check_delta
#'
#' EITrans::check_delta checks whether the left/right deltas are valid EITrans arguments.
#'
#' @author Weiming Hu <weiming@@psu.edu>
#'
#' @param left_delta A left-edge delta value to validate
#' @param left_delta A left-edge delta value to validate
#' @param num_members The number of ensemble members
#' @return NULL if validate; otherwise an error will be raised.
#'
#' @md
#' @export
check_delta <- function(left_delta, right_delta, num_members) {
  sample <- seq(0 + left_delta, 1 + right_delta, length.out =  num_members + 2)
  sample <- sample[2:(num_members+1)]
  if (any(sample <= 0)) stop(paste('Left delta (', left_delta, ') too small!'))
  if (any(sample >= 1)) stop(paste('Right delta (', right_delta, ') too large!'))
  return(NULL)
}
