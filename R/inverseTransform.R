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

#' EITrans::inverseTransform
#'
#' EITrans::inverseTransform resample the analog ensemble members
#' based on the distribution match information from the heuristic rank
#' histogram. This function keeps a subset of members from the original
#' forecast ensembles and the goal is to improve the reliability of
#' the ensemble subset.
#'
#' @author Weiming Hu \email{weiming@@psu.edu}
#'
#' @param heuristic.rank The frequency or count for each ranked bin.
#' An example of this parameter is the list member *rank* generated
#' using the function \code{\link{verifyRankHist}}.
#' @param analogs A 4-dimensional array for analog values. An exmple
#' of this parameter is the list member *analogs* with only the *value*
#' column generated using the function \code{\link{generateAnalogs}}.
#' @param members.to.keep The number of ensemble members to keep.
#' @param digits The number of decimals when computing the cumulous summation.
#'
#' @md
#' @export
inverseTransform <- function(
  heuristic.rank, analogs, members.to.keep, digits = 6) {

  # Sanity checks
  stopifnot(length(dim(analogs)) == 4)

  if (length(heuristic.rank) != dim(analogs)[4] + 1) {
    heuristic.rank <- approx(
      x = 1:length(heuristic.rank), y = heuristic.rank,
      xout = seq(from = 1, to = length(heuristic.rank),
                 length.out = dim(analogs)[4] + 1))$y

    heuristic.rank <- heuristic.rank/sum(heuristic.rank)
  }

  if (members.to.keep >= dim(analogs)[4]) {
    stop("Too many members to keep. Check your members.to.keep and analogs.")
  }

  analogs.order <- apply(analogs, c(1, 2, 3), sort, na.last = T)
  analogs.order <- aperm(analogs.order, c(2, 3, 4, 1))

  # Compute the cumulous rank values
  heuristic.rank.cumsum <- round(cumsum(heuristic.rank), digits = digits)

  # Equally sample
  sample <- seq(0, 1, length.out = members.to.keep + 2)
  sample <- sample[-c(1, members.to.keep + 2)]

  # Which bin does each sample correspond to
  selected <- sapply(sample, function(x) {
    vec <- abs(x - heuristic.rank.cumsum)
    i <- which(vec == min(vec))
  
    if (length(i) == 1) {return(i)}
    else {return(sample(i, 1))}
  })

  if (any(duplicated(selected))) {
    stop("Same ranks end up been selected multiple times. Results might be corrupted.")
  }

  stopifnot(length(selected) == members.to.keep)

  return(analogs.order[, , , selected, drop = F])
}
