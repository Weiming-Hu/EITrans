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

#' EITrans::index_to_value
#'
#' EITrans::index_to_value
#'
#' @author Weiming Hu <weiming@@psu.edu>
#'
#' @import progress
#'
#' @md
#' @export
index_to_value <- function(similarity_time_index, ens, obs = NULL, verbose = F) {

  stopifnot(length(dim(ens)) == 4)
  stopifnot(length(dim(obs)) == 3)
  stopifnot(all.equal(dim(obs), dim(ens)[-4]))
  stopifnot(length(dim(similarity_time_index)) == 4)
  stopifnot(dim(similarity_time_index)[4] == 1)

  num_stations <- dim(similarity_time_index)[1]
  num_test_times <- dim(similarity_time_index)[2]
  num_flts <- dim(similarity_time_index)[3]
  num_members <- dim(ens)[4]

  ens_similar <- array(NA, dim = c(num_stations, num_test_times, num_flts, num_members))
  obs_similar <- array(NA, dim = c(num_stations, num_test_times, num_flts))

  # Initilize a progress bar
  pb <- progress_bar$new(format = "[:bar] :percent eta: :eta",
                         total = num_stations * num_test_times, clear = F)

  # Convert test ensembles to the most similar historical ensembles
  if (verbose) cat('Index analog values ...\n')
  for (station_index in 1:num_stations) {
    for (time_index in 1:num_test_times) {
      for (flt_index in 1:num_flts) {
        most_similar <- similarity_time_index[station_index, time_index, flt_index, 1]

        ens_similar[station_index, time_index, flt_index, ] <-
          ens[station_index, most_similar, flt_index, ]

        if (!is.null(obs)) {
          obs_similar[station_index, time_index, flt_index] <-
            obs[station_index, most_similar, flt_index]
        }
      }

      if (verbose) pb$tick()
    }
  }

  return(list(ens = ens_similar, obs = obs_similar))
}
