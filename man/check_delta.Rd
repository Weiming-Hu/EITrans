% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/check_delta.R
\name{check_delta}
\alias{check_delta}
\title{EITrans::check_delta}
\usage{
check_delta(left_delta, right_delta, num_members)
}
\arguments{
\item{left_delta}{A left-edge delta value to validate}

\item{num_members}{The number of ensemble members}
}
\value{
NULL if validate; otherwise an error will be raised.
}
\description{
EITrans::check_delta checks whether the left/right deltas are valid EITrans arguments.
}
\author{
Weiming Hu \href{mailto:weiming@psu.edu}{weiming@psu.edu}
}
