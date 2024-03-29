% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Utility.R
\name{prepwide}
\alias{prepwide}
\title{Prepares data in a wide format for [mcem].}
\usage{
prepwide(data, studentid, passageid, numwords.p, wrc, time)
}
\arguments{
\item{data}{A data frame.}

\item{studentid}{The column name in the data that represents the unique student identifier.}

\item{passageid}{The column name in the data that represents the unique passage identifier.}

\item{numwords.p}{The column name in the data that represents the number of words in a passage.}

\item{wrc}{The column name in the data that represents the words read correctly for each case.}

\item{time}{The column name in the data that represents the time, in seconds, for each case.}
}
\description{
This function will return a list with 5 elements:
Y: a matrix of words read correctly, where rows represent cases (student and occasion) and columns represent passages
logt10: a [tibble::tibble()] of words read correctly
N: the number of cases per passage
I: the number of passages
}
\examples{
data("passage")

prepwide(passage,
 studentid = "id.student",
 passageid = "id.passage",
 numwords.p = "numwords.pass",
 wrc = "wrc",
 time = "sec")

}
