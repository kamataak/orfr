#' @title Passage reading data
#' @description A longitudinal data set containing passage reading information of 300 students.
#'
#' @format A data frame with 5304 rows and 7 variables:
#' \describe{
#'   \item{id.student}{unique student identifier}
#'   \item{occasion}{identifier for longitudinal assessment occasions;
#'   here a triannual assessment administered in the fall, winter, and spring of a school year}
#'   \item{grade}{student grade level}
#'   \item{id.passage}{unique passage identifier}
#'   \item{numwords.pass}{total number of words in the passage}
#'   \item{wrc}{words read correct}
#'   \item{sec}{seconds to read the passage}
#' }
#' @source \url{https://jnese.github.io/core-blog/}
"passage"

#' @title Passage-level student data set
#'
#' @description A data set consisted of reading accuracy and time data for 12 passages from 85 students.
#'
#' @format 847 rows and 7 variables:
#' \describe{
#'   \item{id.student}{unique student identifier}
#'   \item{occasion}{identifier for longitudinal assessment occasions;
#'   here a triannual assessment administered in the fall, winter, and spring of a school year}
#'   \item{grade}{student grade level}
#'   \item{id.passage}{unique passage identifier}
#'   \item{numwords.pass}{total number of words in the passage}
#'   \item{wrc}{words read correct}
#'   \item{sec}{seconds to read the passage}
#' }
#' @source \url{https://jnese.github.io/core-blog/}
"passage2"
