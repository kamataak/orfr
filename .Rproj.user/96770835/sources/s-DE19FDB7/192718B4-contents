#' This file includes utilities of orfr package.
#'
#' Copyright (C) 2021 The ORF Project Team
#'
#' This program is free software; you can redistribute it and/or modify
#' it under the terms of the GNU General Public License as published by
#' the Free Software Foundation; either version 3 of the License, or
#' (at your option) any later version.
#
#' This program is distributed in the hope that it will be useful,
#' but WITHOUT ANY WARRANTY; without even the implied warranty of
#' MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#' GNU General Public License for more details.
#
#' A copy of the GNU General Public License is available at
#' http://www.gnu.org/licenses/
#'
#'
#' prep function prepares input data for mcem function
#'
#' @param data = student response data
#'
#' @import tidyr
#' @import dplyr
#' @import tidyverse
#'
#' @return data list (Y, logT10, logT, N, I)
#'
#' @export
prep <- function(data){
  # loading logger
  log.initiating()
  flog.info("Begin preparing data process", name = "mcemlog")

  dat <- data
  Y <- dat %>%
    select(student_id, passage_id, wrc) %>%
    spread(key = passage_id, value = wrc) %>%
    select(-student_id)
  Y <- as.matrix(Y)
  for (i in 1:ncol(Y)) {
    Y[,i]<-ifelse(is.na(Y[,i]),NaN,Y[,i])
  }
  logT <- dat %>%
    mutate(logsecs=log(secs)) %>%
    select(student_id, passage_id, logsecs) %>%
    spread(key = passage_id, value = logsecs) %>%
    select(-student_id)
  N <- dat %>%
    group_by(passage_id) %>%
    summarise(numwords.pass=max(nwords.p)) %>%
    select(-passage_id)
  N <- pull(N)
  I <- length(N)
  N.matrix <- matrix(rep(as.matrix(N),dim(Y)[1]),nrow = dim(Y)[1], byrow = TRUE)
  logT10 <- logT - log(N.matrix) + log(10)
  data.in <- list(Y = Y, logT10 = logT10, logT = logT, N = N, I = I)

  flog.info("End preparing data process", name = "mcemlog")
  #class(data.in) <- "prep"
  return(data.in)
}
#' The function get return cases used for wcpm function
#'
#' @param data = student response data
#'
#' @return cases vector
#'
#' @export
get.cases <- function(data) {
  cases <-
    data %>%
    .$stu_season_id2 %>% unique()
  return(invisible(cases))
}
#' The function get perfect accurate cases
#'
#' @param data  = student response data
#'
#' @return perfect accurate case vector
#' @export
#'
get.perfectcases <- function(data) {
  perfect_cases <- data %>% group_by(stu_season_id2) %>%
    summarise(wrc_sum=sum(wrc),
              numwords.pass_sum=sum(nwords.p)) %>%
    filter(wrc_sum == numwords.pass_sum) %>%
    select(stu_season_id2)
  return(invisible(perfect_cases))
}
