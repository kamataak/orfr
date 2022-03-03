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
#' @return data list (Y, logT10, N, I)
#'
#' @export
prep <- function(data, vars="") {
  # loading logger
  log.initiating()
  flog.info("Begin preparing data process", name = "orfrlog")
  check_set <- c("id.student","occasion","id.passage","grade","numwords.pass","wrc","sec")
  anydiff <- setdiff(vars,check_set)
  if (length(anydiff) > 0) {
    flog.info(paste("Variable incorrect:", anydiff), name = "orfrlog")
    return
  }

  dat <- data

  if (length(vars) > 0) {
    # create_data
    for (i in 1:length(vars)) {
      #      print(i)
      if (vars[i] == "id.student") {
        c1 <- unname(dat[i])
      } else if (vars[i] == 'id.passage') {
        c2 <- unname(dat[i])
      } else if (vars[i] == 'numwords.pass') {
        c3 <- unname(dat[i])
      } else if (vars[i] == 'occasion') {
        c4 <- unname(dat[i])
      } else if (vars[i] == 'grade') {
        c5 <- unname(dat[i])
      } else if (vars[i] == 'sec') {
        c6 <- unname(dat[i])
        lgsec <- log(c6)
      } else if (vars[i] == 'wrc') {
        c7 <- unname(dat[i])
      }
    }

    dat <- data.frame(student.id=c1,
                      passage.id=c2,
                      nwords.p=c3,
                      occasion=c4,
                      grade=c5,
                      sec=c6,
                      wrc=c7,
                      lgsec)

  }

  tryCatch(
    expr = {
      Y <- dat %>% select(student.id, passage.id, wrc) %>%
        pivot_wider(names_from = passage.id, values_from = wrc) %>%
        select(-student.id)
      Y <- Y[ , order(names(Y))] # sort by passage.id
      Y <- as.matrix(Y)
      for (i in 1:ncol(Y)) {
        Y[,i]<-ifelse(is.na(Y[,i]),NaN,Y[,i])
      }
      logT <- dat %>%
        mutate(lgsec=log(sec)) %>%
        select(student.id, passage.id, lgsec) %>%
        pivot_wider(names_from = passage.id, values_from = lgsec) %>%
        select(-student.id)
      logT <- logT[ , order(names(logT))] # sort by passage.id
      N <- dat %>%
        group_by(passage.id) %>% arrange(passage.id) %>% # sort by passage.id
        summarise(numwords.pass=max(nwords.p)) %>%
        select(-passage.id)
      N <- pull(N)
      I <- length(N)
      N.matrix <- matrix(rep(as.matrix(N),dim(Y)[1]),nrow = dim(Y)[1], byrow = TRUE)
      logT10 <- logT - log(N.matrix) + log(10)
      logT10 <- logT10[ , order(names(logT10))]
      # data.in <- list(Y = Y, logT10 = logT10, logT = logT, N = N, I = I)
      data.in <- list(Y = Y, logT10 = logT10, N = N, I = I)
    },
    warning = function(w) {
      flog.info("There was a warning message. Something is wrong!", name = "orfrlog")
      flog.info(w, name = "orfrlog")
    },
    error = function(w) {
      flog.info("There was an error message. Something is wrong!", name = "orfrlog")
      flog.info(w, name = "orfrlog")
    }
  )


  output <- list(data.raw=dat,
                 data.in=data.in)
  flog.info("End preparing data process", name = "orfrlog")

  return(output)

}
#' The function get return cases used for wcpm function
#'
#' @param data = student response data
#'
#' @return cases vector
#'
#' @export
get.cases <- function(data) {
  cases <- data %>% select(student.id,occasion) %>% unique() %>%
    unite("cases", student.id:occasion, sep = "_", remove = TRUE, na.rm = FALSE) %>%
    select(cases)
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
  perfect.cases <- data %>% group_by(student.id,occasion) %>%
    summarise(wrc_sum=sum(wrc),
              nwords.p_sum=sum(nwords.p), .groups = "drop_last") %>%
    filter(wrc_sum == nwords.p_sum) %>%
    unite("perfect.cases", student.id:occasion, sep = "_", remove = TRUE, na.rm = FALSE) %>%
    select(perfect.cases)
  return(invisible(perfect.cases))
}
