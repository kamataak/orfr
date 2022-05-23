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
#' @return data list (data.long: data frame,
#'                    data.wide: list of Y, logT10, N, I)
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
                      numwords.p=c3,
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
        summarise(numwords.pass=max(numwords.p)) %>%
        select(-passage.id)
      N <- pull(N)
      I <- length(N)
      N.matrix <- matrix(rep(as.matrix(N),dim(Y)[1]),nrow = dim(Y)[1], byrow = TRUE)
      logT10 <- logT - log(N.matrix) + log(10)
      logT10 <- logT10[ , order(names(logT10))]
      # data.in <- list(Y = Y, logT10 = logT10, N = N, I = I)
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


  output <- list(data.long=dat,
                 data.wide=data.in)
  flog.info("End preparing data process", name = "orfrlog")

  return(output)

}
#' Returns cases (student and occasion) applied in [wcpm] function.
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
#' Returns perfect cases (student and occasion) in which every word was read correctly.
#'
#' @param data  = student response data
#'
#' @return perfect accurate case vector
#' @export
#'
get.perfectcases <- function(data) {
  perfect.cases <- data %>% group_by(student.id,occasion) %>%
    summarise(wrc_sum=sum(wrc),
              numwords.p_sum=sum(numwords.p), .groups = "drop_last") %>%
    filter(wrc_sum == numwords.p_sum) %>%
    unite("perfect.cases", student.id:occasion, sep = "_", remove = TRUE, na.rm = FALSE) %>%
    select(perfect.cases)
  return(invisible(perfect.cases))
}

#' Prepares data in a long format for [wcpm]
#'
#' @param data = student response data
#'
#' @return data frame
#'
#' @export
preplong <- function(data,studentid,passageid,season,grade,numwords.p,wrc,time){
  vars <- c(studentid,passageid,season,grade,numwords.p,wrc,time)
  dat <- data %>%
    select(all_of(vars)) %>%
    rename(student.id=1,passage.id=2,
           occasion=3,grade=4,
           numwords.p=5,wrc=6,sec=7) %>%
    mutate(lgsec=log(.[[7]]))
  #lgsec10 = log(.[[7]] - log(.[[5]]) + log(10)
  #           stu_season_id2=paste(.[[1]],.[[3]],sep="_"))
  return(dat)
}
#' Prepares data in a wide format for [mcem].
#'
#' This function will return a list with 5 elements:
#' Y: a matrix of words read correctly, where rows represent cases (student and occasion) and columns represent passages
#' logt10: a [tibble::tibble()] of words read correctly
#' N: the number of cases per passage
#' I: the number of passages
#'
#' @param data A data frame.
#' @param studentid The column name in the data that represents the unique student identifier.
#' @param passageid The column name in the data that represents the unique passage identifier.
#' @param numwords.p The column name in the data that represents the number of words in a passage.
#' @param wrc The column name in the data that represents the words read correctly for each case.
#' @param time The column name in the data that represents the time, in seconds, for each case.
#'
#' @example
#' data("passage")
#'
#' prepwide(passage,
#'  studentid = "id.student",
#'  passageid = "id.passage",
#'  numwords.p = "numwords.pass",
#'  wrc = "wrc",
#'  time = "sec")
#'
#' @export
prepwide <- function(data, studentid, passageid, numwords.p, wrc, time){
  vars <- c(studentid,passageid,numwords.p,wrc,time)
  dat <- data %>%
    select(all_of(vars))
  Y <- dat %>%
    select(vars[1], vars[2], vars[4]) %>%
    spread(key = vars[2], value = vars[4]) %>%
    select(-vars[1])
  Y <- as.matrix(Y)
  for (i in 1:ncol(Y)) {
    Y[,i]<-ifelse(is.na(Y[,i]),NaN,Y[,i])
  }
  logT <- dat %>%
    mutate(logsecs=log(.[[5]])) %>%
    select(vars[1], vars[2], logsecs) %>%
    spread(key = vars[2], value = logsecs) %>%
    select(-vars[1])
  N <- dat %>%
    group_by_at(2) %>%
    #    summarise_at(3,max) %>%
    summarise_at(.vars = names(.)[3],max) %>%
    select(-vars[2])
  N <- pull(N)
  I <- length(N)
  N.matrix <- matrix(rep(as.matrix(N),dim(Y)[1]),nrow = dim(Y)[1], byrow = TRUE)
  logT10 <- tibble(logT - log(N.matrix) + log(10))
  data.in <- list(Y = Y, logT10 = logT10, N = N, I = I)
  return(data.in)
}
