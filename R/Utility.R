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
#' prep function prepares input data for fit.model function
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
prep <- function(data=data,person.id="",task.id="",occasion="",group="",max.counts="",obs.counts="",time="") {
  # loading logger
  log.initiating()
  flog.info("Begin preparing data process", name = "orfrlog")
  #  col.names = c(studentid,passageid,numwords.p,season,grade,wrc,time)
  col.labels <- c("person.id","task.id","max.counts","occasion","group","obs.counts","time","lgsec")

  dat <- data
  tryCatch (
    expr = {
      # create_data
      c1 <- dat[person.id] # person.id
      c2 <- dat[task.id] # task.id
      c3 <- dat[max.counts] # max.counts
      c4 <- dat[occasion] # occasion
      c5 <- dat[group] # group
      c6 <- dat[obs.counts] # obs.counts
      c7 <- dat[time] # time
      lgsec <- log(c7) # lgsec

      dat <- data.frame(c1,c2,c3,c4,c5,c6,c7,lgsec)
      colnames(dat) <- col.labels

      tp <- as.data.frame(dat %>% select(person.id, task.id, obs.counts) %>%
                            pivot_wider(names_from = task.id, values_from = obs.counts))

      rownames(tp) <- as.character(tp$person.id)
      Y <- tp %>% select(-person.id)
      Y <- Y[ , order(names(Y))] # sort by passage.id
      Y <- as.matrix(Y)
      for (i in 1:ncol(Y)) {
        Y[,i]<-ifelse(is.na(Y[,i]),NA,Y[,i]) #NaN
      }
      logT <-  as.data.frame(dat %>%
                               mutate(lgsec=log(time)) %>%
                               select(person.id, task.id, lgsec) %>%
                               pivot_wider(names_from = task.id, values_from = lgsec) %>%
                               select(-person.id))
      rownames(logT) <- as.character(tp$person.id)
      logT <- logT[ , order(names(logT))] # sort by person.id
      N <- as.data.frame(dat %>%
                           group_by(task.id) %>% arrange(task.id) %>% # sort by person.id
                           summarise(max.counts=max(max.counts)) %>% # numwords.pass
                           select(-task.id))
      rownames(N) <- colnames(Y)
      N <- pull(N)
      I <- length(N)
      N.matrix <- matrix(rep(as.matrix(N),dim(Y)[1]),nrow = dim(Y)[1], byrow = TRUE)
      logT10 <- logT - log(N.matrix) + log(10)
      logT10 <- logT10[ , order(names(logT10))]
      # data.in <- list(Y = Y, logT10 = logT10, N = N, I = I)
      data.in <- list(Y = Y, logT10 = logT10, N = N, I = I)


      output <- list(data.long=dat,
                     data.wide=data.in)
      flog.info("End preparing data process", name = "orfrlog")

      return(output)
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
}
#' Returns cases (person and occasion) applied in [fit.model] function.
#'
#' @param data = person response data
#'
#' @return cases vector
#'
#' @export
get.cases <- function(data) {
  cases <- data %>% select(person.id,occasion) %>% unique() %>%
    unite("cases", person.id:occasion, sep = "_", remove = TRUE, na.rm = FALSE) %>%
    select(cases)
  print(cases)
  return(invisible(cases))
}
#' Returns perfect cases (student and occasion) in which every word was read correctly.
#'
#' @param data  = person response data
#'
#' @return perfect accurate case vector
#' @export
#'
get.perfectcases <- function(data) {
  perfect.cases <- data %>% group_by(person.id,occasion) %>%
    summarise(obs.counts.sum=sum(obs.counts),
              max.counts.sum=sum(max.counts), .groups = "drop_last") %>%
    filter(obs.counts.sum == max.counts.sum) %>%
    unite("perfect.cases", person.id:occasion, sep = "_", remove = TRUE, na.rm = FALSE) %>%
    select(perfect.cases)
  return(invisible(perfect.cases))
}

#' Prepares data in a long format for [wcpm]
#'
#' @param data = person response data
#'
#' @return data frame
#'
preplong <- function(data,
                     person.id,
                     task.id,
                     occasion,
                     group,
                     max.counts,
                     obs.counts,
                     time){
  vars <- c(person.id,
            task.id,
            occasion,
            group,
            max.counts,
            obs.counts,
            time)
  dat <- data %>%
    select(all_of(vars)) %>%
    rename(person.id=1,task.id=2,
           occasion=3,group=4,
           max.counts=5,obs.counts=6,time=7) %>%
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
#'
prepwide <- function(data,
                     person.id,
                     task.id,
                     max.counts,
                     obs.counts,
                     time){
  vars <- c(person.id,
            task.id,
            max.counts,
            obs.counts,
            time)
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
#' To exclude error passages
#'
#' @param passage
#'
#' @return passage data set without error passages
#' @export
exclude_passages <- function(passage) {
  err_list <- get_errlist(passage)
  return (passage %>% filter(!(id.passage %in% err_list)))

}

#' To get a string of error passages that have no at least two readers
#'
#' @param passage
#'
#' @return a string of error passages
get_errlist <- function(passage) {
  # get unique passage list
  passage_ids <- as.matrix(passage %>% select(id.passage) %>% unique())

  flag <- 0
  err_list <- c()
  for (i in 1:length(passage_ids)) {
    flag <- 0
    #print(passage_ids[i])
    set_a <- passage %>% filter(id.passage==passage_ids[i]) %>% select(id.student)
    for (j in 1:length(passage_ids)) {
      if (j != i) {
        set_b <- passage %>% filter(id.passage==passage_ids[j]) %>% select(id.student)
        if (nrow(intersect(set_a,set_b)) > 2) {
          flag <- 1
          #print("break")
          break
        }
      }
    }
    if (flag == 0) {
      print(paste("got ",passage_ids[i]))
      ll <- length(err_list)+1
      err_list[ll] <- passage_ids[i]
    }

  }
  return (err_list)
}
