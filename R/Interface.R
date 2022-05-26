#' This is an interface function to call and run mcem.
#'
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
#' @param object - reading data object comes from prep function
#' @param k.in - number of passages, default is 5
#' @param reps.in - repeats, default is 2
#' @param ests.in - if not given, mom function will be called and get est.in output
#' @param est - estimator keyword, mcem or mcmc
#' @param se - standard error keyword, default is analytic
#' @param verbose - boolean, if shows the summary, default is FALSE
#'
#' @return MCEM list, MCMC list
#' @export
mcem <- function(object, studentid="",passageid="",numwords.p="",wrc="",time="", k.in=5,reps.in=2,ests.in,
                 est="mcem",se="analytic",verbose=FALSE) {
  # loading logger
  log.initiating()

  if (studentid != "") {
    #create wide data
    object <- prepwide(object,studentid,passageid,numwords.p,wrc,time)
  }

  if (est == "mcem") {
    dat <- object
    return(
      run.mcem(dat$Y,dat$logT10,dat$N,dat$I,k.in,reps.in,ests.in,verbose=verbose)
    )
  } else { # for MCMC, mcem parameters are necessary
    # Check MCEM object
    if (class(object)[1] == "mcem") {
      flog.info("Using mcem parameter for MCMC", name = "orfrlog")

      MCEM <- object
      return(
        #        runBayes(mcem=MCEM,cases)
      )
    } else { #run the whole process
      flog.info("Wrong parameters.", name = "orfrlog")
      #      runBayes(object=dat,mcem=MCEM,wcpm=WCPM,cases)
    }
  }
}

#' This is an interface function to call and run wcpm or bootstrap.
#'
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
#' @param object - mcem class object
#' @param stu.data - student reading data
#' @param cases - student id vectors
#' @param est - estimator keyword / c("mle", "map", "eap")
#' @param se - standard error keyword / c("analytic", "bootstrap")
#' @param wo - wcpm option / c("internal", "external"), default is internal
#' @param failsafe - retry time for bootstrap / default 0, can set to 5 ~ 50
#' @param bootstrp - set K number of bootstrap / default 100
#' @param external - if not NULL, will use not student read passages for estimating
#'
#' @return WCPM list or Bootstrap dataset
#' @export
wcpm <- function(object, studentid="",passageid="",season="",grade="",numwords.p="",wrc="",time="", stu.data=data, cases=NA,
                 est="map", se="analytic", wo="internal", failsafe=0, bootstrap=100, external=NULL) {
  # loading logger
  log.initiating()
  # if (studentid != "" & class(object)[1] != "mcem" ) {

  if (class(object)[1] != "mcem" ) {
    #call mcem
    MCEM <- mcem(object,studentid,passageid,numwords.p,wrc,time,est="mcem")
    #create long data
    stu.data <- preplong(object,studentid,passageid,season,grade,numwords.p,wrc,time)
    #    pass.data <- MCEM$pass.param
    object <- MCEM
  }

  # Check MCEM object
  if (wo=="internal") { # internal, object must be mcem object
    if (class(object)[1] == "mcem") {
      MCEM <- object
      # assign pass.data
      pass.data <- MCEM$pass.param
    } else { # if no MCEM object stop running
      flog.info("Missed MCEM object, end wcpm process", name = "orfrlog")
      return(NULL)
    }
  } else { # external,
    #check and prepare data
    #1. check if list of parameter includes a,b,alpha,beta
    #2. check response data
    #3. prepare data for running
    print("Use user supplied list of parameters and response data")
    return(NULL)
  }

  # Check if there is a perfect accurate case
  perfect.cases <- get.perfectcases(stu.data)

  if (count(perfect.cases) != 0) {
    flog.info(paste("The perfect accurate case: ", perfect.cases$perfect.cases), name="orfrlog")
  } else {
    flog.info("There is no perfect accurate case.", name="orfrlog")
  }

  bootstrap.out <- tibble()
  error_case <- tibble()
  if (se == "analytic") {
    run.wcpm(object, stu.data, pass.data, cases, perfect.cases, est, lo = -4, hi = 4, q = 100, kappa = 1, external=external)
  } else if (se == "bootstrap"){ #for bootstrap

    RE_TRY <- failsafe # Define retry, if 0, no retry
    j <- 0 # index for retry time
    i <- 1 # index for case loop

    t_size <- nrow(cases)

    while (i <= t_size) {
      temp <- tibble()
      flog.info(paste("Boostrap running for case:", cases$cases[i]), name = "orfrlog")
      t_case = data.frame(cases=cases$cases[i])
      tryCatchLog(
        temp <- getBootstrapSE(object, stu.data, case=t_case, perfect.cases, est, kappa=1,bootstrap=bootstrap),
        error=function(e) {
          flog.info(paste("Running error:", e), name = "orfrlog")
        }
      )

      if (length(temp) > 2) { #without error
        bootstrap.out <- rbind(bootstrap.out,temp)
        j <- 0 # reset index of retry
        i <- i + 1 # go to next case
      } else { # with error
        if (j == RE_TRY) { # after RE_TRY if still error
          error_temp <- tibble(case_id = cases$cases[i])
          error_case <- rbind(error_case, error_temp)
          i <- i + 1 # go to next case
        } else {
          # let it retry a couple of times
          j <- j + 1
          flog.info(paste("Boostrap try again: time", j), name = "orfrlog")
        }
      }
    }
    result.list <- list(bootstrap.out = bootstrap.out,
                        error_case = error_case)
    class(result.list) <- "bootstrap"
    return(result.list)
  }
}
#' This is an interface function to call and run mcem.
#'
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
#' @param object - reading data object comes from prep function
#' @param k.in - number of passages, default is 5
#' @param reps.in - repeats, default is 2
#' @param ests.in - if not given, mom function will be called and get est.in output
#' @param est - estimator keyword, mcem or mcmc
#' @param se - standard error keyword, default is analytic
#' @param verbose - boolean, if shows the summary, default is FALSE
#'
#' @return MCEM list, MCMC list
#' @export
mcem <- function(object, studentid="",passageid="",numwords.p="",wrc="",time="", k.in=5,reps.in=2,ests.in,
                 est="mcem",se="analytic",verbose=FALSE) {
  # loading logger
  log.initiating()

  if (studentid != "") {
    #create wide data
    object <- prepwide(object,studentid,passageid,numwords.p,wrc,time)
  }

  if (est == "mcem") {
    dat <- object
    return(
      run.mcem(dat$Y,dat$logT10,dat$N,dat$I,k.in,reps.in,ests.in,verbose=verbose)
    )
  } else { # for MCMC, mcem parameters are necessary
    # Check MCEM object
    if (class(object)[1] == "mcem") {
      flog.info("Using mcem parameter for MCMC", name = "orfrlog")

      MCEM <- object
      return(
        #        runBayes(mcem=MCEM,cases)
      )
    } else { #run the whole process
      flog.info("Wrong parameters.", name = "orfrlog")
      #      runBayes(object=dat,mcem=MCEM,wcpm=WCPM,cases)
    }
  }
}

#' This is an interface function to call and run wcpm or bootstrap.
#'
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
#' @param object - mcem class object
#' @param stu.data - student reading data
#' @param cases - student id vectors
#' @param est - estimator keyword / c("mle", "map", "eap")
#' @param se - standard error keyword / c("analytic", "bootstrap")
#' @param wo - wcpm option / c("internal", "external"), default is internal
#' @param failsafe - retry time for bootstrap / default 0, can set to 5 ~ 50
#' @param bootstrp - set K number of bootstrap / default 100
#' @param external - if not NULL, will use not student read passages for estimating
#'
#' @return WCPM list or Bootstrap dataset
#' @export
wcpm <- function(object, studentid="",passageid="",season="",grade="",numwords.p="",wrc="",time="", stu.data=data, cases=NA,
                 est="map", se="analytic", wo="internal", failsafe=0, bootstrap=100, external=NULL) {
  # loading logger
  log.initiating()
  # if (studentid != "" & class(object)[1] != "mcem" ) {

  if (class(object)[1] != "mcem" ) {
    #call mcem
    MCEM <- mcem(object,studentid,passageid,numwords.p,wrc,time,est="mcem")
    #create long data
    stu.data <- preplong(object,studentid,passageid,season,grade,numwords.p,wrc,time)
    #    pass.data <- MCEM$pass.param
    object <- MCEM
  }

  # Check MCEM object
  if (wo=="internal") { # internal, object must be mcem object
    if (class(object)[1] == "mcem") {
      MCEM <- object
      # assign pass.data
      pass.data <- MCEM$pass.param
    } else { # if no MCEM object stop running
      flog.info("Missed MCEM object, end wcpm process", name = "orfrlog")
      return(NULL)
    }
  } else { # external,
    #check and prepare data
    #1. check if list of parameter includes a,b,alpha,beta
    #2. check response data
    #3. prepare data for running
    print("Use user supplied list of parameters and response data")
    return(NULL)
  }

  # Check if there is a perfect accurate case
  perfect.cases <- get.perfectcases(stu.data)

  if (count(perfect.cases) != 0) {
    flog.info(paste("The perfect accurate case: ", perfect.cases$perfect.cases), name="orfrlog")
  } else {
    flog.info("There is no perfect accurate case.", name="orfrlog")
  }

  bootstrap.out <- tibble()
  error_case <- tibble()
  if (se == "analytic") {
    run.wcpm(object, stu.data, pass.data, cases, perfect.cases, est, lo = -4, hi = 4, q = 100, kappa = 1, external=external)
  } else if (se == "bootstrap"){ #for bootstrap

    RE_TRY <- failsafe # Define retry, if 0, no retry
    j <- 0 # index for retry time
    i <- 1 # index for case loop

    t_size <- nrow(cases)

    while (i <= t_size) {
      temp <- tibble()
      flog.info(paste("Boostrap running for case:", cases$cases[i]), name = "orfrlog")
      t_case = data.frame(cases=cases$cases[i])
      tryCatchLog(
        temp <- getBootstrapSE(object, stu.data, case=t_case, perfect.cases, est, kappa=1,bootstrap=bootstrap, external=external),
        error=function(e) {
          flog.info(paste("Running error:", e), name = "orfrlog")
        }
      )

      if (length(temp) > 2) { #without error
        bootstrap.out <- rbind(bootstrap.out,temp)
        j <- 0 # reset index of retry
        i <- i + 1 # go to next case
      } else { # with error
        if (j == RE_TRY) { # after RE_TRY if still error
          error_temp <- tibble(case_id = cases$cases[i])
          error_case <- rbind(error_case, error_temp)
          i <- i + 1 # go to next case
        } else {
          # let it retry a couple of times
          j <- j + 1
          flog.info(paste("Boostrap try again: time", j), name = "orfrlog")
        }
      }
    }
    result.list <- list(bootstrap.out = bootstrap.out,
                        error_case = error_case)
    class(result.list) <- "bootstrap"
    return(result.list)
  }
}
