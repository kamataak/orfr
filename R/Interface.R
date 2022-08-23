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
#' @param data - reading data object comes from prep function
#' @param stu.data - student reading data
#' @param studentid The column name in the data that represents the unique student identifier.
#' @param passageid The column name in the data that represents the unique passage identifier.
#' @param numwords.p The column name in the data that represents the number of words in a passage.
#' @param wrc The column name in the data that represents the words read correctly for each case.
#' @param time The column name in the data that represents the time, in seconds, for each case.
#' @param k.in    - number of imputations, default is 5
#' @param reps.in - number of Monte-Carlo iterations, default is 2
#' @param ests.in - if not given, mom function will be called and get est.in output
#' @param est - estimator keyword, mcem or mcmc
#' @param se - standard error keyword / c("none","analytic", "bootstrap"), default is none
#' @param verbose - boolean, if shows the summary, default is FALSE
#'
#' @return MCEM list, MCMC list
#' @export
mcem <- function(data=NA, stu.data=NA, studentid="",passageid="",numwords.p="",wrc="",time="", k.in=5,reps.in=2,ests.in=NA,
                 est="mcem",se="none",verbose=FALSE) {
  # loading logger
  log.initiating()

  if (studentid != "") {
    #create wide data
    data <- prepwide(stu.data,studentid,passageid,numwords.p,wrc,time)
  }

  if (est == "mcem") {

    if (se == "none") {
      flog.info("Begin mcem process without se", name = "orfrlog")
      MCEMests <- run.mcem(data$Y,data$logT10,data$N,data$I, k.in=k.in,reps.in=reps.in,ests.in=ests.in,verbose=verbose)

      pass.param <-  tibble(
        a = MCEMests$a,
        b = MCEMests$b,
        alpha = MCEMests$alpha,
        beta = MCEMests$beta,
        passage.id = as.numeric(colnames(data$Y)),
        numwords.p = data$N)
      hyper.param <- tibble(vartau = MCEMests$vartau,
                            rho = MCEMests$rho)
      MCEM.ests <- list(pass.param = pass.param,
                        hyper.param = hyper.param)
    } else {
      # test
      # if (is.na(ests.in)) {
      #   ests.in <- mom(data$Y, data$logT10, data$N, data$I)
      # }
      flog.info(paste("Begin mcem process with se",se), name = "orfrlog")

      MCEMests <- run.mcem(data$Y, data$logT10, data$N, data$I,
                           k.in=k.in, reps.in=reps.in,ests.in=ests.in)

      MCEMout <- c(MCEMests$a, MCEMests$b, MCEMests$alpha,
                   MCEMests$beta, MCEMests$vartau, MCEMests$rho)
      if (se == "analytic") {
        CV.analytic <- numerical.cov(data$Y, data$logT10, data$N, data$I,
                                     MCEMout,h.val=1e-10, M=100)
      } else { # bootstrap
        # CV.analytic <- boot.cov(data$Y, data$logT10, data$N, data$I,
        #                         k.in=8, reps.in=3, B=10)
        # an alternative function
        CV.analytic <- bootmodel.cov(data$Y, data$logT10, data$N, data$I,
                                     MCEMout,k.in=k.in, reps.in=reps.in,B=10)

      }
      SE.analytic <- sqrt(diag(CV.analytic))

      pass.param <-  tibble(
        a = MCEMests$a,
        b = MCEMests$b,
        alpha = MCEMests$alpha,
        beta = MCEMests$beta,
        se_a = SE.analytic[1:length(MCEMests$a)],
        se_b = SE.analytic[(length(MCEMests$a)+1):(length(MCEMests$a)+length(MCEMests$b))],
        se_alpha = SE.analytic[(length(MCEMests$a)+length(MCEMests$b)+1):(length(MCEMests$a)+length(MCEMests$b)+length(MCEMests$alpha))],
        se_beta = SE.analytic[(length(MCEMests$a)+length(MCEMests$b)+length(MCEMests$alpha)+1):(length(MCEMests$a)+length(MCEMests$b)+length(MCEMests$alpha)+length(MCEMests$beta))],
        passage.id = as.numeric(colnames(data$Y)),
        numwords.p = data$N)
      hyper.param <- tibble(vartau = MCEMests$vartau,
                            rho = MCEMests$rho,
                            se_vartau = SE.analytic[(length(MCEMests$a)+length(MCEMests$b)+length(MCEMests$alpha)+length(MCEMests$beta)+1)],
                            se_rho = SE.analytic[(length(MCEMests$a)+length(MCEMests$b)+length(MCEMests$alpha)+length(MCEMests$beta)+2)])
      MCEM.ests <- list(pass.param = pass.param,
                        hyper.param = hyper.param)
    }


    # check if shows the summary
    if (verbose == TRUE) {
      summary.mcem(MCEM.ests)
    }
    class(MCEM.ests) <- "mcem" # define class
    flog.info("End mcem process", name = "orfrlog")
    return(invisible(MCEM.ests))


  } else { # for MCMC, mcem parameters are necessary
    # Check MCEM object
    if (class(data)[1] == "mcem") {
      flog.info("Using mcem parameter for MCMC", name = "orfrlog")

      MCEM <- data
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
#' @param calib.data - mcem class object
#' @param stu.data - student reading data
#' @param studentid The column name in the data that represents the unique student identifier.
#' @param passageid The column name in the data that represents the unique passage identifier.
#' @param numwords.p The column name in the data that represents the number of words in a passage.
#' @param wrc The column name in the data that represents the words read correctly for each case.
#' @param time The column name in the data that represents the time, in seconds, for each case.
#' @param cases - student id vectors, will directly use passage data if no calib.data provided
#' @param est - estimator keyword / c("mle", "map", "eap")
#' @param se - standard error keyword / c("analytic", "bootstrap"), default is analytic
#' @param wo - wcpm option / c("internal", "external"), default is internal
#' @param failsafe - retry time for bootstrap / default 0, can set to 5 ~ 50
#' @param bootstrp - set K number of bootstrap / default 100
#' @param external - if not NULL, will use not student read passages for estimating
#'
#' @return WCPM list or Bootstrap dataset
#' @export
wcpm <- function(calib.data=NA, stu.data=NA, studentid="",passageid="",season="",grade="",numwords.p="",wrc="",time="", cases=NULL,
                 est="map", se="analytic", wo="internal", failsafe=0, bootstrap=100, external=NULL) {
  # loading logger
  log.initiating()

  # if (class(calib.data)[1] != "mcem" ) {
  #   #call mcem
  #   calib.data <- mcem(stu.data=stu.data,studentid=studentid,passageid=passageid,
  #                      numwords.p=numwords.p,wrc=wrc,time=time,est="mcem")
  #   #create long data
  #   stu.data <- preplong(stu.data,studentid,passageid,season,grade,numwords.p,wrc,time)
  #   #    pass.data <- MCEM$pass.param
  #   #    calib.data <- MCEM
  # }
  if (studentid != "") {
    stu.data <- preplong(stu.data,studentid,passageid,season,grade,numwords.p,wrc,time)
  }
  # Check MCEM object
  if (class(calib.data)[1] == "mcem") {
    #      MCEM <- calib.data
    # assign pass.data
    pass.data <- calib.data$pass.param
  } else { # if no MCEM object stop running
    flog.info("Missed MCEM object, end wcpm process", name = "orfrlog")
    return(NULL)
  }

  if (length(external) != 0) { # external,
    print(paste("Use external passage:", paste(external, collapse = ",")))
  }

  # Check cases
  if (length(cases) == 0) {
    print("Cases: ")
    cases <- get.cases(stu.data)
  }
  # Check if there is a perfect accurate case
  perfect.cases <- get.perfectcases(stu.data)

  if (count(perfect.cases) != 0) {
    flog.info(paste("The perfect accurate case: ", paste(perfect.cases$perfect.cases, collapse = ", ")), name="orfrlog")
  } else {
    flog.info("There is no perfect accurate case.", name="orfrlog")
  }

  bootstrap.out <- tibble()
  error_case <- tibble()
  if (se == "analytic") {
    run.wcpm(calib.data, stu.data, pass.data, cases, perfect.cases, est, lo = -4, hi = 4, q = 100, kappa = 1, external=external)
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
        temp <- getBootstrapSE(calib.data, stu.data, case=t_case, perfect.cases, est, kappa=1,bootstrap=bootstrap),
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
