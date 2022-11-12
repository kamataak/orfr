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
#' @param person.data - individual reading data
#' @param person.id The column name in the data that represents the unique personal identifier.
#' @param task.id The column name in the data that represents the unique passage identifier.
#' @param max.counts The column name in the data that represents the number of words in a passage.
#' @param obs.counts The column name in the data that represents the words read correctly for each case.
#' @param time The column name in the data that represents the time, in seconds, for each case.
#' @param k.in    - number of imputations, default is 5
#' @param reps.in - number of Monte-Carlo iterations, default is 2
#' @param ests.in - if not given, mom function will be called and get est.in output
#' @param est - estimator keyword, mcem or mcmc, default is mcem
#' @param se - standard error keyword / c("none","analytic", "bootstrap"), default is none
#' @param verbose - boolean, if shows the summary, default is FALSE
#'
#' @return MCEM list, MCMC list
#' @export
fit.model <- function(data=NA, person.data=NA, person.id="",task.id="",max.counts="",obs.counts="",time="", k.in=5,reps.in=2,ests.in=NA,
                      est="mcem",se="none",verbose=FALSE) {
  # loading logger
  log.initiating()

  if (person.id != "") {
    #create wide data
    data <- prepwide(person.data,person.id,task.id,max.counts,obs.counts,time)
  }

  if (est == "mcem") {
    if (se == "none") {

      MCEMests <- run.mcem(data$Y,data$logT10,data$N,data$I,k.in,reps.in,ests.in,verbose=verbose)

      task.param <-  tibble(
        a = MCEMests$a,
        b = MCEMests$b,
        alpha = MCEMests$alpha,
        beta = MCEMests$beta,
        task.id = as.numeric(colnames(data$Y)),
        max.counts = data$N)
      hyper.param <- tibble(vartau = MCEMests$vartau,
                            rho = MCEMests$rho)
      MCEM.ests <- list(task.param = task.param,
                        hyper.param = hyper.param)
    } else {
      # test
      # if (is.na(ests.in)) {
      #   ests.in <- mom(data$Y, data$logT10, data$N, data$I)
      # }
      flog.info("Begin mcem process", name = "orfrlog")
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

      task.param <-  tibble(
        a = MCEMests$a,
        b = MCEMests$b,
        alpha = MCEMests$alpha,
        beta = MCEMests$beta,
        se_a = SE.analytic[1:length(MCEMests$a)],
        se_b = SE.analytic[(length(MCEMests$a)+1):(length(MCEMests$a)+length(MCEMests$b))],
        se_alpha = SE.analytic[(length(MCEMests$a)+length(MCEMests$b)+1):(length(MCEMests$a)+length(MCEMests$b)+length(MCEMests$alpha))],
        se_beta = SE.analytic[(length(MCEMests$a)+length(MCEMests$b)+length(MCEMests$alpha)+1):(length(MCEMests$a)+length(MCEMests$b)+length(MCEMests$alpha)+length(MCEMests$beta))],
        task.id = as.numeric(colnames(data$Y)),
        max.counts = data$N)
      hyper.param <- tibble(vartau = MCEMests$vartau,
                            rho = MCEMests$rho,
                            se_vartau = SE.analytic[(length(MCEMests$a)+length(MCEMests$b)+length(MCEMests$alpha)+length(MCEMests$beta)+1)],
                            se_rho = SE.analytic[(length(MCEMests$a)+length(MCEMests$b)+length(MCEMests$alpha)+length(MCEMests$beta)+2)])
      MCEM.ests <- list(task.param = task.param,
                        hyper.param = hyper.param)
    }


    # check if shows the summary
    if (verbose == TRUE) {
      summary.fit.model(MCEM.ests)
    }
    class(MCEM.ests) <- "fit.model" # define class
    flog.info("End mcem process", name = "orfrlog")
    return(invisible(MCEM.ests))


  } else { # for MCMC, mcem parameters are necessary
    flog.info("Begin mcmc process", name = "orfrlog")
    MCMC.ests <- bayes(person.data,
                       person.id,
                       task.id,
                       max.counts,
                       obs.counts,
                       time,
                       parallel=T, #logical, run in parallel? "T" or "F"
                       n.chains=NA, # pos. int., number of the chains
                       thin=1, #pos. int, thinning interval, a.k.a, period of saving samples
                       iter=NA,  # pos. int., number of the iterations after the burn-in period
                       burn=NA  # pos. int., number of the burn-in iterations)
    )
    # check if shows the summary
    if (verbose == TRUE) {
      summary.fit.model(MCMC.ests)
    }

    class(MCMC.ests) <- "fit.model" # define class
    flog.info("End mcmc process", name = "orfrlog")
    return(invisible(MCMC.ests))
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
#' @param calib.data - fit.model class object
#' @param person.data - individual reading data
#' @param person.id The column name in the data that represents the unique individual identifier.
#' @param task.id The column name in the data that represents the unique task identifier.
#' @param max.counts The column name in the data that represents the number of words in a task.
#' @param occasion The column name in the data that represents the unique occasion.
#' @param group The column name in the data that represents the unique group.
#' @param obs.counts The column name in the data that represents the words read correctly for each case.
#' @param time The column name in the data that represents the time, in seconds, for each case.
#' @param cases - individual id vectors, will directly use task data if no calib.data provided
#' @param est - estimator keyword / c("mle", "map", "eap", "mcmc"), default is mcmc
#' @param se - standard error keyword / c("analytic", "bootstrap"), default is analytic
#' @param failsafe - retry time for bootstrap / default 0, can set to 5 ~ 50
#' @param bootstrp - set K number of bootstrap / default 100
#' @param external - if not NULL, will use not student read passages for estimating
#' @param type - output type, "general" and "orf", default "general" only output tau & theta. "orf" will output wcpm
#'
#' @return scoring list or Bootstrap dataset
#' @export
scoring <- function(calib.data=NA, person.data=NA, person.id="",task.id="",occasion="",group="",max.counts="",obs.counts="",time="", cases=NULL,
                    est="map", se="analytic",failsafe=0, bootstrap=100, external=NULL, type="general") {
  # loading logger
  log.initiating()

  bootstrap.out <- tibble()
  error_case <- tibble()
  if (se == "analytic") {
    if (est != "mcmc") { #not mcmc
      if (person.id != "") {
        person.data <- preplong(person.data,person.id,task.id,occasion,group,max.counts,obs.counts,time)
      }
      # Check MCEM object
      if (class(calib.data)[1] == "fit.model") {
        #      MCEM <- calib.data
        # assign task.data
        task.data <- calib.data$task.param
      } else { # if no MCEM object stop running
        flog.info("Missed fit.model object, end scoring process", name = "orfrlog")
        return(NULL)
      }

      if (length(external) != 0) { # external,
        print(paste("Use external task:", paste(external, collapse = ",")))
      }

      # Check cases
      if (length(cases) == 0) {
        print("Cases: ")
        cases <- get.cases(person.data)
      }
      # Check if there is a perfect accurate case
      perfect.cases <<- get.perfectcases(person.data)

      if (count(perfect.cases) != 0) {
        flog.info(paste("The perfect accurate case: ", paste(perfect.cases$perfect.cases, collapse = ", ")), name="orfrlog")
      } else {
        flog.info("There is no perfect accurate case.", name="orfrlog")
      }

      result.list <- run.scoring(calib.data, person.data, task.data, cases, perfect.cases, est, lo = -4, hi = 4, q = 100, kappa = 1, external=external,type=type)
    } else  { # mcmc
      if (person.id == "") { # if without columns' names
        result.list <- bayes.wcpm(
          calib.data = calib.data,
          person.data = person.data,
          cases = cases,
          external=external,
          parallel=T, #logical, run in parallel? "T" or "F"
          n.chains=NA, # pos. int., number of the chains
          iter=NA,  # pos. int., number of the iterations after the burn-in period
          burn=NA,  # pos. int., number of the burn-in iterations
          thin=1 #pos. int, thinning interval, a.k.a, period of saving samples
        )
      } else {
        result.list <- bayes.wcpm(
          calib.data = calib.data,
          person.data = person.data,
          person.id = person.id,
          task.id = task.id,
          occasion = occasion,
          group = group,
          max.counts = max.counts,
          obs.counts = obs.counts,
          time = time,
          cases = cases,
          external=external,
          parallel=T, #logical, run in parallel? "T" or "F"
          n.chains=NA, # pos. int., number of the chains
          iter=NA,  # pos. int., number of the iterations after the burn-in period
          burn=NA,  # pos. int., number of the burn-in iterations
          thin=1 #pos. int, thinning interval, a.k.a, period of saving samples
        )
      }
      class(result.list) <- "scoring"
      return(result.list)
    }
  } else if (se == "bootstrap"){ #for bootstrap
    # Check if there is a perfect accurate case
    perfect.cases <<- get.perfectcases(person.data)

    RE_TRY <- failsafe # Define retry, if 0, no retry
    j <- 0 # index for retry time
    i <- 1 # index for case loop

    t_size <- nrow(cases)

    while (i <= t_size) {
      temp <- tibble()
      flog.info(paste("Boostrap running for case:", cases$cases[i]), name = "orfrlog")
      t_case = data.frame(cases=cases$cases[i])
      tryCatchLog(
        temp <- getBootstrapSE(calib.data, person.data, case=t_case, perfect.cases, est, kappa=1,bootstrap=bootstrap),
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
