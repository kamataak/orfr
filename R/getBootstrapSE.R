#' Get bootstrap SE
#' case is a single stu_season_id
#'
#' Added MAP function 07/14/2021
#' Modified a bug of MLE 07/23/2021
#' Modified EAP 10/28/2021
#'
#' @param object - mcem class object
#' @param person.data - student response data
#' @param case - case number
#' @param est - SE type.(MLE, EAP, and MAP.) default MAP
#' @param perfect.cases - perfect accurate case
#' @param kappa - Default kappa = 1, better be 5
#' @param bootstrap - K number of bootstrap, default is 100
#' @param external - if not NULL, will use not student read passages for estimating
#'
#' @return SE dataset
getBootstrapSE <- function (object, person.data, case=NA, perfect.cases, est="map", kappa=1, bootstrap=100, external=NULL) {
  log.initiating()
  flog.info("Begin getBootstrapSE process", name = "orfrlog")

  # datasim_fixedZ is a modified version of the simulation code
  # that lets you specify the latent values (theta and tau).
  # It generates as many parametric bootstrap samples as
  # needed (parameter reps = how many to simulate).
  datasim.fixedZ <- function(a,b,alpha,beta,var_tau,rho,N,I,Z.latent,reps=100) {
    s12 <- rho*sqrt(var_tau)
    sigma <- matrix(c(1,s12,s12,var_tau), ncol=2)
    Y <- matrix(nrow = reps, ncol = I)
    logT10 <- matrix(nrow = reps, ncol = I)
    # logT <- matrix(nrow = reps, ncol = I)
    p.store <- matrix(nrow = reps, ncol = I)
    for (i in 1:reps) {
      p <- pnorm(a*Z.latent[1]-b)
      Y.in <- mapply(function(n,p) rbinom(1,n,p), N, p)
      logT10.in <- 1/alpha*rnorm(I)+beta-Z.latent[2]
      # logT.in <- logT10.in - log(10) + log(N)
      Y[i,] <- Y.in
      logT10[i,] <- logT10.in
      # logT[i,] <- logT.in
    }
    data <- list(Y = Y,logT10 = logT10)
    # data <- list(Y = Y,logT10 = logT10, logT = logT)

    return(data)
  }

  MCEM <- object
  Estimator <- est
  # Run wcpm function and get ALL estimator
  task.data <- MCEM$task.param
  WCPM <- MCEM %>% run.scoring(person.data, task.data=task.data, cases=case, perfect.cases, est=Estimator, lo = -4, hi = 4, q = 100, kappa = 1, external=external, type="orf")
  # run.scoring <- function(object, person.data, task.data, cases, perfect.cases, est="map", lo = -4, hi = 4, q = 100, kappa = 1, external=NULL, type=NULL) {

  # Extract relevant parameters for given case
  # person.dat01 <- person.data %>% filter(stu_season_id2==case)
  # task.read <- person.dat01 %>% select(passage_id)
  # task.dat01 <- task.data %>% semi_join(task.read, by = "passage_id")
  # task.n <- nrow(task.dat01)
  # max.counts.total <- person.dat01 %>% select(max.counts) %>% c() %>% unlist() %>% sum()
  # group <- person.dat01 %>% select(group) %>% c() %>% unlist %>% unique()

  case_split <- unlist(str_split(case, "_"))
  person.dat01 <- person.data %>% filter(person.data$person.id==case_split[1], person.data$occasion==case_split[2])
  task.read <- person.dat01 %>% select(task.id)
  task.dat01 <- task.data %>% semi_join(task.read, by = "task.id")
  task.n <- nrow(task.dat01)
  max.counts.total <- person.dat01 %>% select(max.counts) %>% c() %>% unlist() %>% sum()
  group <- person.dat01 %>% select(group) %>% c() %>% unlist %>% unique()


  max.counts.task <- person.dat01 %>% select(max.counts) %>% c() %>% unlist()

  a.par <- task.dat01 %>% select(a) %>% c() %>% unlist()
  b.par <- task.dat01 %>% select(b) %>% c() %>% unlist()
  alpha.par <- task.dat01 %>% select(alpha) %>% c() %>% unlist()
  beta.par <- task.dat01 %>% select(beta) %>% c() %>% unlist()

  if (!is.null(external))  { # When external passages

    # get a, b, alpha, beta from MCEM with specific task.id
    a.par.external <- task.data %>% filter(task.id %in% external) %>% select(a) %>% c() %>% unlist()
    b.par.external <- task.data %>% filter(task.id %in% external) %>% select(b) %>% c() %>% unlist()
    alpha.par.external <- task.data %>% filter(task.id %in% external) %>% select(alpha) %>% c() %>% unlist()
    beta.par.external <- task.data %>% filter(task.id %in% external) %>% select(beta) %>% c() %>% unlist()
    max.counts.task.external <- task.data %>% filter(task.id %in% external) %>% select(max.counts) %>% c() %>% unlist()
  }

  # Using MCEM to calculate rho and vartau
  rho <- mean(MCEM$hyper.param$rho)
  vartau <- mean(MCEM$hyper.param$vartau)

  flog.info(paste(paste("Output", est),"Bootstrap"), name = "orfrlog")
  if (bootstrap != 100)
    flog.info(paste(paste("Set bootstrap K =", bootstrap),"bootstrap"), name = "orfrlog")

  if (est == "mle") {
    # Now, consider MLE as an example
    # Extract relevant latent param estimates
    Z.in <- c(WCPM$theta.mle,WCPM$tau.mle)

    I <- length(max.counts.task)

    K <- bootstrap
    Z.est <- matrix(rep(0,4*K),ncol = 4)

    new.data <- datasim.fixedZ(a.par,b.par,alpha.par,beta.par,vartau,rho,max.counts.task,I,Z.in,K)
    for (k in 1:K) {
      obs.counts <- as.array(new.data$Y[k,])
      # lgsec <- as.array(new.data$logT[k,])
      lgsec10 <- as.array(new.data$logT10[k,])

      tau.mle <- sum(alpha.par^2*(beta.par - lgsec10))/sum(alpha.par^2)
      mod.pd1 <- function(theta) {
        eta <- a.par*theta - b.par
        term1 <- sum(a.par*obs.counts*exp(dnorm(eta,log = TRUE)-pnorm(eta,log.p = TRUE)))
        term2 <- sum(a.par*(max.counts.task-obs.counts)*exp(dnorm(eta,log = TRUE)-pnorm(eta, lower.tail = FALSE, log.p = TRUE)))
        pd1 <- term1 - term2
        return(pd1)
      }

      theta.mle = Inf # for perfect case
      obs.counts.mle <- NA
      wcpm.mle <- NA
      k.theta <- NA
      se.wcpm.mle <- NA
      if (!is.null(Z.in)) {
        if (length(Z.in) > 0 & !is.infinite(Z.in[1])) { #for non-perfect case
          theta.mle <- uniroot(mod.pd1, c(-12, 12))$root
          # MLE WCPM score
          if (is.null(external)) { #internal
            obs.counts.mle <- sum(max.counts.task*pnorm(a.par*theta.mle - b.par))
            secs.mle <- sum(exp(beta.par - log(10) + log(max.counts.task) - tau.mle + ((1/alpha.par)^2)/2))
            wcpm.mle <- obs.counts.mle/secs.mle*60
            k.theta <- sum(a.par*max.counts.task*dnorm( a.par*theta.mle - b.par ))/sum(max.counts.task*pnorm( a.par*theta.mle - b.par ))
          } else { #external
            obs.counts.mle <- sum(max.counts.task.external*pnorm(a.par.external*theta.mle - b.par.external))
            secs.mle <- sum(exp(beta.par.external - log(10) + log(max.counts.task.external) - tau.mle + ((1/alpha.par.external)^2)/2))
            wcpm.mle <- obs.counts.mle/secs.mle*60
            k.theta <- sum(a.par.external*max.counts.task.external*dnorm( a.par.external*theta.mle - b.par.external ))/sum(max.counts.task.external*pnorm( a.par.external*theta.mle - b.par.external ))
          }
          Z.est[k,] <- c(theta.mle, tau.mle, wcpm.mle, k.theta)
        }
      }
    }
    se.mle <- apply(Z.est[,1:2],2,sd)
    wcpm.mle <- mean(Z.est[,3])
    k.theta <- mean(Z.est[,4])
    se.wcpm.mle <- wcpm.mle*(k.theta^2*se.mle[1]^2 + se.mle[2]^2)^0.5

    if (!is.null(WCPM)) {
      SE <- as.data.frame(cbind(do.call(cbind, WCPM),
                                bse.theta.mle=se.mle[1],
                                bse.tau.mle=se.mle[2],
                                bse.wcpm.mle=se.wcpm.mle))
      SE <- SE %>% select(person.id,
                          occasion,
                          group,
                          task.n,
                          max.counts.total,
                          obs.counts.obs,
                          secs.obs,
                          wcpm.obs,
                          tau.mle,
                          theta.mle,
                          se.tau.mle,
                          se.theta.mle,
                          obs.counts.mle,
                          secs.mle,
                          wcpm.mle,
                          se.wcpm.mle,
                          bse.theta.mle,
                          bse.tau.mle,
                          bse.wcpm.mle
      )
    } else {
      SE <- NULL
    }
  } else if (est == "eap") {
    # For QUAD
    # Extract relevant latent param estimates
    Z.in <- c(WCPM$theta.eap,WCPM$tau.eap)
    I <- length(max.counts.task)

    K <- bootstrap # for QUAD 500 should be default?
    Z.est <- matrix(rep(0,4*K),ncol = 4)

    new.data <- datasim.fixedZ(a.par,b.par,alpha.par,beta.par,vartau,rho,max.counts.task,I,Z.in,K)

    # Bivariate EAP for theta and tau
    cov <- rho*sqrt(vartau)
    prior <- list(mu = c(0,0), Sigma = matrix(c(1,cov,cov,vartau),2,2))
    #grid <- init.quad(Q = 2, prior, ip = 100, prune = TRUE)
    # ip should be 500, but set as 100 for test
    grid <- MultiGHQuad::init.quad(Q = 2, prior, ip = 500, prune = F)

    loglik <- function(z) {
      theta <- z[1]
      tau <- z[2]
      loglik.bi <- sum(dbinom(obs.counts, max.counts.task, pnorm((a.par*theta)-b.par), log = T)) +
        sum(dnorm(lgsec10, beta.par-tau, 1/alpha.par, log = T))
    }

    for (k in 1:K) {

      # get values
      obs.counts <- as.array(new.data$Y[k,])
      lgsec10 <- as.array(new.data$logT10[k,])

      ests.quad <- MultiGHQuad::eval.quad(loglik, grid)
      # QUAD WCPM score
      if (is.null(external)) { #internal
        obs.counts.quad <- sum(max.counts.task*pnorm(a.par*ests.quad[1] - b.par))
        secs.quad <- sum(exp(beta.par - log(10) + log(max.counts.task) - ests.quad[2] + ((1/alpha.par)^2)/2))
        wcpm.quad <- obs.counts.quad/secs.quad*60
        k.theta.quad <- sum(a.par*max.counts.task*dnorm( a.par*ests.quad[1] - b.par ))/sum(max.counts.task*pnorm( a.par*ests.quad[1] - b.par ))
      } else {
        obs.counts.quad <- sum(max.counts.task.external*pnorm(a.par.external*ests.quad[1] - b.par.external))
        secs.quad <- sum(exp(beta.par.external - log(10) + log(max.counts.task.external) - ests.quad[2] + ((1/alpha.par.external)^2)/2))
        wcpm.quad <- obs.counts.quad/secs.quad*60
        k.theta.quad <- sum(a.par.external*max.counts.task.external*dnorm( a.par.external*ests.quad[1] - b.par.external ))/sum(max.counts.task.external*pnorm( a.par.external*ests.quad[1] - b.par.external ))
      }
      # End of BiEAP
      Z.est[k,] <- c(ests.quad[1], ests.quad[2], wcpm.quad, k.theta.quad)
    }
    se.quad <- apply(Z.est[,1:2],2,sd)
    wcpm.quad <- mean(Z.est[,3])
    k.theta.quad <- mean(Z.est[,4])
    se.wcpm.quad <- wcpm.quad*(k.theta.quad^2*se.quad[1]^2 + se.quad[2]^2)^0.5

    SE <- as.data.frame(cbind(do.call(cbind, WCPM),
                              bse.theta.eap=se.quad[1],
                              bse.tau.eap=se.quad[2],
                              bse.wcpm.eap=se.wcpm.quad))
    SE <- SE %>% select(person.id,
                        occasion,
                        group,
                        task.n,
                        max.counts.total,
                        obs.counts.obs,
                        secs.obs,
                        wcpm.obs,
                        tau.eap,
                        theta.eap,
                        se.tau.eap,
                        se.theta.eap,
                        obs.counts.eap,
                        secs.eap,
                        wcpm.eap,
                        se.wcpm.eap,
                        bse.theta.eap,
                        bse.tau.eap,
                        bse.wcpm.eap
    )


  } else if (est == "map") {
    Z.in <- c(WCPM$theta.map,WCPM$tau.map)

    # get theta.mle
    WCPM_mle <- MCEM %>% run.scoring(person.data, task.data=task.data, cases=case, perfect.cases, est='mle', lo = -4, hi = 4, q = 100, kappa = 1, external=external, type="orf")
    # run.scoring(person.data, task.data=task.data, cases=case, perfect.cases, est=Estimator, lo = -4, hi = 4, q = 100, kappa = 1, external=external, type="orf")

    I <- length(max.counts.task)
    K <- bootstrap
    Z.est <- matrix(rep(0,4*K),ncol = 4)

    new.data <- datasim.fixedZ(a.par,b.par,alpha.par,beta.par,vartau,rho,max.counts.task,I,Z.in,K)
    #print(paste("kappa=", kappa))

    est.eqs <- function(latent.parms) {
      theta <- latent.parms[1]
      tau <- latent.parms[2]
      eta <- a.par*theta - b.par

      ee1 <- -1/(kappa^2*(1-rho^2))*(theta-rho/sqrt(vartau)*tau) +
        sum(a.par*obs.counts*exp(dnorm(eta,log = TRUE)-pnorm(eta,log.p = TRUE))) -
        sum(a.par*(max.counts.task-obs.counts)*exp(dnorm(eta,log = TRUE)-pnorm(eta, lower.tail = FALSE, log.p = TRUE)))
      ee2 <- -1/(kappa^2*(1-rho^2))*(tau/vartau-rho/sqrt(vartau)*theta) -
        sum(alpha.par^2*(lgsec10 - beta.par + tau))
      ee <- c(ee1,ee2)
      return(ee)
    }

    for (k in 1:K) {
      obs.counts <- as.array(new.data$Y[k,])
      # lgsec <- as.array(new.data$logT[k,])
      lgsec10 <- as.array(new.data$logT10[k,])

      tau.mle <- sum(alpha.par^2*(beta.par - lgsec10))/sum(alpha.par^2)
      mod.pd1 <- function(theta) {
        eta <- a.par*theta - b.par
        term1 <- sum(a.par*obs.counts*exp(dnorm(eta,log = TRUE)-pnorm(eta,log.p = TRUE)))
        term2 <- sum(a.par*(max.counts.task-obs.counts)*exp(dnorm(eta,log = TRUE)-pnorm(eta, lower.tail = FALSE, log.p = TRUE)))
        pd1 <- term1 - term2
        return(pd1)
      }

      if (!is.infinite(WCPM_mle$theta.mle)) #for non-perfect case
        theta.mle <- uniroot(mod.pd1, c(-12, 12))$root
      else #for perfect case
        theta.mle = Inf

      in.vals <- c(max(-5,min(5,theta.mle)),max(-5*sqrt(vartau),min(5*sqrt(vartau),tau.mle)))
      ests.map <- rootSolve::multiroot(est.eqs, in.vals)$root
      # MAP WCPM score
      if (is.null(external)) { #internal
        obs.counts.map <- sum(max.counts.task*pnorm(a.par*ests.map[1] - b.par))
        secs.map <- sum(exp(beta.par - log(10) + log(max.counts.task) - ests.map[2] + ((1/alpha.par)^2)/2))
        wcpm.map <- obs.counts.map/secs.map*60
        k.theta.map <- sum(a.par*max.counts.task*dnorm( a.par*ests.map[1] - b.par ))/sum(max.counts.task*pnorm( a.par*ests.map[1] - b.par ))
      } else {
        obs.counts.map <- sum(max.counts.task.external*pnorm(a.par.external*ests.map[1] - b.par.external))
        secs.map <- sum(exp(beta.par.external - log(10) + log(max.counts.task.external) - ests.map[2] + ((1/alpha.par.external)^2)/2))
        wcpm.map <- obs.counts.map/secs.map*60
        k.theta.map <- sum(a.par.external*max.counts.task.external*dnorm( a.par.external*ests.map[1] - b.par.external ))/sum(max.counts.task.external*pnorm( a.par.external*ests.map[1] - b.par.external ))
      }
      Z.est[k,] <- c(ests.map[1], ests.map[2], wcpm.map, k.theta.map)
    }
    se.map <- apply(Z.est,2,sd)
    wcpm.map <- mean(Z.est[,3])
    k.theta.map <- mean(Z.est[,4])
    se.wcpm.map <- wcpm.map*(k.theta.map^2*se.map[1]^2 + se.map[2]^2)^0.5

    SE <- as.data.frame(cbind(do.call(cbind, WCPM),
                              bse.theta.map=se.map[1],
                              bse.tau.map=se.map[2],
                              bse.wcpm.map=se.wcpm.map))
    SE <- SE %>% select(person.id,
                        occasion,
                        group,
                        task.n,
                        max.counts.total,
                        obs.counts.obs,
                        secs.obs,
                        wcpm.obs,
                        tau.map,
                        theta.map,
                        se.tau.map,
                        se.theta.map,
                        obs.counts.map,
                        secs.map,
                        wcpm.map,
                        se.wcpm.map,
                        bse.theta.map,
                        bse.tau.map,
                        bse.wcpm.map
    )

  }

  flog.info("End getBootstrapSE process", name = "orfrlog")
  return(SE)
}
