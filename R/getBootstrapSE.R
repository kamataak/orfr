#' Get bootstrap SE
#' case is a single stu_season_id
#'
#' Added MAP function 07/14/2021
#' Modified a bug of MLE 07/23/2021
#' Modified EAP 10/28/2021
#'
#' @param object - mcem class object
#' @param stu.data - student response data
#' @param case - case number
#' @param est - SE type.(MLE, EAP, and MAP.) default MAP
#' @param perfect.cases - perfect accurate case
#' @param kappa - Default kappa = 1, better be 5
#' @param bootstrap - K number of bootstrap, default is 100
#'
#' @return SE dataset
getBootstrapSE <- function (object, stu.data, case=NA, perfect.cases, est="map", kappa=1, bootstrap=100) {
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
  pass.data <- MCEM$pass.param
  WCPM <- MCEM %>% run.wcpm(stu.data, pass.data=pass.data, cases=case, perfect.cases, est=Estimator, lo = -4, hi = 4, q = 100, kappa = 1)

  # Extract relevant parameters for given case
  # stu.dat01 <- stu.data %>% filter(stu_season_id2==case)
  # pass.read <- stu.dat01 %>% select(passage_id)
  # pass.dat01 <- pass.data %>% semi_join(pass.read, by = "passage_id")
  # n.pass <- nrow(pass.dat01)
  # numwords.total <- stu.dat01 %>% select(nwords.p) %>% c() %>% unlist() %>% sum()
  # grade <- stu.dat01 %>% select(grade) %>% c() %>% unlist %>% unique()

  case_split <- unlist(str_split(case, "_"))
  stu.dat01 <- stu.data %>% filter(stu.data$student.id==case_split[1], stu.data$occasion==case_split[2])
  pass.read <- stu.dat01 %>% select(passage.id)
  pass.dat01 <- pass.data %>% semi_join(pass.read, by = "passage.id")
  n.pass <- nrow(pass.dat01)
  numwords.total <- stu.dat01 %>% select(nwords.p) %>% c() %>% unlist() %>% sum()
  grade <- stu.dat01 %>% select(grade) %>% c() %>% unlist %>% unique()


  numwords.pass <- stu.dat01 %>% select(nwords.p) %>% c() %>% unlist()

  a.par <- pass.dat01 %>% select(a) %>% c() %>% unlist()
  b.par <- pass.dat01 %>% select(b) %>% c() %>% unlist()
  alpha.par <- pass.dat01 %>% select(alpha) %>% c() %>% unlist()
  beta.par <- pass.dat01 %>% select(beta) %>% c() %>% unlist()

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

    I <- length(numwords.pass)

    K <- bootstrap
    Z.est <- matrix(rep(0,4*K),ncol = 4)

    new.data <- datasim.fixedZ(a.par,b.par,alpha.par,beta.par,vartau,rho,numwords.pass,I,Z.in,K)
    for (k in 1:K) {
      wrc <- as.array(new.data$Y[k,])
      # lgsec <- as.array(new.data$logT[k,])
      lgsec10 <- as.array(new.data$logT10[k,])

      tau.mle <- sum(alpha.par^2*(beta.par - lgsec10))/sum(alpha.par^2)
      mod.pd1 <- function(theta) {
        eta <- a.par*theta - b.par
        term1 <- sum(a.par*wrc*exp(dnorm(eta,log = TRUE)-pnorm(eta,log.p = TRUE)))
        term2 <- sum(a.par*(numwords.pass-wrc)*exp(dnorm(eta,log = TRUE)-pnorm(eta, lower.tail = FALSE, log.p = TRUE)))
        pd1 <- term1 - term2
        return(pd1)
      }

      theta.mle = Inf # for perfect case
      wrc.mle <- NA
      wcpm.mle <- NA
      k.theta <- NA
      se.wcpm.mle <- NA
      if (!is.infinite(Z.in[1])) { #for non-perfect case
        theta.mle <- uniroot(mod.pd1, c(-12, 12))$root
        # MLE WCPM score
        wrc.mle <- sum(numwords.pass*pnorm(a.par*theta.mle - b.par))
        secs.mle <- sum(exp(beta.par - log(10) + log(numwords.pass) - tau.mle + ((1/alpha.par)^2)/2))
        wcpm.mle <- wrc.mle/secs.mle*60
        k.theta <- sum(a.par*numwords.pass*dnorm( a.par*theta.mle - b.par ))/sum(numwords.pass*pnorm( a.par*theta.mle - b.par ))
      }
      Z.est[k,] <- c(theta.mle, tau.mle, wcpm.mle, k.theta)
    }
    se.mle <- apply(Z.est[,1:2],2,sd)
    wcpm.mle <- mean(Z.est[,3])
    k.theta <- mean(Z.est[,4])
    se.wcpm.mle <- wcpm.mle*(k.theta^2*se.mle[1]^2 + se.mle[2]^2)^0.5

    # Original SE ests
    #    orgSE <- c(WCPM$se.theta.mle,WCPM$se.tau)
    # SE <- tibble(stu_season_id=case, se.theta.mle=se.mle[1], se.tau.mle=se.mle[2], se.wcpm.mle=se.wcpm.mle,
    #              org.se.theta.mle=WCPM$se.theta.mle, org.se.tau.mle=WCPM$se.tau.mle, org.se.wcpm.mle=WCPM$se.wcpm.mle)
    SE <- as.data.frame(cbind(do.call(cbind, WCPM),
                              bse.theta.mle=se.mle[1],
                              bse.tau.mle=se.mle[2],
                              bse.wcpm.mle=se.wcpm.mle))
    SE <- SE %>% select(student.id,
                        occasion,
                        grade,
                        n.pass,
                        numwords.total,
                        wrc.obs,
                        secs.obs,
                        wcpm.obs,
                        tau.mle,
                        theta.mle,
                        se.tau.mle,
                        se.theta.mle,
                        wrc.mle,
                        secs.mle,
                        wcpm.mle,
                        se.wcpm.mle,
                        bse.theta.mle,
                        bse.tau.mle,
                        bse.wcpm.mle
    )

  } else if (est == "eap") {
    # For QUAD
    # Extract relevant latent param estimates
    Z.in <- c(WCPM$theta.eap,WCPM$tau.eap)
    I <- length(numwords.pass)

    K <- bootstrap # for QUAD 500 should be default?
    Z.est <- matrix(rep(0,4*K),ncol = 4)

    new.data <- datasim.fixedZ(a.par,b.par,alpha.par,beta.par,vartau,rho,numwords.pass,I,Z.in,K)

    # Bivariate EAP for theta and tau
    cov <- rho*sqrt(vartau)
    prior <- list(mu = c(0,0), Sigma = matrix(c(1,cov,cov,vartau),2,2))
    #grid <- init.quad(Q = 2, prior, ip = 100, prune = TRUE)
    # ip should be 500, but set as 100 for test
    grid <- MultiGHQuad::init.quad(Q = 2, prior, ip = 500, prune = F)

    loglik <- function(z) {
      theta <- z[1]
      tau <- z[2]
      loglik.bi <- sum(dbinom(wrc, numwords.pass, pnorm((a.par*theta)-b.par), log = T)) +
        sum(dnorm(lgsec10, beta.par-tau, 1/alpha.par, log = T))
    }

    for (k in 1:K) {

      # get values
      wrc <- as.array(new.data$Y[k,])
      lgsec10 <- as.array(new.data$logT10[k,])

      ests.quad <- MultiGHQuad::eval.quad(loglik, grid)
      # QUAD WCPM score
      wrc.quad <- sum(numwords.pass*pnorm(a.par*ests.quad[1] - b.par))
      secs.quad <- sum(exp(beta.par - log(10) + log(numwords.pass) - ests.quad[2] + ((1/alpha.par)^2)/2))
      wcpm.quad <- wrc.quad/secs.quad*60
      k.theta.quad <- sum(a.par*numwords.pass*dnorm( a.par*ests.quad[1] - b.par ))/sum(numwords.pass*pnorm( a.par*ests.quad[1] - b.par ))
      Z.est[k,] <- c(ests.quad[1], ests.quad[2], wcpm.quad, k.theta.quad)
      # End of BiEAP
    }
    se.quad <- apply(Z.est[,1:2],2,sd)
    wcpm.quad <- mean(Z.est[,3])
    k.theta.quad <- mean(Z.est[,4])
    se.wcpm.quad <- wcpm.quad*(k.theta.quad^2*se.quad[1]^2 + se.quad[2]^2)^0.5

    # SE <- tibble(stu_season_id=case, se.theta.eap=se.quad[1],se.tau.eap=se.quad[2], se.wcpm.eap=se.wcpm.quad,
    #              org.se.theta.eap=WCPM$se.theta.eap, org.se.tau.eap=WCPM$se.tau.eap, org.se.wcpm.eap=WCPM$se.wcpm.eap)
    SE <- as.data.frame(cbind(do.call(cbind, WCPM),
                              bse.theta.eap=se.quad[1],
                              bse.tau.eap=se.quad[2],
                              bse.wcpm.eap=se.wcpm.quad))
    SE <- SE %>% select(student.id,
                        occasion,
                        grade,
                        n.pass,
                        numwords.total,
                        wrc.obs,
                        secs.obs,
                        wcpm.obs,
                        tau.eap,
                        theta.eap,
                        se.tau.eap,
                        se.theta.eap,
                        wrc.eap,
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
    WCPM_mle <- MCEM %>% run.wcpm(stu.data, pass.data=pass.data, cases=case, perfect.cases, est='mle', lo = -4, hi = 4, q = 100, kappa = 1)

    I <- length(numwords.pass)
    K <- bootstrap
    Z.est <- matrix(rep(0,4*K),ncol = 4)

    new.data <- datasim.fixedZ(a.par,b.par,alpha.par,beta.par,vartau,rho,numwords.pass,I,Z.in,K)
    #print(paste("kappa=", kappa))

    est.eqs <- function(latent.parms) {
      theta <- latent.parms[1]
      tau <- latent.parms[2]
      eta <- a.par*theta - b.par

      ee1 <- -1/(kappa^2*(1-rho^2))*(theta-rho/sqrt(vartau)*tau) +
        sum(a.par*wrc*exp(dnorm(eta,log = TRUE)-pnorm(eta,log.p = TRUE))) -
        sum(a.par*(numwords.pass-wrc)*exp(dnorm(eta,log = TRUE)-pnorm(eta, lower.tail = FALSE, log.p = TRUE)))
      ee2 <- -1/(kappa^2*(1-rho^2))*(tau/vartau-rho/sqrt(vartau)*theta) -
        sum(alpha.par^2*(lgsec10 - beta.par + tau))
      ee <- c(ee1,ee2)
      return(ee)
    }

    for (k in 1:K) {
      wrc <- as.array(new.data$Y[k,])
      # lgsec <- as.array(new.data$logT[k,])
      lgsec10 <- as.array(new.data$logT10[k,])

      tau.mle <- sum(alpha.par^2*(beta.par - lgsec10))/sum(alpha.par^2)
      mod.pd1 <- function(theta) {
        eta <- a.par*theta - b.par
        term1 <- sum(a.par*wrc*exp(dnorm(eta,log = TRUE)-pnorm(eta,log.p = TRUE)))
        term2 <- sum(a.par*(numwords.pass-wrc)*exp(dnorm(eta,log = TRUE)-pnorm(eta, lower.tail = FALSE, log.p = TRUE)))
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
      wrc.map <- sum(numwords.pass*pnorm(a.par*ests.map[1] - b.par))
      secs.map <- sum(exp(beta.par - log(10) + log(numwords.pass) - ests.map[2] + ((1/alpha.par)^2)/2))
      wcpm.map <- wrc.map/secs.map*60
      k.theta.map <- sum(a.par*numwords.pass*dnorm( a.par*ests.map[1] - b.par ))/sum(numwords.pass*pnorm( a.par*ests.map[1] - b.par ))
      Z.est[k,] <- c(ests.map[1], ests.map[2], wcpm.map, k.theta.map)
    }
    se.map <- apply(Z.est,2,sd)
    wcpm.map <- mean(Z.est[,3])
    k.theta.map <- mean(Z.est[,4])
    se.wcpm.map <- wcpm.map*(k.theta.map^2*se.map[1]^2 + se.map[2]^2)^0.5

    # SE <- tibble(stu_season_id=case, se.theta.map=se.map[1],se.tau.map=se.map[2], se.wcpm.map,
    #              org.se.theta.mle=WCPM$se.theta.mle, org.se.tau.mle=WCPM$se.tau.mle, org.se.wcpm.map=WCPM$se.wcpm.map)
    SE <- as.data.frame(cbind(do.call(cbind, WCPM),
                              bse.theta.map=se.map[1],
                              bse.tau.map=se.map[2],
                              bse.wcpm.map=se.wcpm.map))
    SE <- SE %>% select(student.id,
                        occasion,
                        grade,
                        n.pass,
                        numwords.total,
                        wrc.obs,
                        secs.obs,
                        wcpm.obs,
                        tau.map,
                        theta.map,
                        se.tau.map,
                        se.theta.map,
                        wrc.map,
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
