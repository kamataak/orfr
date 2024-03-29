#' This file includes the Base functions of orfr package.
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
#' This is run.mcem function.
#' Update Memo:
#' 04/29/2021 Modified the mcem function
#'            based on Nelis's updated.
#' 10/28/2021 Modified the mcem output
#'
#' @param Y       = n x I matrix of reading scores -- missingness allowed
#' @param logT10  = n x I matrix of log10(reading times) -- missingness allowed
#' @param N       = vector of passage lengths
#' @param I       = number of passages
#' @param k.in    = number of imputations, default is 5
#' @param reps.in = number of Monte-Carlo iterations, default is 2
#' @param ests.in = if not give, mom function will be called and get est.in output
#' @param verbose - boolean, if shows the summary, default is FALSE
#'
#' @import mvtnorm
#' @import tidyverse
#' @import nleqslv
#' @import greekLetters
#'
#'
#' @return mcem list
#' a,b = parameters controlling binomial success probabilities, each length I
#' alpha,beta = parameters controlling reading times, each length I
#' var_tau = variance of latent reading ability tau
#' rho = correlation between two latent variables
run.mcem <- function(Y,logT10,N,I,k.in=5,reps.in=2,ests.in=NA,verbose=FALSE) {

  # loading logger
  log.initiating()

  neg_logratio <- function(data,par) {
    Y <- data[1,]
    N <- data[2,]
    a <- data[3,]
    b <- data[4,]
    alpha <- data[5,]
    Ik <- dim(data)[2]
    z <- par
    LP <- 0
    for (k in 1:Ik) {
      # LP <- LP + log(alpha[k]) + lchoose(N[k],Y[k]) + Y[k]*log(pnorm(a[k]*z-b[k])) + (N[k]-Y[k])*log(1-pnorm(a[k]*z-b[k]))
      LP <- LP + log(alpha[k]) + dbinom(Y[k],N[k],pnorm(a[k]*z-b[k]),log = T)
    }
    LP <- -LP
    return(LP)
  }
  gamma_multiplier <- function(Y,N,a,b,alpha,z.in) {
    if (sum(N-Y)<0.01) { #Condition for someone with a perfect score on all items
      gammax <- exp(sum(log(alpha)))
      zmax <- Inf
    } else {
      if (missing(z.in)) {z.in = 0}
      #      print("here: ")
      #      print(Y)
      #      print(z.in)
      zgam <- optim(z.in, fn = neg_logratio, method="Nelder-Mead", control = list(warn.1d.NelderMead = FALSE), data = rbind(Y,N,a,b,alpha))
      #zgam <- optim(z.in, fn = neg_logratio, method="BFGS", data = rbind(Y,N,a,b,alpha))
      # Note on above line: BFGS tends to be faster, but about
      # 1 out of every 100 cases doesn't work so for now the
      # code uses the slower Nelder-Mead algorithm
      zmax <- zgam$par
      gammax <- exp(-zgam$value)
    }
    gamval <- list(zmax = zmax, gammax = gammax)
    return(gamval)
  }
  imputation_code <- function(Y,logT10,N,a,b,alpha,beta,sigma_tau2,sigma_theta_tau,K,z.in) {
    n <- dim(Y)[1]
    if (missing(z.in)) {z.in = matrix(rep(0,n),nrow = n)}
    zmax <- matrix(rep(0,n),nrow = n)
    gammax <- matrix(rep(0,n),nrow = n)
    mu1 <- matrix(rep(0,n), nrow = n)
    A <- matrix(rep(0,n), nrow = n)
    sigma1_2 <- matrix(rep(0,n), nrow = n)
    sigma2_2 <- matrix(rep(0,n), nrow = n)
    # Parameter setup for imputation step
    for (i in 1:n) {
      index <- which(is.na(Y[i,])==FALSE)
      gamval <- gamma_multiplier(Y[i,index],N[index],a[index],b[index],alpha[index],z.in[i])
      gammax[i] <- gamval$gammax
      zmax[i] <- gamval$zmax
      A[i] <- sum(alpha[index]^2)
      mu1[i] <- -sigma_theta_tau/(1+A[i]*sigma_tau2)*sum(alpha[index]^2.*(logT10[i,index]-beta[index]))
      sigma1_2[i] <- 1/(1+(A[i]*sigma_theta_tau^2)/(1+A[i]*(sigma_tau2-sigma_theta_tau^2)))
      sigma2_2[i] = (sigma_tau2-sigma_theta_tau^2)/((sigma_tau2-sigma_theta_tau^2)*A[i]+1)
    }
    # Imputation step
    Z1 <- matrix(rep(0,n*K),nrow = n, ncol = K)
    mu2 <- matrix(rep(0,n*K), nrow = n, ncol = K)
    tau <- matrix(rep(0,n*K), nrow = n, ncol = K)
    for (i in 1:n) {
      index <- which(is.na(Y[i,])==FALSE)
      for (k in 1:K) {
        check <- 0
        while (check<1) {
          Z1[i,k] <- mu1[i] + sqrt(sigma1_2[i])*rnorm(1)
          U <- runif(1)
          LP <- exp(-neg_logratio(rbind(Y[i,index],N[index],a[index],b[index],alpha[index]),Z1[i,k]))/gammax[i]
          if (U<LP) {
            check <- 2
          }
        }
        mu2[i,k] = -((sigma_tau2-sigma_theta_tau^2)*(sum(alpha[index]^2*(logT10[i,index]-beta[index])))-sigma_theta_tau*Z1[i,k])/((sigma_tau2-sigma_theta_tau^2)*A[i]+1)
        tau[i,k] = mu2[i,k] + sqrt(sigma2_2[i])*rnorm(1)
      }
    }
    theta <- Z1
    imputes <- list(theta = theta, tau = tau, z.opt = zmax)
  }
  est_eq_ab_missing <- function(x,Y,Ind,N,theta){
    a <- x[1]
    b <- x[2]
    Ratio1 <- dnorm(a*theta-b)/pnorm(a*theta-b)
    Ratio2 <- dnorm(a*theta-b)/(1-pnorm(a*theta-b))
    # Conditions where function may become ill-defined numerically
    Ratio1[Ratio1==Inf] <- 0
    Ratio1[Ratio1==0] <- max(Ratio1)
    Ratio2[Ratio2==Inf] <- 0
    Ratio2[Ratio2==0] <- max(Ratio2)

    E1 <- sum(Ind*Y*rowMeans(theta*Ratio1) - Ind*(N-Y)*rowMeans(theta*Ratio2))
    E2 <- sum(Ind*Y*rowMeans(Ratio1) - Ind*(N-Y)*rowMeans(Ratio2))
    E <-c(E1,E2)
    return(E)
  }
  function_s12s22_to_min <- function(data,par) {
    theta <- data$theta
    tau <- data$tau
    s12 <- par
    s22 <- s12^2*(1+mean(theta^2)) - 2*s12*mean(theta*tau) + mean(tau^2)
    F1 <- 0.5*log(s22-s12^2)
    F2 <- 0.5/sqrt(s22-s12^2)*(s22*mean(theta^2)-2*s12*mean(theta*tau)+mean(tau^2))
    Fval <- F1+F2
    return(Fval)
  }
  MCEM_algorithm_one_iteration <- function(Y,logT10,N,I,ests,K,z.in) {
    if (missing(z.in)) {
      n <- dim(Y)[1]
      z.in <- matrix(rep(0,n),nrow = n)
    }
    EMimps <- imputation_code(Y,logT10,N,ests$a,ests$b,ests$alpha,ests$beta,ests$vartau,(ests$rho)*sqrt(ests$vartau),K)
    Ind <- 1 - is.na(Y)
    Y1 <- Y
    logT1 <- logT10
    Y1[Ind==0] <- 0
    logT1[Ind==0] <- 0
    m1 <- colSums(Ind)
    EM.a <- rep(0,length.out = I)
    EM.b <- rep(0,length.out = I)
    EM.beta <- rep(0,length.out = I)
    EM.alpha <- rep(0,length.out = I)
    for (k0 in 1:I) {
      Y0 <- Y1[,k0]
      logT0 <- logT1[,k0]
      Ind0 <- Ind[,k0]
      m0 <- m1[k0]
      f.out <- nleqslv(c(ests$a[k0],ests$b[k0]), fn = est_eq_ab_missing, Y = Y0, Ind = Ind0, N = N[k0], theta = EMimps$theta)
      EM.a[k0] <- f.out$x[1]
      EM.b[k0] <- f.out$x[2]
      EM.beta[k0] <- (sum(Ind0*logT0)+sum(Ind0*rowMeans(EMimps$tau)))/m0
      alp2_inv <- sum(Ind0*rowMeans((logT0-EM.beta[k0]+EMimps$tau)^2))/m0
      EM.alpha[k0] <- 1/sqrt(alp2_inv)
    }
    s12.find <- optim(ests$rho*ests$vartau, fn = function_s12s22_to_min, method = "BFGS", data = EMimps)
    s12.min <- s12.find$par
    EM.vartau <- s12.min^2*(1+mean(EMimps$theta^2)) - 2*s12.min*mean(EMimps$theta*EMimps$tau) + mean(EMimps$tau^2)
    EM.rho <- s12.min/sqrt(EM.vartau)
    EM.ests <- list(a = EM.a, b = EM.b, alpha = EM.alpha, beta = EM.beta, vartau = EM.vartau, rho = EM.rho, z.opt = EMimps$z.opt)
    return(EM.ests)
  }

  nK <- length(k.in)
  #  if (!is.null(Z.in))
  if (length(ests.in) == 1) {
    if (is.na(ests.in)) {
      ests.in <- mom(Y,logT10,N,I)
    }
  }

  # MCEM algorithm can't initiate if any alpha values = inf
  # Somewhat artificial solution to give starting values to EM:
  alpha.check <- ests.in$alpha
  infIndex0 <- which(is.infinite(alpha.check)==0)
  infIndex1 <- which(is.infinite(alpha.check)==1)
  alpha.check[infIndex1] <- max(alpha.check[infIndex0])*10
  ests.in$alpha <- alpha.check

  n <- dim(Y)[1]
  z.in <- matrix(rep(0,n), nrow = n)

  total.K <- rep(k.in[1],reps.in[1])

  if (nK > 1) {
    for (jj in 2:nK) {
      total.K <- c(total.K,rep(k.in[jj],reps.in[jj]))
    }
  }
  JJ <- length(total.K)

  a.store <- matrix(rep(0,JJ*I),nrow = JJ)
  alpha.store <- matrix(rep(0,JJ*I),nrow = JJ)
  b.store <- matrix(rep(0,JJ*I),nrow = JJ)
  beta.store <- matrix(rep(0,JJ*I),nrow = JJ)
  rho.store <- matrix(rep(0,JJ),nrow = JJ)
  vartau.store <- matrix(rep(0,JJ),nrow = JJ)
  se_a.store <- matrix(rep(0,JJ),nrow = JJ)
  se_b.store <- matrix(rep(0,JJ),nrow = JJ)
  se_alpha.store <- matrix(rep(0,JJ),nrow = JJ)
  se_beta.store <- matrix(rep(0,JJ),nrow = JJ)

  for (jj in 1:JJ) {

    EM.iter <- MCEM_algorithm_one_iteration(Y,logT10,N,I,ests.in,total.K[jj],z.in)
    ests.in <- EM.iter
    a.store[jj,] <- ests.in$a
    b.store[jj,] <- ests.in$b
    alpha.store[jj,] <- ests.in$alpha
    beta.store[jj,] <- ests.in$beta
    vartau.store[jj] <- ests.in$vartau
    rho.store[jj] <- ests.in$rho
    z.in <- ests.in$z.opt
  }

  mean_a = a.store[JJ,]
  mean_b = b.store[JJ,]
  mean_alpha = alpha.store[JJ,]
  mean_beta = beta.store[JJ,]

  MCEM.ests <- list(
    a = mean_a,
    b = mean_b,
    alpha = mean_alpha,
    beta = mean_beta,
    numwords.p = N,
    vartau = vartau.store[JJ,],
    rho = rho.store[JJ,]
  )

  return(MCEM.ests)
}
mom <- function(Y,logT10,N,I) {

  evaluate_rho <- function(data,par) {
    R <- data[1]
    Q <- data[2]
    r <- par
    sigma <- matrix(c(1,r,r,1), ncol=2)
    MD <- (R - pmvnorm(upper = c(Q,Q), mean = rep(0,2), sigma = sigma))^2
    return(MD)
  }

  reading_data_parms_MOM <- function(Y,N) {
    Y.bar <- mean(Y)
    S2 <- mean((Y-Y.bar)^2)
    if (Y.bar<N) {
      Q <- qnorm(Y.bar/N)
    } else {
      Q <- qnorm((Y.bar+0.05)/(N+0.1))
    }
    R <- (S2 + Y.bar^2 - Y.bar)/(N*(N-1))
    parms.in <- 0.5
    rho.min <- optim(parms.in, fn = evaluate_rho, method = "Brent", data = c(R,Q), lower = 0, upper = 1)
    a.out <- sqrt(rho.min$par/(1-rho.min$par))
    b.out <- -Q*sqrt(1+a.out^2)
    parms <- list(a = a.out, b = b.out)
  }

  n <- dim(Y)[1]
  nancov_T <- cov(logT10,use = "pairwise.complete.obs")
  a.out <- matrix(nrow = 1, ncol = I)
  b.out <- matrix(nrow = 1, ncol = I)
  for (i in 1:I) {
    ab.out <- reading_data_parms_MOM(na.omit(Y[,i]),N[i])
    a.out[i] <- ab.out$a
    b.out[i] <- ab.out$b
  }
  beta.out <- apply(logT10,2,mean,na.rm = TRUE)
  count <- 1
  alpha.inv <- numeric(0)
  s22.comp <- numeric(0)
  for (i in 1:I) {
    if (i == I) {
      break
    }
    for (j in (i+1):I) {
      alpha.inv <- rbind(alpha.inv,matrix(rep(0,I+1),nrow = 1, ncol = (I+1)))
      alpha.inv[count,i] <- max(0,nancov_T[i,i]-nancov_T[i,j])
      alpha.inv[count,j] <- max(0,nancov_T[j,j]-nancov_T[i,j])
      s22.comp <- rbind(s22.comp,nancov_T[i,j])
      count <- count + 1
    }
  }
  alpha.inv_out <- apply(alpha.inv[,1:I],2,sum,na.rm=TRUE)/(I-apply(is.na(alpha.inv[,1:I]),2,sum)-1)
  alpha.out <- 1/sqrt(alpha.inv_out)
  vartau.out <- max(0,mean(s22.comp,na.rm = TRUE))
  s12.calc <- matrix(rep(0,I),nrow = I)
  for (j in 1:I) {
    CV <- cov(Y[,j],logT10[,j],use = "pairwise.complete.obs")
    s12.calc[j] <- -CV/N[j]*sqrt(2*pi)*sqrt(a.out[j]^2+1)/a.out[j]*exp(0.5*b.out[j]^2/(1+a.out[j]^2))
    s12.calc[j] <- max(-sqrt(vartau.out),min(sqrt(vartau.out),s12.calc[j]))
  }
  s12.out <- mean(s12.calc)
  if (vartau.out>0) {s12.cor <- s12.out/sqrt(vartau.out)} else {s12.cor <- 0}
  MOMests <- list(a = a.out, b = b.out, alpha = alpha.out, beta = beta.out, vartau = vartau.out, rho = s12.cor)
}
sim_reading_data <- function(a,b,alpha,beta,var_tau,rho,N,I,n,missing_prop,max_missing){
  # beta is already on the logT10 scale
  s12 <- rho*sqrt(var_tau)
  sigma <- matrix(c(1,s12,s12,var_tau), ncol=2)
  Z.latent <- rmvnorm(n=n, mean=c(0,0), sigma=sigma)
  Y <- matrix(nrow = n, ncol = I)
  logT10 <- matrix(nrow = n, ncol = I)
  S <- matrix(nrow = n, ncol = I)
  p.store <- matrix(nrow = n, ncol = I)
  for (i in 1:n) {
    is.missing <- as.double(runif(1)<missing_prop)
    num.missing <- ceiling(max_missing*runif(1))
    Missing <- is.missing*num.missing
    Index <- sample(x=1:I, size = Missing, replace = FALSE)
    p <- pnorm(a*Z.latent[i,1]-b)
    Y.in <- mapply(function(n,p) rbinom(1,n,p), N, p)
    Y.in[Index] <- NA
    logT10.in <- 1/alpha*rnorm(I)+beta-Z.latent[i,2]
    logT10.in[Index] <- NA
    S.in <- rep(1, times = I)
    S.in[Index] <- 0
    Y[i,] <- Y.in
    logT10[i,] <- logT10.in
    S[i,] = S.in
    p.store[i,] <- p
  }
  data <- list(Y = Y, logT10 = logT10, S = S, Z = Z.latent, p = p.store)
  return(data)
}
numerical.cov <- function(Y,logT10,N,I,parms,h.val,M) {

  if (missing(M)) {M <- 50}
  if (missing(h.val)) {h.val <- 1e-8}

  logT10 <- as.matrix(logT10)

  cond.distribution <- function(Y,logT10,N,I,parms,Z) {

    a <- parms[1:I]
    b <- parms[(I+1):(2*I)]
    alpha <- parms[(2*I+1):(3*I)]
    beta <- parms[(3*I+1):(4*I)]
    vartau <- parms[4*I+1]
    rho <- parms[4*I+2]

    n <- dim(Y)[1]

    Ind <- 1 - is.na(Y)
    m1 <- colSums(Ind)
    Y[Ind==0] <- 0
    logT10[Ind==0] <- 0

    M <- dim(Z)[1]
    Z1 <- Z[,1]
    Z2 <- Z[,2]

    a.mat <- matrix(a,nrow=n,ncol=I,byrow=T)
    b.mat <- matrix(b,nrow=n,ncol=I,byrow=T)
    alpha.mat <- matrix(alpha,nrow=n,ncol=I,byrow=T)
    beta.mat <- matrix(beta,nrow=n,ncol=I,byrow=T)
    N.mat <- matrix(N,nrow=n,ncol=I,byrow=T)

    log.f <- matrix(0,n,M)

    here.is <- Ind > 0

    for (m in 1:M) {

      logL1 <- rowSums(dbinom(Y,N.mat,pnorm(a.mat*Z1[m]-b.mat),log=T)*(Ind>0))
      #logL2B <- rowSums((log(alpha.mat)+dnorm(alpha.mat*(logT10-beta.mat+rho*sqrt(vartau)*Z1[m]+sqrt(vartau*(1-rho^2))*Z2[m]),log=T))*(Ind>0))
      result <- logT10-beta.mat+rho*sqrt(vartau)*Z1[m]+sqrt(vartau*(1-rho^2))*Z2[m]

      sd <-  1/alpha.mat
      #      print(sapply(result, class))
      #      print(sapply(sd, class))

      dd <- dnorm(result,sd=sd,log=T)

      logL2 <- rowSums(dd*(Ind>0))
      log.f[,m] <- logL1 + logL2

    }

    f <- rowMeans(exp(log.f))

    logL <- log(f)
    logL[is.infinite(logL)] <- -1e-10

    return(logL)

  }

  z <- qnorm((1:M)/(M+1))
  Z <- expand.grid(z,z)

  n <- dim(Y)[1]

  score.mat <- matrix(0,nrow=n,ncol=(4*I+2))

  for (i in 1:(4*I+2)) {
    h <- rep(0,4*I+2)
    h[i] <- h.val
    score.mat[,i] <- (cond.distribution(Y,logT10,N,I,parms+h,Z)-cond.distribution(Y,logT10,N,I,parms-h,Z))/(2*h.val)
  }

  I.incomp <- (n-1)*cov(score.mat) + colSums(score.mat)%*%t(colSums(score.mat))

  I.eigen <- eigen(I.incomp)$values
  if (sum(abs(Im(I.eigen)))>0) {pd.check <- -1}
  if (sum(abs(Im(I.eigen)))==0) {pd.check <- min(eigen(I.incomp)$values)}
  if (pd.check<0) {
    I.incomp <- nearPD(I.incomp,keepDiag = TRUE)$mat
    print("Estimated Information Matrix not Positive Definite")
    print("Finding nearest Positive Definite Matrix")
  }
  EigenD <- eigen(I.incomp)
  I.eigen <- EigenD$values
  Q <- EigenD$vectors
  CV <- Q%*%diag(1/I.eigen)%*%t(Q)
  #CV2 <- solve(I.incomp)
  return(CV)
}

boot.cov <- function(Y,logT10,N,I,k.in,reps.in,B,alpha.inv) {
  # Important -- the inputs Y and logT10 should still include the "NA" values
  # Not one of the "edited" datasets replacing with zeros

  n <- dim(Y)[1]

  if (missing(B)) {B <- 50}
  if (missing(alpha.inv)) {alpha.inv <- FALSE}
  #  logT10 <- as.matrix(logT10)

  boot.parms <- matrix(0,nrow=B,ncol=(4*I+2))

  for (b in 1:B) {
    Boot.start <- Sys.time()
    index <- sample(1:n,n,replace=T)
    Y.boot <- Y[index,]
    logT10.boot <- logT10[index,]
    MOMboot <- mom(Y.boot,logT10.boot,N,I)
    MCEMboot <- run.mcem(Y.boot,logT10.boot,N,I,k.in,reps.in,ests.in=MOMboot)

    boot.parms[b,] <- c(MCEMboot$a,MCEMboot$b,MCEMboot$alpha,
                        MCEMboot$beta,MCEMboot$vartau,MCEMboot$rho)
    if (b%%5==0) {
      print(b)
      Boot.end <- Sys.time()
      Boot.time <- Boot.end - Boot.start
      print(Boot.time)
    }
  }
  if (alpha.inv==T) {
    boot.parms[,(2*I+1):(3*I)] <- 1/boot.parms[,(2*I+1):(3*I)]
  }
  CV.boot <- cov(boot.parms)
  return(CV.boot)

}

bootmodel.cov <- function(Y,logT10,N,I,parms,k.in,reps.in,B) {
  # Important -- the inputs Y and logT10 should still include the "NA" values
  # Not one of the "edited" datasets replacing with zeros

  if (missing(B)) {B <- 50}

  n <- dim(Y)[1]
  a <- parms[1:I]
  b <- parms[(I+1):(2*I)]
  alpha <- parms[(2*I+1):(3*I)]
  beta <- parms[(3*I+1):(4*I)]
  var_tau <- parms[4*I+1]
  rho <- parms[4*I+2]

  Empty <- is.na(Y)
  boot.parms <- matrix(0,nrow=B,ncol=(4*I+2))
  for (bb in 1:B) {
    boot.sample <- sim_reading_data(a,b,alpha,beta,var_tau,rho,N,I,n,missing_prop=0,max_missing=0)
    Y.boot <- boot.sample$Y
    Y.boot[Empty==1] <- NA
    logT10.boot <- boot.sample$logT10
    logT10.boot[Empty==1] <- NA
    MOMboot <- mom(Y.boot,logT10.boot,N,I)
    MCEMboot <- run.mcem(Y.boot,logT10.boot,N,I,k.in,reps.in,ests.in=MOMboot)
    boot.parms[bb,] <- c(MCEMboot$a,MCEMboot$b,MCEMboot$alpha,
                         MCEMboot$beta,MCEMboot$vartau,MCEMboot$rho)
  }

  CV.boot <- cov(boot.parms)
  return(CV.boot)

}

#' This is run.wcpm function.
#'
#' Update Memo:
#' 04/29/2021 Modified the wcpm function
#'            based on Nelis's updated.
#' 06/01/2021 Modified based on Nelis's MAP function.
#' 06/01/2021 Modified based on Sarunya's BiEAP function.
#' 06/20/2021 Modified based on Nelis's updated for MLE and EAP
#' 06/21/2021 Modified a bug of MAP function.
#' 07/12/2021 Modified est.eqs function based on Nelis's code.
#' 07/13/2021 Added Map function for bootstrap.
#' 07/30/2021 Modified wcpm function based on Sarunya's update
#'
#' @param object = MCEM object, if not be given will occur error and stop running
#' @param stu.data = student response passage data
#' @param pass.data = estimate parameters data
#' @param cases = student season id vector
#' @param est = estimator, c("mle", "map", "eap", "mcmc"), default "map"
#' @param perfect.cases = perfect accurate case
#' @param lo = default -4
#' @param hi = default 4
#' @param q  = default 100
#' @param kappa  = default 1
#' @param external = if not NULL, will use not student read passages for estimating
#'
#' @import rootSolve
#' @import doParallel
#' @import parallel
#' @import foreach
#' @import tidyverse
#' @import dplyr
#' @import MultiGHQuad
#'
#' @return wcpm list
run.wcpm <- function(object, stu.data, pass.data, cases, perfect.cases, est="map", lo = -4, hi = 4, q = 100, kappa = 1, external=NULL) {
  # loading logger
  log.initiating()
  flog.info("Begin wcpm process", name = "orfrlog")

  # Check MCEM object
  if (class(object)[1] == "mcem") {
    MCEM <- object
  } else { # if no MCEM object stop running
    flog.info("Missed MCEM object, end wcpm process", name = "orfrlog")
    return(NULL)
  }


  # get estimator type
  Estimator <- est
  flog.info(paste(paste("Output", est),"WCPM score"), name = "orfrlog")

  est.theta.tau <- function(stu.data, pass.data, case, lo = -4, hi = 4, q = 100, external=NULL) {

    case_split <- unlist(str_split(case, "_"))
    stu.dat01 <- stu.data %>% filter(stu.data$student.id==case_split[1], stu.data$occasion==case_split[2])
    if (nrow(stu.dat01) == 0) {
      # flog.info(paste("No data for:", case), name = "orfrlog")
      return(NULL)
    }
    pass.read <- stu.dat01 %>% select(passage.id)

    # passage.id should be included in MCEM object
    pass.dat01 <- pass.data %>% semi_join(pass.read, by = "passage.id")
    n.pass <- nrow(pass.dat01)

    numwords.total <- stu.dat01 %>% select(numwords.p) %>% c() %>% unlist() %>% sum()
    grade <- stu.dat01 %>% select(grade) %>% c() %>% unlist %>% unique()

    wrc <- stu.dat01 %>% select(wrc) %>% c() %>% unlist()
    lgsec <- stu.dat01 %>% select(lgsec) %>% c() %>% unlist()
    numwords.p <- stu.dat01 %>% select(numwords.p) %>% c() %>% unlist()
    lgsec10 <- lgsec-log(numwords.p) + log(10)

    a.par <- pass.dat01 %>% select(a) %>% c() %>% unlist()
    b.par <- pass.dat01 %>% select(b) %>% c() %>% unlist()
    alpha.par <- pass.dat01 %>% select(alpha) %>% c() %>% unlist()
    beta.par <- pass.dat01 %>% select(beta) %>% c() %>% unlist()

    if (!is.null(external))  { # When external passages

      # get a, b, alpha, beta from MCEM with specific passage.id
      a.par.external <- pass.data %>% filter(passage.id %in% external) %>% select(a) %>% c() %>% unlist()
      b.par.external <- pass.data %>% filter(passage.id %in% external) %>% select(b) %>% c() %>% unlist()
      alpha.par.external <- pass.data %>% filter(passage.id %in% external) %>% select(alpha) %>% c() %>% unlist()
      beta.par.external <- pass.data %>% filter(passage.id %in% external) %>% select(beta) %>% c() %>% unlist()
      numwords.p.external <- pass.data %>% filter(passage.id %in% external) %>% select(numwords.p) %>% c() %>% unlist()
    }

    # Compute n.pass.wcpm and numwords.total.wcpm
    if (is.null(external)) {
      n.pass.wcpm <- n.pass
      numwords.total.wcpm <- sum(numwords.p)
    } else {
      n.pass.wcpm <- length(external)
      numwords.total.wcpm <- sum(numwords.p.external)
    }

    # Using MCEM to calculate rho and vartau
    rho <- mean(MCEM$hyper.param$rho)
    vartau <- mean(MCEM$hyper.param$vartau)

    # Compute observed accuracy (wrc), speed (secs), and fluency (wcpm).
    wrc.obs <- stu.dat01 %>% select(wrc) %>% sum()
    secs.obs <- stu.dat01 %>% select(sec) %>% sum()
    wcpm.obs <- wrc.obs/secs.obs*60

    mod.pd1 <- function(theta) {
      eta <- a.par*theta - b.par
      #term1 <- sum(a.par*wrc*dnorm(eta)/pnorm(eta))
      #term2 <- sum(a.par*(numwords.p-wrc)*dnorm(eta)/(1-pnorm(eta)))
      term1 <- sum(a.par*wrc*exp(dnorm(eta,log = TRUE)-pnorm(eta,log.p = TRUE)))
      term2 <- sum(a.par*(numwords.p-wrc)*exp(dnorm(eta,log = TRUE)-pnorm(eta, lower.tail = FALSE, log.p = TRUE)))
      pd1 <- term1 - term2
    }

    # Following logic will not work on the perfect accurate cases
    # So check it first
    # Initiating the variables
    theta.mle <- Inf # for perfect case
    eta <- NA
    se.theta.mle <- NA
    wrc.mle0 <- NA
    wcpm.mle0 <- NA
    k.theta0 <- NA
    se.wcpm.mle0 <- NA
    secs.mle0 <- NA

    # MLE for tau & theta
    tau.mle <- sum(alpha.par^2*(beta.par - lgsec10))/sum(alpha.par^2)
    se.tau <- sum(alpha.par^2)^(-0.5)

    # only for non-perfect season case
    if (!(case %in% perfect.cases$perfect.cases)) {
      theta.mle <- uniroot(mod.pd1, c(-12, 12))$root
      eta <- a.par*theta.mle - b.par
      #se.theta.mle <- sum((a.par*numwords.p*dnorm(eta))/(pnorm(eta)*(1-pnorm(eta))))^(-0.5)
      I.theta <- sum(a.par^2*numwords.p*dnorm(eta)^2/(pnorm(eta)*(1-pnorm(eta))))
      se.theta.mle <- 1/sqrt(I.theta)

      if (is.null(external)) { #internal
        #secs.mle0 <- sum(exp(beta.par + log(numwords.p) - tau.mle + ((1/alpha.par)^2)/2))
        secs.mle0 <- sum(exp(beta.par - log(10) + log(numwords.p) - tau.mle + ((1/alpha.par)^2)/2))
        wrc.mle0 <- sum(numwords.p*pnorm(a.par*theta.mle - b.par))
        wcpm.mle0 <- wrc.mle0/secs.mle0*60
        k.theta0 <- sum(a.par*numwords.p*dnorm( a.par*theta.mle - b.par ))/sum(numwords.p*pnorm( a.par*theta.mle - b.par ))
        se.wcpm.mle0 <- wcpm.mle0*(k.theta0^2*se.theta.mle^2 + se.tau^2)^0.5
      } else {
        # if external, will calculate with external a, b, alpha, and beta
        secs.mle0 <- sum(exp(beta.par.external - log(10) + log(numwords.p.external) - tau.mle + ((1/alpha.par.external)^2)/2))
        wrc.mle0 <- sum(numwords.p.external*pnorm(a.par.external*theta.mle - b.par.external))
        wcpm.mle0 <- wrc.mle0/secs.mle0*60
        k.theta0 <- sum(a.par.external*numwords.p.external*dnorm( a.par.external*theta.mle - b.par.external ))/sum(numwords.p.external*pnorm( a.par.external*theta.mle - b.par.external ))
        se.wcpm.mle0 <- wcpm.mle0*(k.theta0^2*se.theta.mle^2 + se.tau^2)^0.5

      }
    }

    if (Estimator == "mle") {
      # flog.info(paste(paste("Output", est),"WCPM score"), name = "orfrlog")

      out <- tibble(student.id=case_split[1], occasion=case_split[2], grade=grade,
                    n.pass=n.pass, numwords.total=numwords.total,
                    wrc.obs, secs.obs, wcpm.obs,
                    tau.mle,
                    theta.mle,
                    se.tau.mle=se.tau,
                    se.theta.mle,
                    wrc.mle=wrc.mle0, secs.mle=secs.mle0,
                    n.pass.wcpm,numwords.total.wcpm,
                    wcpm.mle=wcpm.mle0, se.wcpm.mle=se.wcpm.mle0)
      return(out)
    } else if (Estimator == "map") {
      # flog.info(paste(paste("Output", est),"WCPM score"), name = "orfrlog")
      # Add map estimation function
      est.eqs <- function(latent.parms) {
        theta <- latent.parms[1]
        tau <- latent.parms[2]
        eta <- a.par*theta - b.par

        ee1 <- -1/(kappa^2*(1-rho^2))*(theta-rho/sqrt(vartau)*tau) +
          sum(a.par*wrc*exp(dnorm(eta,log = TRUE)-pnorm(eta,log.p = TRUE))) -
          sum(a.par*(numwords.p-wrc)*exp(dnorm(eta,log = TRUE)-pnorm(eta, lower.tail = FALSE, log.p = TRUE)))
        # Modified the bug here
        # ee2 <- -1/(kappa^2*(1-rho^2))*(theta/vartau-rho/sqrt(vartau)*tau) -
        #   sum(alpha.par*(lgsec10 - beta.par + tau))
        ee2 <- -1/(kappa^2*(1-rho^2))*(tau/vartau-rho/sqrt(vartau)*theta) -
          sum(alpha.par^2*(lgsec10 - beta.par + tau))
        ee <- c(ee1,ee2)
        return(ee)
      }

      # MAP estimation: Updated 6/7/2022
      ests.map <- NA
      in.vals <- c(max(-5,min(5,theta.mle)),max(-5*sqrt(vartau),min(5*sqrt(vartau),tau.mle)))
      ests.map <- rootSolve::multiroot(est.eqs, in.vals)$root
      # MAP Standard Errors for tau and theta
      se.tau.map <- sum(alpha.par^2)^(-0.5)
      eta <- a.par*ests.map[1] - b.par
      I.theta <- sum(a.par^2*numwords.p*dnorm(eta)^2/(pnorm(eta)*(1-pnorm(eta))))
      se.theta.map <- 1/sqrt(I.theta)
      # MAP WCPM score
      if (is.null(external)) { #internal
        wrc.map <- sum(numwords.p*pnorm(a.par*ests.map[1] - b.par))
        secs.map <- sum(exp(beta.par - log(10) + log(numwords.p) - ests.map[2] + ((1/alpha.par)^2)/2))
        wcpm.map <- wrc.map/secs.map*60
        k.theta.map <- sum(a.par*numwords.p*dnorm( a.par*ests.map[1] - b.par ))/sum(numwords.p*pnorm( a.par*ests.map[1] - b.par ))
        se.wcpm.map <- wcpm.map*(k.theta.map^2*se.theta.map^2 + se.tau^2)^0.5
      } else {
        # if external, will calculate with external a, b, alpha, and beta
        wrc.map <- sum(numwords.p.external*pnorm(a.par.external*ests.map[1] - b.par.external))
        secs.map <- sum(exp(beta.par.external - log(10) + log(numwords.p.external) - ests.map[2] + ((1/alpha.par.external)^2)/2))
        wcpm.map <- wrc.map/secs.map*60
        k.theta.map <- sum(a.par.external*numwords.p.external*dnorm( a.par.external*ests.map[1] - b.par.external ))/sum(numwords.p.external*pnorm( a.par.external*ests.map[1] - b.par.external ))
        se.wcpm.map <- wcpm.map*(k.theta.map^2*se.theta.map^2 + se.tau^2)^0.5
      }
      # End of MAP estimation

      out <- tibble(student.id=case_split[1], occasion=case_split[2], grade=grade,
                    n.pass=n.pass, numwords.total=numwords.total,
                    wrc.obs, secs.obs, wcpm.obs,
                    tau.map=ests.map[2],
                    theta.map=ests.map[1],
                    # add these two columns, similar to EAP output
                    se.tau.map=se.tau.map,
                    se.theta.map=se.theta.map,
                    wrc.map, secs.map,
                    n.pass.wcpm,numwords.total.wcpm,
                    wcpm.map, se.wcpm.map)
      return(out)
    } else if (Estimator == "eap") {
      # flog.info(paste(paste("Output", est),"WCPM score"), name = "orfrlog")
      # EAP for theta
      Q <- seq(lo, hi, length = q)
      width <- (Q[2] - Q[1])/2
      Qh <- Q + width
      cw <- pnorm(Qh)
      cw2 <- c(0, cw[ - q])
      WQ <- cw - cw2

      LQ <- rep(0, q)
      for(i in 1:q) {
        eta <- a.par*Q[i] - b.par
        binom.lik <- dbinom(wrc, numwords.p, pnorm(eta), log = T)
        LQk <- exp(sum(binom.lik))
        LQ[i] <- LQk
      }

      theta.eap <- sum(Q * LQ * WQ)/sum(LQ * WQ)
      se.theta.eap <- sqrt(sum((Q - theta.eap)^2 * LQ * WQ)/sum(LQ * WQ))
      #se.theta.eap <- sum((Q - theta.eap)^2 * LQ * WQ)/sum(LQ * WQ)

      # Bivariate EAP for theta and tau
      cov <- rho*sqrt(vartau)
      prior <- list(mu = c(0,0), Sigma = matrix(c(1,cov,cov,vartau),2,2))
      #grid <- init.quad(Q = 2, prior, ip = 100, prune = TRUE)
      grid <- MultiGHQuad::init.quad(Q = 2, prior, ip = 500, prune = F)

      loglik <- function(z) {
        theta <- z[1]
        tau <- z[2]
        loglik.bi <- sum(dbinom(wrc, numwords.p, pnorm((a.par*theta)-b.par), log = T)) +
          sum(dnorm(lgsec10, beta.par-tau, 1/alpha.par, log = T))
      }

      ests.quad <- MultiGHQuad::eval.quad(loglik, grid)
      varmat <- attr(ests.quad, "variance")
      se.quad <- c(sqrt(varmat[1,1]), sqrt(varmat[2,2]))
      # QUAD WCPM score
      if (is.null(external)) { #internal
        wrc.quad <- sum(numwords.p*pnorm(a.par*ests.quad[1] - b.par))
        secs.quad <- sum(exp(beta.par - log(10) + log(numwords.p) - ests.quad[2] + ((1/alpha.par)^2)/2))
        wcpm.quad <- wrc.quad/secs.quad*60
        k.theta.quad <- sum(a.par*numwords.p*dnorm( a.par*ests.quad[1] - b.par ))/sum(numwords.p*pnorm( a.par*ests.quad[1] - b.par ))
        se.wcpm.quad <- wcpm.quad*(k.theta.quad^2*se.quad[1]^2 + se.quad[2]^2)^0.5
      } else {
        # if external, will calculate with external a, b, alpha, and beta
        wrc.quad <- sum(numwords.p.external*pnorm(a.par.external*ests.quad[1] - b.par.external))
        secs.quad <- sum(exp(beta.par.external - log(10) + log(numwords.p.external) - ests.quad[2] + ((1/alpha.par.external)^2)/2))
        wcpm.quad <- wrc.quad/secs.quad*60
        k.theta.quad <- sum(a.par.external*numwords.p.external*dnorm( a.par.external*ests.quad[1] - b.par.external ))/sum(numwords.p.external*pnorm( a.par.external*ests.quad[1] - b.par.external ))
        se.wcpm.quad <- wcpm.quad*(k.theta.quad^2*se.quad[1]^2 + se.quad[2]^2)^0.5
      }
      # End of BiEAP

      out <- tibble(student.id=case_split[1], occasion=case_split[2], grade=grade,
                    n.pass=n.pass, numwords.total=numwords.total,
                    wrc.obs, secs.obs, wcpm.obs,
                    tau.eap = ests.quad[2],
                    theta.eap = ests.quad[1],
                    se.tau.eap = se.quad[2],
                    se.theta.eap = se.quad[1],
                    wrc.eap = wrc.quad,
                    secs.eap = secs.quad,
                    n.pass.wcpm,numwords.total.wcpm,
                    wcpm.eap = wcpm.quad,
                    se.wcpm.eap = se.wcpm.quad)
      return(out)
    }

  }

  numCores <- detectCores() - 1


  cl <- makeCluster(numCores)
  registerDoParallel(cl)

  seq_id_all <- nrow(cases)

  theta.tau <- foreach(i=1:seq_id_all, .combine = 'rbind', .packages = c("tidyverse", "MultiGHQuad")) %dopar%
    { est.theta.tau(stu.data,
                    pass.data,
                    cases$cases[i],
                    lo, hi, q, external=external)
    }
  if (length(theta.tau[1]) != 0) {
    class(theta.tau) <- "wcpm" #define wcpm class
  }

  on.exit(stopCluster(cl))

  flog.info("End wcpm process", name = "orfrlog")
  return(invisible(theta.tau))
}
