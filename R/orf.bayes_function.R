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
#' This is orf.bayes function.
#' Estimates Oral Reading Fluency (ORF) model parameters
#' with the fully Bayesian approach.
#'
#' @param time.data   #matrix, examinees X passages, raw time scale (in seconds)
#' @param count.data   #matrix, examinees X passages, number of the words read correctly per passage
#' @param numwords.pass      #int. vec., vector of number of the words per passage!
#' @param param     #chr. vec., which parameters to estimate. Default is only passage parameters!
#' @param bayes.soft #chr., which software to use "jags" or "stan"?
#' @param parallel #logical, run in parallel? "T" or "F"
#' @param n.chains #pos. int., number of the chains
#' @param iter #pos. int., number of the iterations after the burn-in period
#' @param burn #pos. int., number of the burn-in iterations
#' @param thin #os. int, thinning interval, a.k.a, period of saving samples
#'
#' @import parallel
#' @import tidyverse
#' @import dplyr
#'
#' @return stan or jags data set ()
#' @export
orf.bayes <- function(
  time.data=NA,  # matrix, examinees X passages, raw time scale (in seconds)
  count.data=NA, # matrix, examinees X passages, number of the words read correctly per passage
  numwords.pass=NA, #int. vec., vector of number of the words per passage!
  param=NA, # chr. vec., which parameters to estimate. Default is only passage parameters!
  bayes.soft="jags", # chr., which software to use "jags" or "stan"?
  parallel=T, #logical, run in parallel? "T" or "F"
  n.chains=1, # pos. int., number of the chains
  iter=3e4,  # pos. int., number of the iterations after the burn-in period
  burn=1e4,  # pos. int., number of the burn-in iterations
  thin=1 #pos. int, thinning interval, a.k.a, period of saving samples
  )
  {

  #Identify number of the examinees and passages from time data matrix
  #And bundle the data as a list
  J <- nrow(time.data)
  I <- ncol(time.data)
  data.list <- list(J=J, I=I, tim=log(time.data), res=count.data, nw=numwords.pass)

  #Create a full name of the parameters as the default parameters to be estimated!
  if(is.na(param)){
  param_est <- c("a", "b", "alpha", "beta")
  }else{
    param_est <- param
  }

  #Now, create a syntax file for the model based on the selected software
  #And follow software-based specifications to run the analyses!

  if(bayes.soft=="jags"){
    sink("model.txt") # ----------------------------------------------------- JAGS Model Syntax
    cat("
    model{
# J students and I passages
    for (j in 1:J) {
      for (i in 1:I) {
        res[j,i] ~ dbin(p[j,i], nw[i])
        cnt_ex[j,i] <- p[j,i] * nw[i]
        probit(p[j,i]) <- a[i] * (theta[j] - b[i])
        tim[j,i] ~ dnorm(mu[j,i], prec.t[i])
        mu[j,i] <- beta_raw[i] - tau[j]
        tim_ex[j,i] <- exp(mu[j,i] + 0.5 * 1/(pow(alpha[i],2)))
      }
      theta[j] ~ dnorm(0,1)
      tau[j] ~ dnorm(mtau[j], ptau)
      mtau[j] <- cvr * theta[j]
      exp_cnt[j] <- sum(cnt_ex[j,])
      exp_tim[j] <- sum(tim_ex[j,])
      orf[j] <- exp_cnt[j]/exp_tim[j]*60
    }
# Priors for passage parameters
    for(i in 1:I) {
      prec.t[i] <- pow(alpha[i], 2)
      alpha[i]  ~ dnorm(0, 0.01) I(0,)
      beta_raw[i] ~ dnorm(0, 0.01)
      a[i] ~ dnorm(0, 0.01) I(0,)
      b[i] ~ dnorm(0, 0.01)
      beta[i] <- beta_raw[i] - log(nw[i]/10)
    }
# Priors for person parameters
    ptau ~ dgamma(0.01, 0.01)
    vtau <- 1/ptau
    var_tau <- vtau + (pow(cvr, 2))
    cvr ~ dnorm(0, 0.01)
    crl <- cvr/sqrt(var_tau)
}
        ")

sink() # ---------------------------------------------------------------------- JAGS Model Syntax

if(isTRUE(parallel)){
  jags_meth <- "parallel"
} else{
  jags_meth <- "rjags"
}

jags_out <- runjags::run.jags(
      model = "model.txt",
      monitor = param_est,
      data = data.list,
      n.chains = n.chains,
      burnin = burn,
      sample = iter,
      thin = thin,
      method = jags_meth
    )

    jags_est <- runjags::add.summary(jags_out)$summaries %>%
      as.data.frame() %>%
      tibble::rownames_to_column(var = "Parameter") %>%
      dplyr::select(Parameter, Mean, SD, Lower95, Upper95, ESS=SSeff, BGR=psrf)


  }else if(bayes.soft=="stan"){
    sink("model.stan")  # -----------------------------------------------------------------
    cat("
data{
  int <lower=0> J; //number of individuals
  int <lower=0> I; //number of passages
  int <lower=0> res[J,I]; //array of counts
  real tim[J,I]; //array of times
  int <lower=0> nw[I]; //vector of number of the words per passage
}

parameters{
  vector <lower=0> [I] alpha; //time discrimintion
  vector[I] beta_raw; //time intensity
  vector <lower=0> [I] a; //accuracy discrimination
  vector[I] b; //accuracy difficulty

  real <lower=0> sigma_alpha; //sd of alpha's lognormal prior
  real mu_beta; //mean of beta's normal prior
  real <lower=0> sigma_beta; //sd of beta's normal prior
  real <lower=0> sigma_a; //sd of a's lognormal prior
  real mu_b; //mean of b's normal prior
  real <lower=0> sigma_b; //sd of b's normal prior

  real <lower=0> stau; //sd of tau[j]'s normal distribution in the model
  real cvr; //covariance between theta and tau
  real <lower=0> sigma_cvr; //sd of cvr's normal prior

  vector[J] theta; //accuracy ability
  vector[J] tau; //speed ability
}

transformed parameters{
  real <lower=0> var_tau; //Actaul variance of tau in the MVN distribution
  real crl; //Correlation between theta and tau
  vector[I] beta; //time intensity per 10 words;
  var_tau=square(stau) + square(cvr);
  crl=cvr/sqrt(var_tau);
  for(i in 1:I){
  beta[i]=beta_raw[i] - log(nw[i]/10.0);
  }

}

model{
  // Priors
  alpha  ~ lognormal(0, sigma_alpha);
  sigma_alpha ~ cauchy(0,5);
  beta_raw ~ normal(mu_beta, sigma_beta);
  mu_beta ~ normal(0,5);
  sigma_beta ~ cauchy(0,5);
  a ~ lognormal(0, sigma_a);
  sigma_a ~ cauchy(0,5);
  b ~ normal(mu_b, sigma_b);
  mu_b ~ normal(0,5);
  sigma_b ~ cauchy(0,5);

  stau ~ cauchy(0,5);
  cvr ~ normal(0, sigma_cvr);
  sigma_cvr ~ cauchy(0,5);
  theta ~ normal(0, 1);


// Likelihood
for(j in 1:J){vector[J] mtau;
              mtau[j]=cvr * theta[j];
              tau[j] ~ normal(mtau[j], stau);
  for(i in 1:I){real p;
              p=Phi(a[i] * (theta[j] - b[i])); //used probit link to match Nelis's approach!
              res[j,i] ~ binomial(nw[i], p); //nw is the vector of number of the words per passage
              tim[j,i] ~ normal(beta_raw[i] - tau[j], 1/alpha[i]); //This is the normal distribution of log-times!
              }}
}


//Model-based WCPM
generated quantities{
  //Expeced time matrix
  real tim_ex[J,I];
  //Expected count matrix
  real <lower=0> cnt_ex[J,I];
  //Model-based WCPM
  vector <lower=0> [J] orf;


  for(j in 1:J){
    for(i in 1:I){real p; real mu;
      p=Phi(a[i] * (theta[j] - b[i]));
      cnt_ex[j,i]=p*nw[i];
      mu=beta_raw[i] - tau[j];
      tim_ex[j,i]=exp(mu + 0.5 * 1/(square(alpha[i])));
    }
    orf[j]=sum(cnt_ex[j,])/sum(tim_ex[j,])*60;
  }
}
")

sink() # -------------------------------------------------------------------------------------------------------

#Specify the parallel running or not!
if(isTRUE(parallel)){
  n.cores <- detectCores()-1
} else{
  n.cores <- 1
}

stan_out <- rstan::stan(file = "model.stan",
  pars = param_est,
  data = data.list,
  chains = n.chains,
  warmup  = burn,
  iter = iter,
  thin = thin,
  cores = n.cores
)

stan_est <- summary(stan_out)$summary %>%
  as.data.frame() %>%
  tibble::rownames_to_column(var = "Parameter") %>%
  dplyr::select(Parameter, Mean=mean, SD=sd, Lower95=`2.5%`, Upper95=`97.5%`, ESS=n_eff, BGR=Rhat) %>%
  dplyr::filter(Parameter!="lp__")
  }
}

runBayes <- function(Y,logT,N,I,K.in,reps.in,ests.in,data_check,est) {

  #--------------prepare data--------------#
  set.seed(1234)
  theta <-


  per_par <- mvtnorm::rmvnorm(n = 100, mean = c(0, 0), sigma = matrix(c(1, 0.07, 0.07, 0.03), nc=2))

  theta <- per_par[,1]
  tau <- per_par[,2]

  a_par <- c(0.4258509, 0.4941767, 0.4701327, 0.5048772, 0.4470414, 0.4760288, 0.5547857, 0.5231372, 0.5746109)
  b_par <- c(-3.159198, -2.831126, -2.943872, -2.855551, -2.948745, -2.877003, -2.419438, -2.399501, -2.425473)
  alpha_par <- c(8.710270, 4.758931, 6.478158, 5.063602, 4.408486, 3.935842, 6.277621, 5.858994, 6.820602)
  beta_par <- c(3.921922, 3.661640, 3.910936, 3.359062, 3.314531, 3.341339, 3.438361, 3.412690, 3.492000) #in raw scale, nor per 10 words!

  numwords.pass <- c(88, 69, 87, 50, 50, 49, 52, 51, 54)


  J <- length(theta) # number of the examinees
  I <- length(numwords.pass) #number of the passages

  a_rep <- rep(a_par, each=J)
  b_rep <- rep(b_par, each=J)
  alpha_rep <- rep(alpha_par, each=J)
  beta_rep <- rep(beta_par, each=J)

  theta_rep <- rep(theta, I)
  tau_rep <- rep(tau, I)

  numwords.pass_rep <- rep(numwords.pass, each=J)

  #Generate Time data
  mu_time <- beta_rep - tau_rep
  sd_time <- 1/alpha_rep
  time_data <- rnorm(n=J*I, mean=mu_time, sd=sd_time)
  time_mat <- matrix(time_data, nc=I) %>% exp() #Make it in original scale!

  #Generate count data
  #prb <- 1/(1+exp(-a_rep*(theta_rep-b_rep))) #Logit model
  prb <- pnorm(a_rep*(theta_rep-b_rep)) #Probit model
  respon <- rbinom(n = J*I, size = numwords.pass_rep, prob = prb)
  respon_mat <- matrix(respon, nc=I)
}

# call bayes functions
runBayes <- function(object=NULL,mcem=NULL, wcpm=NULL, cases) {
  dat <- object
  MCEM <- mcem
  WCPM <- wcpm

  a_par <- dat$a
  b_par
  alpha_par
  beta_par

}
