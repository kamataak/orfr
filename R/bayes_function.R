##################################################################################################
##### THIS IS THE FUNCTION TO ESTIMATE PASSAGE PARAMETERS AND PERSON-LEVEL HYPERPARAMETERS########
##################################################################################################
#' Bayes function when running mcem with mcmc setting
#'
#' @param person.data - student reading data
#' @param person.id The column name in the data that represents the unique individual identifier.
#' @param task.id The column name in the data that represents the unique task identifier.
#' @param max.counts The column name in the data that represents the number of words in a task.
#' @param obs.counts The column name in the data that represents the words read correctly for each case.
#' @param time The column name in the data that represents the time, in seconds, for each case.
#' @param parallel parallel=T, #logical, run in parallel? "T" or "F"
#' @param n.chains int., number of the chains
#' @param thin int, thinning interval, a.k.a, period of saving samples
#' @param iter int., number of the iterations after the burn-in period
#' @param burn int., number of the burn-in iterations
#'
#'
#' @return list
bayes <- function(
    person.data=NA,   # data frame, long format, required columns: person.id, task.id, max.counts, obs.counts, time
    person.id = "",
    task.id = "",
    max.counts = "",
    obs.counts = "",
    time = "",
    parallel=T, #logical, run in parallel? "T" or "F"
    n.chains=NA, # pos. int., number of the chains
    thin=1, #pos. int, thinning interval, a.k.a, period of saving samples
    iter=NA,  # pos. int., number of the iterations after the burn-in period
    burn=NA  # pos. int., number of the burn-in iterations
)
{

  person.data <- person.data[,c(person.id, task.id, max.counts, obs.counts, time)]
  colnames(person.data) <- c("person.id", "task.id", "max.counts", "obs.counts", "time")

  #Extract the sub-components of the data for analyses
  time.data <- person.data %>%
    select(person.id, task.id, time) %>%
    pivot_wider(names_from = task.id, values_from = time) %>%
    column_to_rownames("person.id") %>%
    select(sort(colnames(.))) %>%
    as.matrix()

  count.data <- person.data %>%
    select(person.id, task.id, obs.counts) %>%
    pivot_wider(names_from = task.id, values_from = obs.counts) %>%
    column_to_rownames("person.id") %>%
    select(sort(colnames(.))) %>%
    as.matrix()

  n.words <- person.data %>%
    select(task.id, max.counts) %>%
    arrange(task.id) %>%
    distinct() %>%
    deframe()


  #Identify number of the examinees and passages from time data matrix
  #And bundle the data as a list (note that time matrix converted to log seconds!)
  J <- nrow(time.data)
  I <- ncol(time.data)

  data.list <- list(J=J, I=I, tim=log(time.data), res=count.data, nw=n.words)

  #Estimate only main parameters
  param_est <- c("a", "b", "alpha", "beta", "vartau", "rho")

  #Specify number of the chains
  if(is.na(n.chains)){
    n.chains <- min(max(4), detectCores()-1)
  }else{
    n.chains <- min(max(4), n.chains)
  }

  #Generate initial values based on the number of chains
  inits <- vector(mode = "list", length = n.chains)

  for(i in 1:n.chains){
    a=runif(I, 0, 1)
    b=runif(I, -4, -1)
    alpha=runif(I, 3, 10)
    beta_raw=runif(I, 1, 2)+log(n.words/10) #beta_raw = beta_10 + log(nwords/10)
    cvr=runif(1, 0, 1)
    ptau=1/(runif(1, 0, 1))

    #Select different RNG names and seeds for each parallel chain as recommended
    rng_names <- c("base::Wichmann-Hill", "base::Marsaglia-Multicarry", "base::Super-Duper", "base::Mersenne-Twister")

    if(n.chains <= length(rng_names)){
      rng_name_sel <- rng_names[i]
    }else{
      rng_names_rep <- rep(rng_names, n.chains)
      rng_name_sel <- rng_names_rep[i]
    }

    gen_init <- list(a=a, b=b, alpha=alpha, beta_raw=beta_raw, cvr=cvr, ptau=ptau,
                     .RNG.name=rng_name_sel, .RNG.seed=i)
    inits[[i]] <- gen_init
  }


  #Now, check if the data have any missing values.
  #If so, we will use JAGS. If not, we will use STAN.
  time.mis <- T %in% is.na(time.data)
  count.mis <- T %in% is.na(count.data)

  if(time.mis==T | count.mis==T){
    bayes.soft="jags"
    cat("\n \n ==== Estimation will be done with JAGS ==== \n \n")
  }else{
    bayes.soft="stan"
    cat("\n \n ==== Estimation will be done with STAN ==== \n \n")
  }


  if(bayes.soft=="jags"){

    runjags::runjags.options(force.summary=T)

    # -------------------------------------------------------- JAGS syntax
    jags.syntax <- "
    model{
# J students and I passages
    for (j in 1:J) {
      for (i in 1:I) {
        res[j,i] ~ dbin(p[j,i], nw[i])
        probit(p[j,i]) <- a[i] * theta[j] - b[i]
        tim[j,i] ~ dnorm(mu[j,i], prec.t[i])
        mu[j,i] <- beta_raw[i] - tau[j]
      }
      theta[j] ~ dnorm(0,1)
      tau[j] ~ dnorm(mtau[j], ptau)
      mtau[j] <- cvr * theta[j]
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
    vartau <- vtau + (pow(cvr, 2))
    cvr ~ dnorm(0, 0.01)
    rho <- cvr/sqrt(vartau)
}
"

  # ---------------------------------------------------------------------- JAGS Model Syntax

  #Check if the user specified parallel simulation or not
  if(isTRUE(parallel)){
    jags_meth <- "parallel"
  } else{
    jags_meth <- "rjags"
  }

  jags_out <- runjags::autorun.jags(
    model = jags.syntax,
    monitor = param_est,
    data = data.list,
    n.chains = n.chains,
    inits = inits,
    thin = thin,
    method = jags_meth,
    modules = "glm"
  )

  par_est <- jags_out$summaries %>%
    as.data.frame() %>%
    rownames_to_column(var = "Parameter") %>%
    select(Parameter, Mean, SD, Lower95, Upper95, ESS=SSeff, BGR=psrf)


  }else if(bayes.soft=="stan"){

    rstan::rstan_options(auto_write = TRUE)

    #Transpose data matrices
    time.data <- time.data %>% t()
    count.data <- count.data %>% t()

    data.list <- list(J=J, I=I, tim=log(time.data), res=count.data, nw=n.words)



    # ------------------------------------------------------------- STAN Syntax
    stan.syntax <- "
data{
  int <lower=0> J; //number of individuals
  int <lower=0> I; //number of passages
  int <lower=0> res[I,J]; //array of counts
  real tim[I,J]; //array of times
  int <lower=0> nw[I]; //vector of number of the words per passage
}

parameters{
  vector <lower=0> [I] alpha; //time discrimintion
  vector[I] beta_raw; //time intensity
  vector <lower=0> [I] a; //accuracy discrimination
  vector[I] b; //accuracy difficulty (threshold style)

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
  real <lower=0> vartau; //Actaul variance of tau in the MVN distribution
  real rho; //Correlation between theta and tau
  vector[I] beta; //time intensity per 10 words;
  vartau=square(stau) + square(cvr);
  rho=cvr/sqrt(vartau);

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

  tau ~ normal(cvr * theta, stau);


// Likelihood
  for(i in 1:I){
  res[i] ~ binomial(nw[i], Phi(a[i] * theta - b[i]));
  tim[i] ~ normal(beta_raw[i] - tau, 1/alpha[i]);
              }
}

"
# -------------------------------------------------------------------------------------------------------

#Specify the parallel running or not!
if(isTRUE(parallel)){
  n.cores <- min(max(4), detectCores()-1)
} else{
  n.cores <- 1
}

stan_out <- rstan::stan(model_code = stan.syntax,
                        pars = param_est,
                        data = data.list,
                        chains = n.chains,
                        warmup  = 2e3, #for now, keep those values as default for stan
                        iter = 1e4,
                        thin = thin,
                        cores = n.cores,
                        init = inits,
                        control = list(adapt_delta = 0.99)
)

par_est <- summary(stan_out)$summary %>%
  as.data.frame() %>%
  rownames_to_column(var = "Parameter") %>%
  select(Parameter, Mean=mean, SD=sd, Lower95=`2.5%`, Upper95=`97.5%`, ESS=n_eff, BGR=Rhat) %>%
  mutate(ESS=round(ESS)) %>%
  filter(Parameter!="lp__")

  }

#Create an output as the same structure as mcem function
par_est_list <- list(
  task.param=tibble(
    a=par_est$Mean[grep(pattern = "^a\\[", x = par_est$Parameter)],
    b=par_est$Mean[grep(pattern = "^b\\[", x = par_est$Parameter)],
    alpha=par_est$Mean[grep(pattern = "^alpha\\[", x = par_est$Parameter)],
    beta=par_est$Mean[grep(pattern = "^beta\\[", x = par_est$Parameter)],
    se.a=par_est$SD[grep(pattern = "^a\\[", x = par_est$Parameter)],
    se.b=par_est$SD[grep(pattern = "^b\\[", x = par_est$Parameter)],
    se.alpha=par_est$SD[grep(pattern = "^alpha\\[", x = par_est$Parameter)],
    se.beta=par_est$SD[grep(pattern = "^beta\\[", x = par_est$Parameter)],
    task.id=names(n.words),
    max.counts=n.words),

  hyper.param=tibble(
    vartau=par_est$Mean[grep(pattern = "vartau", x = par_est$Parameter)],
    rho=par_est$Mean[grep(pattern = "rho", x = par_est$Parameter)],
    se.vartau=par_est$SD[grep(pattern = "vartau", x = par_est$Parameter)],
    se.rho=par_est$SD[grep(pattern = "rho", x = par_est$Parameter)]))

par_est_list

}
