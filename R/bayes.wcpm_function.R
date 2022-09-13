##################################################################################################
################### THIS IS THE FUNCTION TO ESTIMATE MODEL-BASED WCPM PARAMETERS #################
##################################################################################################

#First, check whether the required packages are installed or not
# req_pkg <- c("tidyverse", "runjags", "rstan", "parallel")
# if(!any(req_pkg %in% installed.packages())) install.packages(pkgs=req_pkg)

#Now, load the required packages
# lapply(req_pkg, library, character.only=T)

#Start writing the "orf.bayes2" function
bayes.wcpm <- function(
    calib.data=NA,
    stu.data=NA,
    studentid=NULL,
    passageid = NULL,
    season = NULL,
    grade = NULL,
    numwords.p = NULL,
    wrc = NULL,
    time = NULL,
    cases = NULL,
    external=NULL,
    parallel=T, #logical, run in parallel? "T" or "F"
    n.chains=NA, # pos. int., number of the chains
    iter=NA,  # pos. int., number of the iterations after the burn-in period
    burn=NA,  # pos. int., number of the burn-in iterations
    thin=1 #pos. int, thinning interval, a.k.a, period of saving samples
)
{
  # loading logger
  log.initiating()
  flog.info("Begin wcpm process with mcmc setting", name = "orfrlog")

  #First, check if calib.data exists. If not, print a warning.
  if(class(calib.data)[1] == "mcem"){
    pass.data <- calib.data$pass.param
  }else{
    flog.info("Missing MCEM object, end wcpm process", name = "orfrlog")
    return(NULL)
  }

  #Now, check whether user entered the column names of the data or
  #just used the output of prep data. For former, rename columns of stu.data
  if(is.null(c(studentid, passageid, season, grade, numwords.p, wrc, time))){
    stu.data <- stu.data %>%
      select(-lgsec)
    colnames(stu.data) <- c("studentid", "passageid", "numwords.p", "season", "grade", "wrc", "time")
  }else{
    stu.data <- stu.data[,c(studentid, passageid, season, grade, numwords.p, wrc, time)]
    colnames(stu.data) <- c("studentid", "passageid", "season", "grade", "numwords.p", "wrc", "time")
  }

  #Now, check whether users supplied cases or not. If they didn't, then WCPMs will be
  #estimated for all unique cases appear in the supplied student data!
  if(is.null(cases)){
    stu.data <- stu.data
    #%>% arrange(studentid, season, grade)
  }else{
    stu.data <- stu.data %>%
      #arrange(studentid, season, grade) %>%
      mutate(case_sel=paste(studentid, season, sep = "_")) %>%
      filter(case_sel %in% cases$cases) %>%
      select(-case_sel)
  }


  #Identify parameters of the passages read from the calibrated pool
  stu_id <- stu.data %>%
    select(studentid) %>%
    distinct()

  pas_param_read <- calib.data$pass.param %>%
    filter(passage.id %in% stu.data$passageid) %>%
    arrange(passage.id)

  #Create descriptive part of the output
  desc_out <- stu.data %>%
    rename(occasion=season) %>%
    group_by(studentid, occasion, grade) %>%
    summarise(n.pass=n_distinct(passageid),
              numwords.total=sum(numwords.p),
              wrc.obs=sum(wrc),
              secs.obs=sum(time)) %>%
    ungroup() %>%
    mutate(wcpm.obs=wrc.obs/secs.obs*60)

  desc_out <- stu_id %>%
    left_join(desc_out) %>%
    select(studentid, everything())

  #Now, create the datasets for MCMC based on selected cases and passages
  time.data <- stu.data %>%
    select(studentid, passageid, time) %>%
    pivot_wider(names_from = passageid, values_from = time) %>%
    column_to_rownames("studentid") %>%
    select(sort(colnames(.))) %>%
    as.matrix()

  count.data <- stu.data %>%
    select(studentid, passageid, wrc) %>%
    pivot_wider(names_from = passageid, values_from = wrc) %>%
    column_to_rownames("studentid") %>%
    select(sort(colnames(.))) %>%
    as.matrix()

  n.words <- stu.data %>%
    select(passageid, numwords.p) %>%
    arrange(passageid) %>%
    distinct() %>%
    deframe()


  #Now, check if the user supplied external passages or not. If not, WCPM will be estimated for
  #all passages they read in the stu.dat!
  if(is.null(external)){
    nonmis_ind <- !is.na(time.data)
    pas_est_ind <- apply(nonmis_ind, 2, as.numeric)
    rownames(pas_est_ind) <- rownames(time.data)
    pas_param_est <- pas_param_read

    desc_out <- desc_out %>%
      mutate(n.pass.wcpm=n.pass,
             numwords.total.wcpm=numwords.total)
  }else{
    pas_param_est <- calib.data$pass.param %>%
      filter(passage.id %in% external) %>%
      arrange(passage.id)

    pas_est_ind <- matrix(1, nrow = nrow(time.data),
                          ncol = length(external),
                          dimnames = list(rownames(time.data), pas_param_est$passage.id))
    desc_out <- desc_out %>%
      mutate(n.pass.wcpm=n_distinct(pas_param_est$passage.id),
             numwords.total.wcpm=sum(pas_param_est$numwords.p))
  }


  #Identify number of the examinees and passages from time data matrix
  #And bundle the data as a list (note that time matrix is converted to log seconds!)
  #Also, use known passage parameters as well as person hyper parameters!
  J <- nrow(time.data)
  I <- ncol(time.data)
  K <- nrow(pas_param_est)

  #Estimate count, time and wcpm.
  param_est <- c("exp_cnt", "exp_tim", "wcpm")

  #Specify number of the chains
  if(is.na(n.chains)){
    n.chains <- min(max(4), detectCores()-1)
  }else{
    n.chains <- min(max(4), n.chains)
  }

  #Generate initial values based on the number of chains
  inits <- vector(mode = "list", length = n.chains)

  for(i in 1:n.chains){
    theta=rnorm(J)

    #Select different RNG names and seeds for each parallel chain as recommended
    rng_names <- c("base::Wichmann-Hill", "base::Marsaglia-Multicarry", "base::Super-Duper", "base::Mersenne-Twister")

    if(n.chains <= length(rng_names)){
      rng_name_sel <- rng_names[i]
    }else{
      rng_names_rep <- rep(rng_names, n.chains)
      rng_name_sel <- rng_names_rep[i]
    }

    gen_init <- list(theta=theta,
                     .RNG.name=rng_name_sel, .RNG.seed=i)
    inits[[i]] <- gen_init
  }

  #Now, check if the data have any missing values.
  #If so, we will use JAGS. If not, we will use STAN.
  time.mis <- T %in% is.na(time.data)
  count.mis <- T %in% is.na(count.data)

  if(time.mis==T | count.mis==T){
    bayes.soft="jags"
    cat("==== Running the analyses with JAGS ==== \n \n")
  }else{
    bayes.soft="stan"
    cat("==== Running the analysis with STAN ==== \n \n")
  }

  #Now, create a syntax file for model based on the selected software
  #And follow software-based specifications to run the analyses!

  if(bayes.soft=="jags"){

    runjags.options(force.summary=T)

    #Check passage parameters from mcem if they exist.
    #Otherwise, user will supply the known parameter values!
    data.list <- list(J=J,
                      I=I,
                      K=K,
                      tim=log(time.data),
                      res=count.data,

                      nw_read=n.words,
                      a_read=pas_param_read$a,
                      b_reg_read=pas_param_read$a*-pas_param_read$b,
                      alpha_read=pas_param_read$alpha,
                      beta_raw_read=pas_param_read$beta + log(pas_param_read$numwords.p/10),

                      pas_est_ind=pas_est_ind,
                      nw_est=pas_param_est$numwords.p,
                      a_est=pas_param_est$a,
                      b_reg_est=pas_param_est$a*-pas_param_est$b,
                      alpha_est=pas_param_est$alpha,
                      beta_raw_est=pas_param_est$beta + log(pas_param_est$numwords.p/10),

                      ptau=1/calib.data$hyper.param$vartau,
                      cvr=calib.data$hyper.param$rho*sqrt(calib.data$hyper.param$vartau))

    # -------------------------------------------------------- JAGS syntax
    jags.syntax <- "
    model{
# J students and I passages
    for (j in 1:J) {
      for (i in 1:I) {
        res[j,i] ~ dbin(p[j,i], nw_read[i])
        probit(p[j,i]) <- a_read[i] * theta[j] + b_reg_read[i]
        tim[j,i] ~ dnorm(mu[j,i], pow(alpha_read[i], 2))
        mu[j,i] <- beta_raw_read[i] - tau[j]
      }
      theta[j] ~ dnorm(0,1)
      tau[j] ~ dnorm(mtau[j], ptau)
      mtau[j] <- cvr * theta[j]
    }

# Estimation of WCPM for target passages
    for(j in 1:J){
    for(k in 1:K){
      cnt_ex[j,k] <- (phi(a_est[k] * theta[j] + b_reg_est[k]) * nw_est[k]) * pas_est_ind[j,k]
      tim_ex[j,k] <- (exp(beta_raw_est[k] - tau[j] + 0.5 * 1/(pow(alpha_est[k],2))))* pas_est_ind[j,k]
    }
      exp_cnt[j] <- sum(cnt_ex[j,])
      exp_tim[j] <- sum(tim_ex[j,])
      wcpm[j] <- exp_cnt[j]/exp_tim[j]*60
    }
}
"

  # ------------------------------------------------------------------ JAGS Syntax

  #Check if the user specified parallel simulation or not
  if(isTRUE(parallel)){
    jags_meth <- "parallel"
  } else{
    jags_meth <- "rjags"
  }


  jags_out <- autorun.jags(
    model = jags.syntax,
    monitor = param_est,
    data = data.list,
    n.chains = n.chains,
    thin = thin,
    inits = inits,
    method = jags_meth,
    modules = "glm"
  )

  par_est <- jags_out$summaries %>%
    as.data.frame() %>%
    select(Mean, SD, Lower95, Upper95)

  rm(jags_out)
  gc()

  } else if(bayes.soft=="stan") {


    #Transpose data matrices
    time.data <- time.data %>% t()
    count.data <- count.data %>% t()
    pas_est_ind <- pas_est_ind %>% t()

    data.list <- list(J=J,
                      I=I,
                      K=K,
                      tim=log(time.data),
                      res=count.data,

                      nw_read=n.words,
                      a_read=pas_param_read$a,
                      b_reg_read=pas_param_read$a*-pas_param_read$b,
                      alpha_read=pas_param_read$alpha,
                      beta_raw_read=pas_param_read$beta + log(pas_param_read$numwords.p/10),

                      pas_est_ind=pas_est_ind,
                      nw_est=pas_param_est$numwords.p,
                      a_est=pas_param_est$a,
                      b_reg_est=pas_param_est$a*-pas_param_est$b,
                      alpha_est=pas_param_est$alpha,
                      beta_raw_est=pas_param_est$beta + log(pas_param_est$numwords.p/10),

                      stau=sqrt(calib.data$hyper.param$vartau),
                      cvr=calib.data$hyper.param$rho*sqrt(calib.data$hyper.param$vartau))

    # ------------------------------------------------------------- STAN Syntax
    stan.syntax <- "
data{
// Data
  int <lower=0> J; //number of individuals
  int <lower=0> I; //number of passages read
  int <lower=0> K; //number of passages external
  int <lower=0> res[I,J]; //array of counts
  real tim[I,J]; //array of times
  real pas_est_ind[K, J]; //array of passage indicators
  int <lower=0> nw_read[I]; //vector of number of the words per passage
  int <lower=0> nw_est[K]; //vector of number of the words per passage

// Known Parameters
  real <lower=0> stau; // SD of tau
  real cvr; // Covariance between theta and tau

  vector <lower=0> [I] alpha_read; //time discrimintion
  vector[I] beta_raw_read; //time intensity
  vector <lower=0> [I] a_read; //accuracy discrimination
  vector[I] b_reg_read; //accuracy difficulty (threshold style)

  vector <lower=0> [K] alpha_est; //time discrimintion
  vector[K] beta_raw_est; //time intensity
  vector <lower=0> [K] a_est; //accuracy discrimination
  vector[K] b_reg_est; //accuracy difficulty (threshold style)

}

parameters{
  vector[J] theta; //accuracy ability
  vector[J] tau; //speed ability
}


model{
  // Priors
  theta ~ normal(0, 1);

  // Likelihood
for(i in 1:I){
  res[i] ~ binomial(nw_read[i], Phi(a_read[i] * theta + b_reg_read[i]));
  tim[i] ~ normal(beta_raw_read[i] - tau, 1/alpha_read[i]);
              }
}

// Estimation of model-based WCPM
generated quantities{

  real tim_ex[K,J]; //Expeced time matrix
  real <lower=0> cnt_ex[K,J]; //Expected count matrix
  vector <lower=0> [J] exp_cnt; //Model-based wrc
  vector <lower=0> [J] exp_tim; //Model-based secs
  vector <lower=0> [J] wcpm; //Model-based WCPM


 for(j in 1:J){
    for(k in 1:K){real p; real mu;
      p=Phi(a_est[k] * theta[j] + b_reg_est[k]);
      cnt_ex[k,j]=(p*nw_est[k]) * pas_est_ind[k, j];
      mu=beta_raw_est[k] - tau[j];
      tim_ex[k,j]=(exp(mu + 0.5 * 1/(square(alpha_est[k])))) * pas_est_ind[k, j];
    }
    exp_cnt[j]=sum(cnt_ex[,j]);
    exp_tim[j]=sum(tim_ex[,j]);
    wcpm[j]=exp_cnt[j]/exp_tim[j]*60;
  }
}

"

# ------------------------------------------------------------------- STAN Syntax

    #Specify the parallel running or not!
    if(isTRUE(parallel)){
      n.cores <- min(max(4), detectCores()-1)
    } else{
      n.cores <- 1
    }


    stan_out <- stan(model_code = stan.syntax,
                     pars = param_est,
                     data = data.list,
                     chains = n.chains,
                     warmup  = 1e3, #Keep as the default values for stan.
                     iter = 5e3,
                     thin = thin,
                     cores = n.cores,
                     init = inits,
                     control = list(adapt_delta = 0.99)
    )

    par_est <- summary(stan_out)$summary %>%
      as.data.frame() %>%
      rownames_to_column(var = "Parameter") %>%
      filter(Parameter!="lp__") %>%
      select(Mean=mean, SD=sd, Lower95=`2.5%`, Upper95=`97.5%`)
  }

  par_est <- par_est %>%
    mutate(studentid=rep(desc_out$studentid, length(param_est)),
           occasion=rep(desc_out$occasion, length(param_est)),
           grade=rep(desc_out$grade, length(param_est)),
           Parameter=rep(param_est, each=J))

  par_est_wide <- par_est %>%
    pivot_wider(names_from = c("Parameter"), values_from = c("Mean", "SD", "Lower95", "Upper95")) %>%
    select(studentid, occasion, grade,
           wrc.est=Mean_exp_cnt,
           secs.est=Mean_exp_tim,
           wcpm.est=Mean_wcpm,
           se.wcpm.est=SD_wcpm,
           low.95.est.wcpm=Lower95_wcpm,
           up.95.est.wcpm=Upper95_wcpm)

  final_out <- desc_out %>%
    left_join(par_est_wide) %>%
    select(student.id=studentid, occasion, grade, n.pass, numwords.total, wrc.obs, secs.obs, wcpm.obs, wrc.est, secs.est, n.pass.wcpm,
           numwords.total.wcpm, wcpm.est, se.wcpm.est, low.95.est.wcpm, up.95.est.wcpm)

  colnames(final_out) <- gsub(pattern = "est", x = colnames(final_out), replacement =bayes.soft)

  flog.info("End wcpm process with mcmc setting", name = "orfrlog")
  return(invisible(final_out))
}
