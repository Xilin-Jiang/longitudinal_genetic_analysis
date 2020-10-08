require(survival)
require(dplyr)
library(mvtnorm)


###############################
# functions for simulate disease 
# interval-censored model
###############################
generate_genetics <- function(maf, n){
  u <- runif(n)
  sapply(u, function(y) ifelse(y<maf^2, 2 , ifelse(y < (2*maf - maf^2), 1,0)))
}

simulate_genetic_profiles_ph <- function(N, num_snp,num_buckets,latent_curve, SIGMA, h_0, flag_multi_curve, gamma_shape){
  # generate risk age profile
  if(flag_multi_curve){
    # the first curve should be a flat curve
    flat_curve <- rep(mean(latent_curve), num_buckets)
    #prf1 <- rmvnorm(floor(num_snp/2), mean = flat_curve, SIGMA)
    #prf2 <- rmvnorm((num_snp - floor(num_snp/2)), mean = latent_curve, SIGMA)
    # simulate with less decreasing effect
    prf1 <- rmvnorm(num_snp - 5, mean = flat_curve, SIGMA)
    prf2 <- rmvnorm(5, mean = latent_curve, SIGMA)
    genetic_profiles <- rbind(prf1, prf2)
  }else{
    genetic_profiles <-  rmvnorm(num_snp, mean = latent_curve, SIGMA)
  }
  
  # simulate genetics for the huge population
  allele_frequency <- runif(num_snp, max = 0.5) # assuming all MAF are risk
  genetics_population <- sapply(allele_frequency, function(x) generate_genetics(x, N))
  
  # simulate unobserved effect
  u <- rgamma(N, shape = gamma_shape, rate = gamma_shape)
  
  # compute the risk age profiles for all individual
  risk_profiles <- genetics_population %*% genetic_profiles
  
  # the procedure is separate into two part: 1. simulate disease onset for each interval 
  # simulate a censoring process, then take the first event
  h_population <- exp(risk_profiles) %*% diag(h_0)
  h_population <- u * h_population # fastest opetion for adding unobserved effect
  # apply function to each one of the cell in a matrix
  failure <- apply(h_population,1:2, function(x) rexp(1, x))
  # simulate censoring
  h_censoring <- .01
  censor <- rexp(N,h_censoring)
  censor <- pmin(censor,40)
  int_length <- 40/num_buckets
  # define a function to extract the faillure time
  extract_failure <- function(x,int_length){
    id <- which(x<int_length)
    if(length(id) == 0){
      return(41) # make it censored, as the censoring will happen <=40
    }else{
      return((id[1]-1)* int_length + x[id[1]])
    }
  }
  failure <- apply(failure, 1, function(x) extract_failure(x, int_length))
  state <- 1*(censor > failure)
  time <- pmin(failure, censor)
  
  # get the survival curve and fit it
  h <- rep(0,40) 
  for(i in 1:39){
    h[i] <- sum(state*time > (i-1) & state*time < i) / sum(time > (i-1)) 
  }
  
  Beta_surv <- matrix(0, num_snp, num_buckets)
  SE_Surv <- matrix(0, num_snp,num_buckets) 
  for(gp in 1:num_buckets){
    cases <- (failure > int_length * (gp-1)) * (failure <= int_length * (gp))*state
    controls <- (failure > int_length * (gp)) * (censor > int_length * (gp-1))
    for(snp in 1:num_snp){
      genetics <- genetics_population[,snp]
      sample <- data_frame(genetics,time,cases,controls) %>% 
        filter(cases == 1 | controls == 1 ) %>%
        mutate(time = pmin(time-int_length * (gp-1), int_length))
      mySurv <- coxph(Surv(time, cases) ~ genetics, data = sample)
      SE_Surv[snp, gp] <- sqrt(mySurv$var)
      Beta_surv[snp, gp] <- mySurv$coefficients
    }
  }
  return(list(Beta_surv, SE_Surv, h))
}
# Simulate a population, while 
simulate_multivariate_ph <-function(N, num_snp,num_buckets,latent_curve, SIGMA, h_0, flag_multi_curve, gamma_shape){
  # generate risk age profile
  if(flag_multi_curve){
    # the first curve should be a flat curve
    flat_curve <- rep(mean(latent_curve), num_buckets)
    #prf1 <- rmvnorm(floor(num_snp/2), mean = flat_curve, SIGMA)
    #prf2 <- rmvnorm((num_snp - floor(num_snp/2)), mean = latent_curve, SIGMA)
    # simulate with less decreasing effect
    prf1 <- rmvnorm(num_snp - 5, mean = flat_curve, SIGMA)
    prf2 <- rmvnorm(5, mean = latent_curve, SIGMA)
    genetic_profiles <- rbind(prf1, prf2)
  }else{
    genetic_profiles <-  rmvnorm(num_snp, mean = latent_curve, SIGMA)
  }
  
  # simulate genetics for the huge population
  allele_frequency <- runif(num_snp, max = .5)
  genetics_population <- sapply(allele_frequency, function(x) generate_genetics(x, N))
  # simulate unobserved effect
  u <- rgamma(N, shape = gamma_shape, rate = gamma_shape)
  # compute the risk age profiles for all individual
  risk_profiles <- genetics_population %*% genetic_profiles
  
  # the procedure is separate into two part: 1. simulate disease onset for each interval 
  # simulate a censoring process, then take the first event
  h_population <- exp(risk_profiles) %*% diag(h_0)
  h_population <- u * h_population # fastest opetion for adding unobserved effect
  # apply function to each one of the cell in a matrix
  failure <- apply(h_population,1:2, function(x) rexp(1, x))
  # simulate censoring
  h_censoring <- .01
  censor <- rexp(N,h_censoring)
  censor <- pmin(censor,40)
  int_length <- 40/num_buckets
  # define a function to extract the faillure time
  extract_failure <- function(x,int_length){
    id <- which(x<int_length)
    if(length(id) == 0){
      return(41) # make it censored, as the censoring will happen <=40
    }else{
      return((id[1]-1)* int_length + x[id[1]])
    }
  }
  failure <- apply(failure, 1, function(x) extract_failure(x, int_length))
  state <- 1*(censor > failure)
  time <- pmin(failure, censor)
  
  # get the survival curve and fit it
  h <- rep(0,40) 
  for(i in 1:39){
    h[i] <- sum(state*time > (i-1) & state*time < i) / sum(time > (i-1)) 
  }
  
  Beta_surv <- matrix(0, num_snp, num_buckets)
  SE_Surv <- matrix(0, num_snp,num_buckets) 
  for(gp in 1:num_buckets){
    cases <- (failure > int_length * (gp-1)) * (failure <= int_length * (gp))*state
    controls <- (failure > int_length * (gp)) * (censor > int_length * (gp-1))
    # get the list of all variables
    g_list <- "V1 "
    for(g in 2:num_snp){
      g_list <- paste0(g_list , " + V", g)
    }
    genetics_population <- as.data.frame(genetics_population)
    sample <- cbind(genetics_population,data_frame(time,cases,controls)) %>% 
      filter(cases == 1 | controls == 1 ) %>%
      mutate(time = pmin(time-int_length * (gp-1), int_length))
    mySurv <- coxph(as.formula(paste0("Surv(time, cases) ~" , g_list) ), data = sample)
    SE_Surv[, gp] <- sqrt(diag(mySurv$var))
    Beta_surv[, gp] <- mySurv$coefficients
  }
  return(list(Beta_surv, SE_Surv, h))
}
# function that infer risk profile with interaction
simulate_interacted_genetic_profiles_ph <- function(N, num_snp,num_buckets,beta_0, variance_g, h_0){
  # generate risk age profile
  gamma_shape <- 1/variance_g
  genetic_profiles <-  matrix(beta_0 * rgamma(num_snp*N, shape = gamma_shape, rate = gamma_shape), nrow = num_snp, ncol = N)
  
  # simulate genetics for the huge population
  allele_frequency <- runif(num_snp, max = 0.5) # assuming all MAF are risk
  genetics_population <- sapply(allele_frequency, function(x) generate_genetics(x, N))
  
  # compute the risk age profiles for all individual
  risk_profiles <- rowSums(genetics_population * t(genetic_profiles) )
  risk_profiles <- matrix(rep(risk_profiles,num_buckets), ncol = num_buckets)
  # the procedure is separate into two part: 1. simulate disease onset for each interval 
  # simulate a censoring process, then take the first event
  h_population <- exp(risk_profiles) %*% diag(h_0)
  # apply function to each one of the cell in a matrix
  failure <- apply(h_population,1:2, function(x) rexp(1, x))
  # simulate censoring
  h_censoring <- .01
  censor <- rexp(N,h_censoring)
  censor <- pmin(censor,40)
  int_length <- 40/num_buckets
  # define a function to extract the faillure time
  extract_failure <- function(x,int_length){
    id <- which(x<int_length)
    if(length(id) == 0){
      return(41) # make it censored, as the censoring will happen <=40
    }else{
      return((id[1]-1)* int_length + x[id[1]])
    }
  }
  failure <- apply(failure, 1, function(x) extract_failure(x, int_length))
  state <- 1*(censor > failure)
  time <- pmin(failure, censor)
  
  # get the survival curve and fit it
  h <- rep(0,40) 
  for(i in 1:39){
    h[i] <- sum(state*time > (i-1) & state*time < i) / sum(time > (i-1)) 
  }
  
  Beta_surv <- matrix(0, num_snp, num_buckets)
  SE_Surv <- matrix(0, num_snp,num_buckets) 
  
  # get the list of all variables
  g_list <- "V1 "
  for(g in 2:num_snp){
    g_list <- paste0(g_list , " + V", g)
  }
  genetics_population <- as.data.frame(genetics_population)
  
  for(gp in 1:num_buckets){
    cases <- (failure > int_length * (gp-1)) * (failure <= int_length * (gp))*state
    controls <- (failure > int_length * (gp)) * (censor > int_length * (gp-1))
    sample <- cbind(genetics_population,data_frame(time,cases,controls)) %>% 
      filter(cases == 1 | controls == 1 ) %>%
      mutate(time = pmin(time-int_length * (gp-1), int_length))
    mySurv <- coxph(as.formula(paste0("Surv(time, cases) ~" , g_list) ), data = sample)
    SE_Surv[, gp] <- sqrt(diag(mySurv$var))
    Beta_surv[, gp] <- mySurv$coefficients
  }
  
  # set parameter for EM
  para <- list()
  para$p <- 3
  para$X <- X_cb_spline(para$p,8)
  para$sigma0inv <- solve(diag(para$p) * 1) 
  para$snp_lst <- as.list(1:num_snp)
  sigmasnp <- 0.0004
  para$sigmaSAMPLE <- lapply(para$snp_lst, function(x) diag(SE_Surv[x,]^2) + sigmasnp) 
  # there is cases when sigmaSAMPLE is singular, need to check
  if(any(lapply(para$sigmaSAMPLE, function(x) det(x) == 0))){
    # for sigular sigmaSAMPLE, there should be some annoying numeric issue, just ignore it
    return(NA)
  }
  para$sigmainverse <- para$sigmaSAMPLE %>%
    lapply(function(x) solve(x))
  para$betaj <- lapply(para$snp_lst, function(x) matrix(Beta_surv[x,]))
  para$sigmabeta <- mapply(function(x,y) x %*% y, para$sigmainverse, para$betaj, SIMPLIFY = FALSE) 
  para$M <- length(para$betaj[[1]])
  para$S <- length(para$betaj)
  
  # do EM
  rslt_alt_hyp <- EM_K(1,100, para)
  
  para$num_snp <- num_snp
  para$num_buckets <- num_buckets
  para$beta_0 <- beta_0
  para$variance_g <- variance_g
  para$h_0 <- h_0
  return(list(para, rslt_alt_hyp))
}

# here is a function for the threshold based model
generate_liability_disease <- function(S0, Drift, Variance, censored = 40){ # N: number of individauls 
  PATH <- matrix(nrow = length(S0), ncol = censored)
  V_matrix <- matrix(.2, censored, censored) + (1-.2) * diag(censored)
  Steps <- sapply(1:length(S0), function(x) rmvnorm(n = 1, mean = rep(Drift[x], censored),sigma =  V_matrix * Variance[x])) %>% t
  PATH <- apply(Steps, 1,  cumsum) %>% t
  PATH <- PATH + S0
  event <- apply(PATH, 1, function(x) ifelse(length(which(x > 0)) > 0, which(x > 0)[1], censored+1) )
  return(list(PATH, event))
}

###############################
# functions for EM
###############################
# polynomial -> cubic spline: with increasing degree of freedom P 
X_cb_spline <- function(P,t){
  X_base <- cbind(rep(1,t), 1:t, (1:t)^2, (1:t)^3)
  if(P <= 4){
    return(as.matrix(X_base[,1:P]))
  }else{
    X <- sapply(1+(1:(P-4))*(t-1)/(P-3), function(x) pmax(((1:t)-x)^3,0))
    return(cbind(X_base, X))
  }
}
# natural cubic spline
d_helper <- function(X, k, epsi_K, epsi_k){
  (pmax((X-epsi_k)^3,0) - pmax((X-epsi_K)^3,0))/(epsi_K - epsi_k)
}
X_nat_cb <- function(P,t){
  X_base <- cbind(rep(1,t), 1:t)
  if(P <= 2){
    return(as.matrix(X_base[,1:P]))
  }else{
    epsi <- 1+(1:P)*(t-1)/(P+1)
    hp <- function(x,k) d_helper(x, k, epsi[P], epsi[k])
    X <- matrix(nrow = t, ncol = P-2)
    for(k in 1:(P-2)){
      X[,k] <- sapply(1:t, function(x) hp(x,k) -  hp(x,P-1))
    }
    return(cbind(X_base, X))
  }
}
# function to compute the log likelihoood
comp_ll <- function(pi_saver, d_beta_hat,para){
  pi_d_beta_hat <- mapply(function(x,y) x * y, pi_saver, d_beta_hat, SIMPLIFY = FALSE)
  vector_for_sum <- sapply(1:para$S, function(j) Reduce("+", lapply(pi_d_beta_hat, function(x) x[j])))
  return(sum(log(vector_for_sum)))
}
# function to update the E step, update z distribtion
comp_z <- function(pi_saver, d_beta_hat){
  pi_d_beta_hat <- mapply(function(x,y) x * y, pi_saver, d_beta_hat, SIMPLIFY = FALSE)
  d_sum <- Reduce('+', pi_d_beta_hat)
  z_prob <- lapply(pi_d_beta_hat, function(x) x/d_sum)
  return(z_prob)
}
# function to update theta
# No prior sigma0 here; frequentist 
comp_theta <- function(z_prob, para){
  z_betaj_sum <- lapply(z_prob, function(z) Reduce("+", lapply(1:para$S, function(j) z[j] * t(para$X) %*% para$sigmabeta[[j]])))
  z_sigma_sum <- lapply(z_prob, function(z) para$sigma0inv + Reduce("+", lapply(1:para$S, function(j) z[j] * t(para$X) %*% para$sigmainverse[[j]] %*% para$X)))
  theta <- mapply(function(x,y) solve(x, y), z_sigma_sum, z_betaj_sum, SIMPLIFY = FALSE)
  return(theta)
}
# function to update pi_saver
comp_pi <- function(z_prob){
  Ni <- lapply(z_prob, sum)
  N <- Reduce("+", Ni)
  pi_saver <- lapply(Ni, function(x) x/N)
  return(pi_saver)
}
EM_K <- function(K, num_itr, para){ # flat flag here refer to whether the first line need to be flat
  # para$X will determine the fitting: if X is constant, it is just fitting a constant effect
  # para$X will determine the fitting: if X is frailty, it is just fitting a frailty effect
  # initialize parameter
  pi_saver <- list()
  theta <- list()
  BETA <- list()
  d_beta_hat <- list()
  for(j in 1:K){
    x <- runif(K)
    pi_saver[[j]] <- x[j]/sum(x)
    # pi_saver[[j]] <- 1/K
    theta[[j]] <-  matrix(rnorm(para$p,sd=.0001),para$p,1) # using a small sd value to avoid NA
    BETA[[j]] <- para$X %*% theta[[j]]
    d_beta_hat[[j]] <- sapply(1:para$S, function(x) 
      dmvnorm(t(para$betaj[[x]]), mean = BETA[[j]], sigma = para$sigmaSAMPLE[[x]]))
  }
  saver_ll <- matrix(0, num_itr, 1)
  for(itr in 1:num_itr){
    if(itr %% 50 ==0) print(itr)
    # compute ll
    ll <- comp_ll(pi_saver, d_beta_hat, para)
    if(is.na(ll)) break
    # assert_that(!is.na(ll),  msg = "NA values found!")
    saver_ll[itr,1] <- ll
    # M step
    z_prob <- comp_z(pi_saver, d_beta_hat)
    # E step
    theta <- comp_theta(z_prob, para)
    pi_saver <- comp_pi(z_prob)
    for(j in 1:K){
      BETA[[j]] <- para$X %*% theta[[j]]
      d_beta_hat[[j]] <- sapply(1:para$S, function(x) 
        dmvnorm(t(para$betaj[[x]]), mean = BETA[[j]], sigma = para$sigmaSAMPLE[[x]]))
    }
  }
  return(list(pi_saver, z_prob, theta, BETA, saver_ll, ll))
}

# this function is used to load a SNP data from csv file, and return a para which could be used directly in EM fitting 
load_SNP_coef <- function(SNP_FILE, FLIP_FILE){
  snp_for_clustering <- read.csv(SNP_FILE)
  flip_snp <- read.csv(FLIP_FILE)
  
  useful_snp_data <- snp_for_clustering %>% 
    # I need to remove the NAs
    mutate(SE = ifelse(is.na(SE), 10, SE)) %>%
    mutate(OR = ifelse(is.na(OR), 1, OR)) %>%
    mutate(OR = log(OR)) %>%
    # I will have to control the overflowing of coef (some -> inf)
    mutate(SE = ifelse(abs(OR) > 10, 10, SE)) %>%
    mutate(OR = ifelse(abs(OR) > 10, 0, OR)) %>% 
    mutate(OR = ifelse(SNP %in% flip_snp$SNP, -OR, OR))
  
  para <- list()
  para$snp_lst <- useful_snp_data %>% 
    group_by(SNP)  %>%
    summarise(n()) %>%
    select(SNP) 
  sl <- para$snp_lst[["SNP"]]
  para$snp_lst <- as.list(sl)
  names(para$snp_lst) <- sl
  
  para$p <- 3
  para$X <- X_cb_spline(para$p,8) # default X
  para$sigma0inv <- solve(diag(para$p) * 1) 
  
  sigmasnp <- 0.0004
  para$sigmaSAMPLE <- para$snp_lst %>%
    lapply(function(x) dplyr::filter(useful_snp_data, SNP == x)) %>%
    lapply(function(x) diag(x[["SE"]]^2) + sigmasnp)
  
  para$sigmainverse <- para$sigmaSAMPLE %>%
    lapply(function(x) solve(x))
  
  para$betaj <- para$snp_lst %>%
    lapply(function(x) dplyr::filter(useful_snp_data, SNP == x)) %>%
    lapply(function(x) matrix(x[["OR"]]))
  
  para$sigmabeta <- mapply(function(x,y) x %*% y, para$sigmainverse, para$betaj, SIMPLIFY = FALSE) 
  para$M <- length(para$betaj[[1]])
  para$S <- length(para$betaj)

  return(para)
}


