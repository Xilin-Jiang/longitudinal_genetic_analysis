library(mvtnorm)
library(dplyr)
library(survival)
library(ggplot2)

cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
# set the commonly used color
blue <- cbPalette[6]
red <- cbPalette[7]
green <- cbPalette[4]
orange <- cbPalette[2]
grey <- cbPalette[1]
yellow <- cbPalette[5]
purple <- cbPalette[8]
skyblue <- cbPalette[3]


args <- commandArgs(trailingOnly = TRUE) 
s_id <- as.numeric(args[1])

# checking the source of the bias
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
  failure <- apply(failure, 1, function(x)extract_failure(x, int_length))
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


# then do the clustering
X_cb_spline <- function(P,t){
  X_base <- cbind(rep(1,t), 1:t, (1:t)^2, (1:t)^3)
  if(P <= 4){
    return(as.matrix(X_base[,1:P]))
  }else{
    X <- sapply(1+(1:(P-4))*(t-1)/(P-3), function(x) pmax(((1:t)-x)^3,0))
    return(cbind(X_base, X))
  }
}

comp_ll <- function(pi_saver, d_beta_hat,para){
  pi_d_beta_hat <- mapply(function(x,y) x * y, pi_saver, d_beta_hat, SIMPLIFY = FALSE)
  vector_for_sum <- sapply(1:para$S, function(j) Reduce("+", lapply(pi_d_beta_hat, function(x) x[j])))
  return(sum(log(vector_for_sum)))
}
comp_z <- function(pi_saver, d_beta_hat){
  pi_d_beta_hat <- mapply(function(x,y) x * y, pi_saver, d_beta_hat, SIMPLIFY = FALSE)
  d_sum <- Reduce('+', pi_d_beta_hat)
  z_prob <- lapply(pi_d_beta_hat, function(x) x/d_sum)
  return(z_prob)
}
comp_theta <- function(z_prob, para){
  sigma0inv <- solve(matrix(.2, para$p, para$p) + (1-.2) * diag(para$p)) # add prior, super unefficient but good for now
  z_betaj_sum <- lapply(z_prob, function(z) rowSums(sapply(1:para$S, function(j) z[j] * t(para$X) %*% para$sigmabeta[[j]])))
  z_sigma_sum <- lapply(z_prob, function(z) sigma0inv + Reduce("+", lapply(1:para$S, function(j) z[j] * t(para$X) %*% para$sigmainverse[[j]] %*% para$X)))
  theta <- mapply(function(x,y) solve(x, y), z_sigma_sum, z_betaj_sum, SIMPLIFY = FALSE)
  return(theta)
}
comp_pi <- function(z_prob){
  Ni <- lapply(z_prob, sum)
  N <- Reduce("+", Ni)
  pi_saver <- lapply(Ni, function(x) x/N)
  return(pi_saver)
}
EM_K <- function(K, num_itr, flat_flag, para){ # flat flag here refer to whether the first line need to be flat
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
    if(j == 1 && flat_flag){
      # always force the first component to be a flat line
      BETA[[j]] <- matrix(0, para$M, 1)
    }
    else BETA[[j]] <- para$X %*% theta[[j]]
    d_beta_hat[[j]] <- sapply(1:para$S, function(x) 
      dmvnorm(t(para$betaj[[x]]), mean = BETA[[j]], sigma = para$sigmaSAMPLE[[x]]))
  }
  saver_ll <- matrix(0, num_itr, 1)
  for(itr in 1:num_itr){
    if(itr %% 10 ==0) print(itr)
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
      if(j == 1 && flat_flag){
        # always force the first component to be a flat line
        z <- z_prob[[1]]
        beta0 <- sum(sapply(1:para$S, function(j) z[j] * para$sigmabeta[[j]]))/sum(sapply(1:para$S, function(j) z[j] * para$sigmainverse[[j]]))
        BETA[[j]] <- matrix(beta0, para$M, 1)
      }
      else BETA[[j]] <- para$X %*% theta[[j]]
      d_beta_hat[[j]] <- sapply(1:para$S, function(x) 
        dmvnorm(t(para$betaj[[x]]), mean = BETA[[j]], sigma = para$sigmaSAMPLE[[x]]))
    }
  }
  return(list(pi_saver, z_prob, theta, BETA, saver_ll, ll))
}
# the first simulation function is to simulate single curve basic case
simulation <- function(N, num_snp,num_buckets,latent_curve, SIGMA, h_0, mv_flag, mtcv_flag, gamma_shape){
  # mv_flag determine whether the estimate is univariate or multivariate
  # mtcv_flag determine what hypothesis this is testing
  # simulate with multivariate genetics effect
  if(mv_flag){
    rslt_ph <- simulate_multivariate_ph(N, num_snp,num_buckets,latent_curve, SIGMA, h_0, flag_multi_curve = mtcv_flag, gamma_shape)
  }else{
    rslt_ph <- simulate_genetic_profiles_ph(N, num_snp,num_buckets,latent_curve, SIGMA, h_0, flag_multi_curve = mtcv_flag, gamma_shape)
  }
  
  Beta_surv <- rslt_ph[[1]]
  SE_Surv <- rslt_ph[[2]]
  h <- rslt_ph[[3]]
  # use para to save all the inference parameters
  para <- list()
  para$p <- 3
  para$X <- X_cb_spline(para$p,8)
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
  if(mtcv_flag){
    rslt_alt_hyp <- EM_K(2,100, FALSE, para)
    rslt_null_hyp <- EM_K(1,100, FALSE, para)
  }else{
    rslt_alt_hyp <- EM_K(1,100, FALSE, para)
    rslt_null_hyp <- EM_K(1,100, TRUE, para)
  }
  # save some important variable
  para$num_snp <- num_snp
  para$num_buckets <- num_buckets
  para$latent_curve <- latent_curve
  para$SIGMA <- SIGMA
  para$h_0 <- h_0
  return(list(para, rslt_ph, rslt_alt_hyp, rslt_null_hyp))
}

set.seed(19940110 + s_id)
num_buckets <- 8
N <- 50000 # for the surviaval analysis, after case matching, we usually only have 50,000 individual left

h_censoring <- .01
h_0 <- (5 + (1:8))/20000 # underlying risk should be increasing over time
# h_0 <- (1:8)*0.0001
beta_0 <- 0.1

num_snp <- 50

# first use different variance
sigmasnp <- 0.0004
SIGMA <- 10^(-4) * diag(num_buckets) + sigmasnp
latent_curve_1 <- rep(beta_0, num_buckets)
######################################################
# simulation for the test of non-flat effect
######################################################
SIM <- list()
s_lst <- (-100:100)/10000
s <- s_lst[s_id]
latent_curve <- latent_curve_1 + s * (-3.5:3.5)
# do many simulation for a single latent curve
for(j in 1:100){
  try(SIM[[j]] <- simulation(N,num_snp,num_buckets,latent_curve, SIGMA, h_0, mv_flag = TRUE, mtcv_flag = FALSE, 10^15))
}

save(SIM, file = paste0("Non_flat_test_jointly/joint_non_flat_test",s_id,".Rdata"))

# get the power test figure for log-likelihood testing
p_non_flat <- function(SIM){
  rslt_alt_hyp <- SIM[[3]]
  rslt_null_hyp <- SIM[[4]]
  ratio_ll <- 2*(rslt_alt_hyp[[6]] - rslt_null_hyp[[6]])
  return(1 - pchisq(ratio_ll, df=2))
}
p_non_flat_test <- matrix(nrow = length(s_lst), ncol = 100)
for(i in 1:201){
  load(paste0("simulation/joint_non_flat_test",i,".Rdata"))
  p_non_flat_test[i,1:length(SIM)] <- sapply(1:length(SIM), function(x) ifelse(length(SIM[[x]]) == 4,p_non_flat(SIM[[x]]), NA))
}
p_thre <- 0.05
power_mean <- rowMeans(p_non_flat_test < p_thre,na.rm = T)
power_var <- sqrt(power_mean*(1-power_mean)/100)
df_non_flat_test <- data.frame("Slope" = s_lst, "Power" = power_mean, "Power_std" = power_var)
plt <- ggplot(df_non_flat_test) + 
  geom_line(aes(x=Slope, y = Power), color = blue) +
  geom_ribbon(aes(x=Slope, ymax = Power + 1.96*Power_std, ymin = Power - 1.96*Power_std ), alpha = 0.3) +
  # geom_smooth(aes(x=Slope, y = P), span = 0.1) + 
  theme(legend.position = "none",panel.background=element_blank()) +
  ylab("Power")

ggsave("~/Desktop/Writting_up_genetic_risk/Joint_estimation_non_flat.png", width = 6, height = 4)

# do the plotting of latent curve
SIM1 <- SIM[[2]]
para <- SIM1[[1]]
rslt_alt_hyp <- SIM1[[3]]
rslt_null_hyp <- SIM1[[4]]
BETA <- rslt_alt_hyp[[4]]
z_prob <- rslt_alt_hyp[[2]]
P <- .2
sigma0 <- matrix(P, para$p, para$p) + (1-P) * diag(para$p)
sigma0inv <- solve(sigma0)
z_asign <- sapply(1:para$S, function(s) which.max(lapply(z_prob, function(x) x[s])))
sigma_i <- list()
for(j in 1:length(z_prob)){
  ZSigma <- mapply(function(x,y) x*y, z_prob[[j]],  para$sigmainverse,SIMPLIFY = FALSE)
  sum_sigmainv <- Reduce('+', ZSigma)
  Aj <- sigma0inv + t(para$X) %*% sum_sigmainv %*% para$X
  sigma_i[[j]] <- para$X %*% solve(Aj) %*% t(para$X)
}
df_mean <- data.frame("PHE"= 1:8)
for(i in 1:length(BETA)){
  df_mean[paste("beta_", i, sep = "")] = BETA[[i]]
  df_mean[paste("variance_", i, sep = "")] = sqrt(diag(sigma_i[[i]]))
}
df_mean["true_curve"] <- para$latent_curve
plt <- ggplot(df_mean, aes(PHE)) + xlab("Age phase (per 5 year from <45 to >75)") +
  theme(legend.position = "none",panel.background=element_blank()) + 
  xlab("Age group") + ylab("Log Odds Ratio") + 
  geom_line(aes(y = beta_1), colour = red) + 
  geom_ribbon(aes(ymin = beta_1 - 2*variance_1, ymax = beta_1 + 2*variance_1), fill = red, alpha = 0.3) +
  geom_line(aes(y = true_curve), colour = blue, linetype = "dashed", size = 3) + ylim(c(0.04, 0.15))
ggsave("~/Desktop/Writting_up_genetic_risk/simulation_fitting_multivariate_1.png", width = 4, height = 4)
#########################################
# 2020-2-17 get the coverage of interval
#########################################
P <- .2
sigma0 <- matrix(P, 3, 3) + (1-P) * diag(3)
sigma0inv <- solve(sigma0)
comp_covrg <- function(SIM, cv){
  if(!length(SIM) == 4) return(matrix(nrow = length(cv), ncol = 1))
  BETA <- SIM[[3]][[4]][[1]]
  para <- SIM[[1]]
  sum_sigmainv <- Reduce('+', para$sigmainverse)
  Aj <- sigma0inv + t(para$X) %*% sum_sigmainv %*% para$X
  invAj <- para$X %*% solve(Aj) %*% t(para$X)
  upbeta <- BETA + 2*sqrt(diag(invAj))
  lwbeta <- BETA - 2*sqrt(diag(invAj))
  return( (cv < upbeta) & (cv > lwbeta) )
}
p_covrg <- matrix(nrow = length(s_lst), ncol = 8)
for(i in 1:201){
  load(paste0("simulation/joint_non_flat_test",i,".Rdata"))
  cv <- ifelse(is.na(SIM[[1]][[1]]), SIM[[2]][[1]]$latent_curve,SIM[[1]][[1]]$latent_curve)
  p_covrg[i,] <- rowMeans(sapply(1:length(SIM), function(x) comp_covrg(SIM[[x]], cv) ), na.rm = TRUE)
}
power_mean <- rowMeans(p_covrg)
power_std <- sqrt(power_mean*(1-power_mean)/800)
df_covrg <- data.frame("Slope" = s_lst, "Power" = power_mean, "Power_std" = power_std)
plt <- ggplot(df_covrg) + 
  geom_line(aes(x=Slope, y = Power), color = blue) +
  geom_ribbon(aes(x=Slope, ymax = Power + 1.96*Power_std, ymin = Power - 1.96*Power_std ), alpha = 0.3) +
  # geom_smooth(aes(x=Slope, y = P), span = 0.1) + 
  theme(legend.position = "none",panel.background=element_blank()) +
  ylab("Coverage")

ggsave("~/Desktop/Writting_up_genetic_risk/coverage_over_slope.png" ,plt, width = 6, height = 4)
######################################################
# do the multiple cluster testing
######################################################
SIM <- list()
s_lst <- (-150:150)/4000
s <- s_lst[s_id]
latent_curve <- latent_curve_1 + s * (-3.5:3.5)
# do many simulation for a single latent curve
for(j in 1:100){
  try(SIM[[j]] <- simulation(N,num_snp,num_buckets,latent_curve, SIGMA, h_0, mv_flag = TRUE, mtcv_flag = TRUE, 10^15))
}

save(SIM, file =  paste0("multi_curve_test_jointly/multi_curve_test_jointly_",s_id,".Rdata"))

# get the power test figure for log-likelihood testing
p_multi_cluster <- function(SIM){
  rslt_alt_hyp <- SIM[[3]]
  rslt_null_hyp <- SIM[[4]]
  ratio_ll <- 2*(rslt_alt_hyp[[6]] - rslt_null_hyp[[6]])
  # degree of freedom is 5
  return(1 - pchisq(ratio_ll, df=4)) 
}
# due to the numeric error
s_lst <- (-150:150)/4000 # could be removed
p_multi_curve_test <- matrix(nrow = length(s_lst), ncol = 100)
for(i in 1:301){ # for(i in 1:201)
  load(paste0("simulation/multi_curve_test_jointly_",i,".Rdata"))
  # p_multi_curve_test[i,] <- sapply(1:length(SIM), function(x) ifelse(anyNA(SIM[[x]]),NA,p_multi_cluster(SIM[[x]])))
  p_multi_curve_test[i,] <- sapply(1:length(SIM), function(x) ifelse(anyNA(SIM[[x]]),NA,p_multi_cluster(SIM[[x]])))
}
p_thre <- 0.05
power_mean <- rowMeans(p_multi_curve_test< p_thre, na.rm = TRUE)
power_var <- sqrt(power_mean*(1-power_mean)/99)
df_multi_curv_test <- data.frame("Slope" = s_lst, "Power" = power_mean, "Power_std" = power_var)
plt <- ggplot(df_multi_curv_test) + 
  geom_line(aes(x=Slope, y = Power), color = red) +
  geom_ribbon(aes(x=Slope, ymax = Power + 1.96*Power_std, ymin = Power - 1.96*Power_std ), alpha = 0.3) +
  # geom_smooth(aes(x=Slope, y = P), span = 0.1) + 
  theme(legend.position = "none",panel.background=element_blank()) +
  ylab("Power")
ggsave("~/Desktop/Writting_up_genetic_risk/Joint_estimation_multiple_curve.png", width = 6, height = 4)
# compute the classification accuracy
compute_accuracy <- function(SIM){
  if(length(SIM) != 4) return(NA)
  rslt_alt_hyp <- SIM[[3]]
  para <- SIM[[1]]
  BETA <- rslt_alt_hyp[[4]]
  num_non_flat <- 5
  latent_curve <- para$latent_curve
  z_prob <- rslt_alt_hyp[[2]]
  z_asign <- sapply(1:para$S, function(s) which.max(lapply(z_prob, function(x) x[s])))
  if(sum((BETA[[1]] - latent_curve)^2) > sum((BETA[[2]] - latent_curve)^2)){
    true_assign <- c(rep(1,para$num_snp-num_non_flat),rep(2,num_non_flat))
    accu <- sum(z_asign == true_assign)/para$num_snp
  }else{
    true_assign <- c(rep(2,para$num_snp-num_non_flat),rep(1,num_non_flat))
    accu <- sum(z_asign == true_assign)/para$num_snp
  }
  return(accu)
}
sensitivity <- matrix(nrow = length(s_lst), ncol = 100)
for(i in 1:201){
  load(paste0("simulation/multi_curve_test_jointly",i,".Rdata"))
  sensitivity[i,] <- sapply(1:length(SIM), function(x) ifelse(anyNA(SIM[[x]]),NA,compute_accuracy(SIM[[x]])))
}
sens_mean <- rowMeans(sensitivity, na.rm = TRUE)
power_var <- sapply(1:length(s_lst), function(x) sqrt(var(sensitivity[x,], na.rm = TRUE)/100))
df_non_flat_test <- data.frame("Slope" = s_lst, "Power" = sens_mean, "Power_std" = power_var)
ggplot(df_non_flat_test) + 
  geom_line(aes(x=Slope, y = Power), color = red) +
  geom_ribbon(aes(x=Slope, ymax = Power + 1.96*Power_std, ymin = Power - 1.96*Power_std ), alpha = 0.3) +
  # geom_smooth(aes(x=Slope, y = P), span = 0.1) + 
  theme(legend.position = "none",panel.background=element_blank()) +
  ylab("Power")


# do some plotting
SIM2 <- SIM[[1]]
para <- SIM2[[1]]
rslt_alt_hyp <- SIM2[[3]]
rslt_null_hyp <- SIM2[[4]]
BETA <- rslt_alt_hyp[[4]]
z_prob <- rslt_alt_hyp[[2]]
P <- .2
sigma0 <- matrix(P, para$p, para$p) + (1-P) * diag(para$p)
sigma0inv <- solve(sigma0)
z_asign <- sapply(1:para$S, function(s) which.max(lapply(z_prob, function(x) x[s])))
sigma_i <- list()
for(j in 1:length(z_prob)){
  ZSigma <- mapply(function(x,y) x*y, z_prob[[j]],  para$sigmainverse,SIMPLIFY = FALSE)
  sum_sigmainv <- Reduce('+', ZSigma)
  Aj <- sigma0inv + t(para$X) %*% sum_sigmainv %*% para$X
  sigma_i[[j]] <- para$X %*% solve(Aj) %*% t(para$X)
}
df_mean <- data.frame("PHE"= 1:8)
for(i in 1:length(BETA)){
  df_mean[paste("beta_", i, sep = "")] = BETA[[i]]
  df_mean[paste("variance_", i, sep = "")] = sqrt(diag(sigma_i[[i]]))
}
df_mean[["true_curve"]] <- para$latent_curve
df_mean[["true_flat"]] <- rep(mean(para$latent_curve), num_buckets)
plt <- ggplot(df_mean, aes(PHE)) + xlab("Age phase (per 5 year from <45 to >75)") +
  theme(legend.position = "none",panel.background=element_blank()) + 
  xlab("Age group") + ylab("Log Odds Ratio") + 
  geom_line(aes(y = beta_1), colour = red) + 
  geom_ribbon(aes(ymin = beta_1 - 2*variance_1, ymax = beta_1 + 2*variance_1), fill = red, alpha = 0.3) +
  geom_line(aes(y = beta_2), colour = blue) + 
  geom_ribbon(aes(ymin = beta_2 - 2*variance_2, ymax = beta_2 + 2*variance_2), fill = blue, alpha = 0.3) +
  geom_line(aes(y = true_curve), colour = blue, linetype = "dashed", size = 3) + 
  geom_line(aes(y = true_flat), colour = red, linetype = "dashed", size = 3) 
ggsave("~/Desktop/Writting_up_genetic_risk/Simulation_two_curve_joint_1.png", width = 4, height = 4)


h_poly <- function(t, b,gamma, k){
  (gamma * t^k )/(1+ gamma/((k+1)*b)*t^(k+1) )
}

# do the multivariate simulation 
sim <- list()
for(i in 1:10){
  print(i)
  sim[[i]] <- simulation(num_snp,num_buckets,latent_curve_2, SIGMA, h_0, mv_flag = TRUE, 10^15)
}

# this is for simple visualization
sim <- list()
for(i in 1:10){
  print(i)
  sim[[i]] <- simulate_genetic_profiles_ph(N, num_snp,num_buckets,latent_curve_2, SIGMA, h_0, flag_multi_curve = FALSE, 10^15)
}
sumsim <- matrix(0, num_snp, num_buckets)
for(i in c(1:10)){
  sumsim <- sumsim + sim[[i]][[1]]
}
plot(colMeans(sumsim)/10)

#####################################
# simulation with different frailty
#####################################
SIM <- list()
u_lst <- exp(0.1*(1:100))
latent_curve <- latent_curve_1 
gamma_shape <- u_lst[s_id]
# do many simulation for a single latent curve
for(j in 1:100){
  try(SIM[[j]] <- simulation(N,num_snp,num_buckets,latent_curve, SIGMA, h_0, mv_flag = TRUE, mtcv_flag = FALSE, gamma_shape))
}

save(SIM, file = paste0("frailty_test/frailty_test",s_id,".Rdata"))

# get the power test figure for log-likelihood testing
p_non_flat <- function(SIM){
  rslt_alt_hyp <- SIM[[3]]
  rslt_null_hyp <- SIM[[4]]
  ratio_ll <- 2*(rslt_alt_hyp[[6]] - rslt_null_hyp[[6]])
  return(1 - pchisq(ratio_ll, df=2))
}
p_frailty <- matrix(nrow = 100, ncol = 100)
for(i in c(1:100)){
  load(paste0("simulation/jointly_frailty_test",i,".Rdata"))
  p_frailty[i,] <- sapply(1:length(SIM), function(x) ifelse(anyNA(SIM[[x]]),NA,p_non_flat(SIM[[x]])))
}
p_thre <- 0.05
power_mean <- rowMeans(p_frailty < p_thre, na.rm = TRUE)
power_var <- sqrt(power_mean*(1-power_mean)/100)
df_frailty <- data.frame("frailty" = 1/u_lst, "Power" = power_mean, "Power_std" = power_var)
plt <- ggplot(df_frailty) + 
  geom_line(aes(x=frailty, y = Power), color = blue) +
  geom_ribbon(aes(x=frailty, ymax = Power + 1.96*Power_std, ymin = Power - 1.96*Power_std ), alpha = 0.3) +
  # geom_smooth(aes(x=Slope, y = P), span = 0.1) + 
  theme(legend.position = "none",panel.background=element_blank()) +
  ylab("Power") + 
  scale_x_continuous(trans='log10') +
  xlab("Gamma scale parameter (variance of risk)")
ggsave("~/Desktop/Writting_up_genetic_risk/frailty_power.png" ,plt, width = 6, height = 4)
###############################
# plot the effect of frailty
###############################
load("simulation/jointly_frailty_test2.Rdata")
SIM1 <- SIM[[1]]
para <- SIM1[[1]]
rslt_alt_hyp <- SIM1[[3]]
rslt_null_hyp <- SIM1[[4]]
BETA <- rslt_alt_hyp[[4]]
z_prob <- rslt_alt_hyp[[2]]
P <- .2
sigma0 <- matrix(P, para$p, para$p) + (1-P) * diag(para$p)
sigma0inv <- solve(sigma0)
z_asign <- sapply(1:para$S, function(s) which.max(lapply(z_prob, function(x) x[s])))
sigma_i <- list()
for(j in 1:length(z_prob)){
  ZSigma <- mapply(function(x,y) x*y, z_prob[[j]],  para$sigmainverse,SIMPLIFY = FALSE)
  sum_sigmainv <- Reduce('+', ZSigma)
  Aj <- sigma0inv + t(para$X) %*% sum_sigmainv %*% para$X
  sigma_i[[j]] <- para$X %*% solve(Aj) %*% t(para$X)
}
df_mean <- data.frame("PHE"= 1:8)
for(i in 1:length(BETA)){
  df_mean[paste("beta_", i, sep = "")] = BETA[[i]]
  df_mean[paste("variance_", i, sep = "")] = sqrt(diag(sigma_i[[i]]))
}
df_mean["true_curve"] <- para$latent_curve
plt <- ggplot(df_mean, aes(PHE)) + xlab("Age phase (per 5 year from <45 to >75)") +
  theme(legend.position = "none",panel.background=element_blank()) + 
  xlab("Age group") + ylab("Log Odds Ratio") + 
  geom_line(aes(y = beta_1), colour = red) + 
  geom_ribbon(aes(ymin = beta_1 - 2*variance_1, ymax = beta_1 + 2*variance_1), fill = red, alpha = 0.3) +
  geom_line(aes(y = true_curve), colour = blue, linetype = "dashed", size = 3)
ggsave("~/Desktop/Writting_up_genetic_risk/Frailty_effect_shape1.png" ,plt, width = 4, height = 4)

# do a simulation with no generating variance, making all effect exactly as the latent profiles
load("simulation/multi_curve_test_jointly1.Rdata")
h1 <- SIM[[1]][[2]][[3]]
load("simulation/multi_curve_test_jointly200.Rdata")
h2 <- SIM[[1]][[2]][[3]]
df_hazard_over_age <- data.frame("AGE" = 1:40, "hazard1" = h1,"hazard2" = h2) 
ggplot(df_hazard_over_age,aes(AGE)) + 
  geom_line(aes(y = hazard1, colour = "hazard of decreasing risk"), size = 2) + 
  geom_line(aes(y = hazard2, colour = "hazard of increasing risk"), size = 2) 


######################################################
# test the effect of a healthier population
######################################################
SIM <- list()
s_lst <- (-100:100)/100
s <- s_lst[s_id]
latent_curve <- latent_curve_1 
h_0 <- (9.5 + s * (-3.5:3.5))/10000 # underlying risk should be increasing over time
# do many simulation for a single latent curve
for(j in 1:100){
  try(SIM[[j]] <- simulation(N,num_snp,num_buckets,latent_curve, SIGMA, h_0, mv_flag = TRUE, mtcv_flag = FALSE, 10^15))
}

save(SIM, file = paste0("healthy_elder_effect/joint_healthy_elder_test",s_id,".Rdata"))

# get the power test figure for log-likelihood testing
p_non_flat <- function(SIM){
  rslt_alt_hyp <- SIM[[3]]
  rslt_null_hyp <- SIM[[4]]
  ratio_ll <- 2*(rslt_alt_hyp[[6]] - rslt_null_hyp[[6]])
  return(1 - pchisq(ratio_ll, df=2))
}
p_non_flat_test <- matrix(nrow = length(s_lst), ncol = 100)
for(i in 1:201){
  load(paste0("simulation/joint_healthy_elder_test",i,".Rdata"))
  p_non_flat_test[i,] <- sapply(1:length(SIM), function(x) p_non_flat(SIM[[x]]))
}
p_thre <- 0.05
power_mean <- rowMeans(p_non_flat_test < p_thre)
power_var <- sqrt(power_mean*(1-power_mean)/100)
df_age_effect <- data.frame("Slope" = s_lst, "Power" = power_mean, "Power_std" = power_var)
plt <- ggplot(df_age_effect) + 
  geom_line(aes(x=Slope, y = Power), color = blue) +
  geom_ribbon(aes(x=Slope, ymax = Power + 1.96*Power_std, ymin = Power - 1.96*Power_std ), alpha = 0.3) +
  # geom_smooth(aes(x=Slope, y = P), span = 0.1) + 
  theme(legend.position = "none",panel.background=element_blank()) +
  ylab("Power")

ggsave("~/Desktop/Writting_up_genetic_risk/Age_effect.png" ,plt, width = 6, height = 4)

SIM1 <- SIM[[3]]
para <- SIM1[[1]]
rslt_alt_hyp <- SIM1[[3]]
rslt_null_hyp <- SIM1[[4]]
BETA <- rslt_alt_hyp[[4]]
z_prob <- rslt_alt_hyp[[2]]
P <- .2
sigma0 <- matrix(P, para$p, para$p) + (1-P) * diag(para$p)
sigma0inv <- solve(sigma0)
z_asign <- sapply(1:para$S, function(s) which.max(lapply(z_prob, function(x) x[s])))
sigma_i <- list()
for(j in 1:length(z_prob)){
  ZSigma <- mapply(function(x,y) x*y, z_prob[[j]],  para$sigmainverse,SIMPLIFY = FALSE)
  sum_sigmainv <- Reduce('+', ZSigma)
  Aj <- sigma0inv + t(para$X) %*% sum_sigmainv %*% para$X
  sigma_i[[j]] <- para$X %*% solve(Aj) %*% t(para$X)
}
df_mean <- data.frame("PHE"= 1:8)
for(i in 1:length(BETA)){
  df_mean[paste("beta_", i, sep = "")] = BETA[[i]]
  df_mean[paste("variance_", i, sep = "")] = sqrt(diag(sigma_i[[i]]))
}
df_mean["true_curve"] <- para$latent_curve
plt <- ggplot(df_mean, aes(PHE)) + xlab("Age phase (per 5 year from <45 to >75)") +
  theme(legend.position = "none",panel.background=element_blank()) + 
  xlab("Age group") + ylab("Log Odds Ratio") + 
  geom_line(aes(y = beta_1), colour = red) + 
  geom_ribbon(aes(ymin = beta_1 - 2*variance_1, ymax = beta_1 + 2*variance_1), fill = red, alpha = 0.3) +
  geom_line(aes(y = true_curve), colour = blue, linetype = "dashed", size = 3) + ylim(c(0.05, 0.15))
plt

#####################################################################################
# test the for any difference of profile in univariate v.s. multivariate estimation
#####################################################################################
s <- -0.0
latent_curve <- latent_curve_1 + s * (-3.5:3.5)
SIM_univariate <- simulation(N,num_snp,num_buckets,latent_curve, SIGMA, h_0, mv_flag = FALSE, mtcv_flag = FALSE, 10^15)
SIM_multivariate <- simulation(N,num_snp,num_buckets,latent_curve, SIGMA, h_0, mv_flag = TRUE, mtcv_flag = FALSE, 10^15)

h_uni <- SIM_univariate[[2]][[3]]
prod(1-h_uni)

h_multi <- SIM_multivariate[[2]][[3]]
prod(1-h_multi)

para<- SIM_univariate[[1]]
rslt_alt_hyp <- SIM_univariate[[3]]
rslt_null_hyp <- SIM_univariate[[4]]
BETA <- rslt_alt_hyp[[4]]
z_prob <- rslt_alt_hyp[[2]]
P <- .2
sigma0 <- matrix(P, para$p, para$p) + (1-P) * diag(para$p)
sigma0inv <- solve(sigma0)
z_asign <- sapply(1:para$S, function(s) which.max(lapply(z_prob, function(x) x[s])))
sigma_i <- list()
for(j in 1:length(z_prob)){
  ZSigma <- mapply(function(x,y) x*y, z_prob[[j]],  para$sigmainverse,SIMPLIFY = FALSE)
  sum_sigmainv <- Reduce('+', ZSigma)
  Aj <- sigma0inv + t(para$X) %*% sum_sigmainv %*% para$X
  sigma_i[[j]] <- para$X %*% solve(Aj) %*% t(para$X)
}
df_mean <- data.frame("PHE"= 1:8)
for(i in 1:length(BETA)){
  df_mean[paste("beta_univariate_", i, sep = "")] = BETA[[i]]
  df_mean[paste("variance_univariate_", i, sep = "")] = sqrt(diag(sigma_i[[i]]))
}
df_mean["true_curve"] <- para$latent_curve

para<- SIM_multivariate[[1]]
rslt_alt_hyp <- SIM_multivariate[[3]]
rslt_null_hyp <- SIM_multivariate[[4]]
BETA <- rslt_alt_hyp[[4]]
z_prob <- rslt_alt_hyp[[2]]
P <- .2
sigma0 <- matrix(P, para$p, para$p) + (1-P) * diag(para$p)
sigma0inv <- solve(sigma0)
z_asign <- sapply(1:para$S, function(s) which.max(lapply(z_prob, function(x) x[s])))
sigma_i <- list()
for(j in 1:length(z_prob)){
  ZSigma <- mapply(function(x,y) x*y, z_prob[[j]],  para$sigmainverse,SIMPLIFY = FALSE)
  sum_sigmainv <- Reduce('+', ZSigma)
  Aj <- sigma0inv + t(para$X) %*% sum_sigmainv %*% para$X
  sigma_i[[j]] <- para$X %*% solve(Aj) %*% t(para$X)
}
for(i in 1:length(BETA)){
  df_mean[paste("beta_multivariate_", i, sep = "")] = BETA[[i]]
  df_mean[paste("variance_multivariate_", i, sep = "")] = sqrt(diag(sigma_i[[i]]))
}

plt <- ggplot(df_mean, aes(PHE)) + xlab("Age phase (per 5 year from <45 to >75)") +
  theme(panel.background=element_blank()) + 
  xlab("Age group") + ylab("Log Odds Ratio") + 
  geom_line(aes(y = beta_univariate_1, colour = "univariate")) + 
  geom_ribbon(aes(ymin = beta_univariate_1 - 2*variance_univariate_1, ymax = beta_univariate_1 + 2*variance_univariate_1), fill = purple, alpha = 0.3) +
  geom_line(aes(y = beta_multivariate_1, colour = "multivariate")) + 
  geom_ribbon(aes(ymin = beta_multivariate_1 - 2*variance_multivariate_1, ymax = beta_multivariate_1 + 2*variance_multivariate_1), fill = green, alpha = 0.3) +
  geom_line(aes(y = true_curve, colour = "True Curve"), linetype = "dashed", size = 3) + ylim(c(0.04, 0.15)) + 
  scale_colour_manual(name="Inference pattern",values=c("univariate" = purple, "multivariate" = green, "True Curve" = blue))
plt

#############################################
# simulate with GxE and GxG, heterogeneity
#############################################
source("Genetic_longitudinal_functions.R")
N <- 50000
num_snp <- 50
num_buckets <- 8
beta_0 <- 0.1
variance_g <- 4
s <- 0
h_0 <- (9.5 + s * (-3.5:3.5))/10000

GxE_simualtion_rslt <- simulate_interacted_genetic_profiles_ph(N, num_snp,num_buckets,beta_0, variance_g , h_0)
para <- GxE_simualtion_rslt[[1]]

BETA <- GxE_simualtion_rslt[[2]][[4]]
z_prob <- GxE_simualtion_rslt[[2]][[2]]
sigma_i <- list()
for(j in 1:length(z_prob)){
  ZSigma <- mapply(function(x,y) x*y, z_prob[[j]],  para$sigmainverse,SIMPLIFY = FALSE)
  sum_sigmainv <- Reduce('+', ZSigma)
  Aj <- t(para$X) %*% sum_sigmainv %*% para$X
  sigma_i[[j]] <- para$X %*% solve(Aj) %*% t(para$X)
}

df_mean <- data.frame("PHE"= 1:8)
df_mean[paste("beta_", 1, sep = "")] = BETA[[1]]
df_mean[paste("variance_", 1, sep = "")] = sqrt(diag(sigma_i[[1]]))
df_mean["true_curve"] <- rep(para$beta_0, num_buckets)
  
plt <- ggplot(df_mean, aes(PHE)) + xlab("Age phase (per 5year from <45 to >75)") +
    ggtitle(paste0("GxE effect"))+
    theme(panel.background=element_blank()) + 
    xlab("Age group") + ylab("Log Odds Ratio") + 
    geom_line(aes(y = beta_1, colour = "Inferred profile")) + 
    geom_ribbon(aes(ymin = beta_1 - 2*variance_1, ymax = beta_1 + 2*variance_1), fill = red, alpha = 0.3) +
    geom_line(aes(y = true_curve, colour = "Central Effect"), linetype = "dashed", size = 3) + 
    scale_colour_manual(name="Type",values=c("Inferred profile" = red, "Central Effect" = blue))

ggsave("~/Desktop/Writting_up_genetic_risk/GxE.png",width = 6, height = 4,plt)

#########################
# threshold simulation
#########################
source("Genetic_longitudinal_functions.R")
# first condition, the starting point is changed by genetics
N <- 50000
allle_freq <- 0.3
genetics_population <- generate_genetics(allle_freq, N)
beta_0 <- 0.1
S0 <- -80 * exp(- beta_0 * genetics_population)
Drift <- rep(1, N)
Variance <- rep(4, N)

rslt_event <- generate_liability_disease(S0, Drift, Variance, censored = 60)

PATH <- rslt_event[[1]][,21:60] # only take after 20 years
event <- pmax(rslt_event[[2]] - 20, 0)

state <- 1*(event < 41)
time <- event

# get the survival curve and fit it
h <- rep(0,40) 
for(i in 1:39){
  h[i] <- sum(state*time > (i-1) & state*time <= i) / sum(time > (i-1)) 
}

num_buckets <- 8
Beta_surv <- matrix(0, 1, num_buckets)
SE_Surv <- matrix(0, 1,num_buckets) 

# get the list of all variables
genetics_population <- as.data.frame(genetics_population)

int_length <- 40/num_buckets

for(gp in 1:num_buckets){
  cases <- (time > int_length * (gp-1)) * (time <= int_length * (gp))*state
  controls <- (time > int_length * (gp)) 
  sample <- cbind(genetics_population,data_frame(time,cases,controls)) %>% 
    filter(cases == 1 | controls == 1 ) %>%
    mutate(time = pmin(time-int_length * (gp-1), int_length))
  mySurv <- coxph(as.formula(paste0("Surv(time, cases) ~ genetics_population" ) ), data = sample)
  SE_Surv[, gp] <- sqrt(diag(mySurv$var))
  Beta_surv[, gp] <- mySurv$coefficients
}

FULL_path <- rslt_event[[1]]
# get the genetics
non_allele <- FULL_path[which(genetics_population == 0)[4],,drop=FALSE]
one_allele <- FULL_path[which(genetics_population == 1)[1],,drop=FALSE]
two_allele <- FULL_path[which(genetics_population == 2)[1],,drop=FALSE]
mean_non_allele <- colMeans(non_allele)
mean_one_allele <- colMeans(one_allele)
mean_two_allele  <- colMeans(two_allele)

mean_path <- colMeans(FULL_path)

df_path <- data_frame(age = 1:60, mean_path = mean_path, mean_non_allele = mean_non_allele, mean_one_allele  = mean_one_allele , mean_two_allele  = mean_two_allele )
plt <- ggplot(df_path, aes(x = age)) + 
  theme(panel.background=element_blank()) + 
  geom_line(aes( y = mean_non_allele, color = "Homogeneous wild type"),size = 1) +
  geom_line(aes( y = mean_one_allele, color = "Heterogeneous"),size = 1) +
  geom_line(aes( y = mean_two_allele, color = "Homogeneous risk alleles"),size = 1) +
  geom_hline(yintercept=0 , linetype="dashed", color = red, size = 1) +
  scale_colour_manual(name="Genotype",values=c("Homogeneous wild type" = green, "Heterogeneous" = blue, "Homogeneous risk alleles" = red))
ggsave("~/Desktop/Writting_up_genetic_risk/Threshold_example.png", width = 7, height = 4,plt)

# plot the estimation
df_profile_liability <- data_frame(age = c(20,25,30,35,40,45,50,55), Beta = as.vector(Beta_surv), SE = as.vector(SE_Surv))
plt <- ggplot(df_profile_liability, aes(age)) + 
  geom_line(aes(y = Beta), color = red) +
  geom_ribbon(aes(ymin=Beta-1.96*SE, ymax=Beta+ 1.96*SE), fill= red, alpha=.3) + 
  theme(legend.position = "none",panel.background=element_blank()) + 
  xlab("Age group") + ylab("Log Odds Ratio")
ggsave("~/Desktop/Writting_up_genetic_risk/effect_size_liability.png", width = 4, height = 4,plt)
# 
#   geom_ribbon(aes( ymax = path_975, ymin = path_025), fill = "blue",alpha = 0.3) + 
#   geom_ribbon(aes( ymax = new_daily + 2*err_daily, ymin = new_daily - 2*err_daily), fill = "red",alpha = 0.3) +
#   geom_ribbon(aes( ymax = ab_new_daily + 2*err_ab, ymin = ab_new_daily- 2*err_ab), fill = "green",alpha = 0.3)

