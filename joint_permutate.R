######################################
# creating permutated dataset
######################################
# file joint_permutate.R
library(ggplot2)
require(dplyr)
library(stringr)
library(pdist)
require(survival)
library("mvtnorm")
# get the id
args <- commandArgs(trailingOnly = TRUE) 
# first argument is the permutation id
ID <- as.numeric(args[1])
# second is the disease id
ds_id <- args[2]
set.seed(19940110+ID)
# using iteration to avoid too many files
iter_rep <- 100

# read table
control <- read.table(paste(ds_id,"/",ds_id,"_control.txt",sep=""), sep="\t") %>%
  rename(eid = V1, time = V2, group = V3, censor = V4)
case <- read.table(paste(ds_id,"/",ds_id,"_case.txt",sep=""), sep="\t") %>%
  rename(eid = V1, time = V2, group = V3, censor = V4)

data_surv <- read.table(paste(ds_id,"/",ds_id,"_pheno_over_age.txt",sep=""), sep="\t") %>%
  rename(eid = V1, time = V2, group = V3, censor = V4)

# get the accumulated bucket size
bucket_sz <- rep(0,8)
for(i in 1:8){
  sz <- case %>% filter(group == i) %>% dim()
  bucket_sz[i] <- sz[1]
}
bucket_sz <- c(0, cumsum(bucket_sz))

# the case data are now permutated, it has to be permutated after computing bucket
for(itr in 1:iter_rep){
  case_idx <- sample(nrow(case), replace = FALSE)
  # the control index is the sampled index for control sample
  control_idx <- c(sapply(case_idx, function(x) 4*(x-1)+c(1,2,3,4)))
  
  case <- case[case_idx,]
  control <- control[control_idx,]
  
  pheno_data <- list()
  
  for(i in 1:8){
    pheno_data[[i]] <- case %>%
      slice((bucket_sz[i]+1):(bucket_sz[i+1])) %>%
      mutate(group = i)
    pheno_data[[i]] <- control %>% 
      slice((bucket_sz[i]*4+1):(bucket_sz[i+1]*4)) %>%
      mutate(group = i) %>%
      rbind(pheno_data[[i]]) 
  }
  pheno_data <- bind_rows(pheno_data)
  DIR <- paste(ds_id,"/permutation",ID, sep = "")
  write.table(pheno_data, paste(DIR, "/",ds_id,"_per_",ID, "_itr_",itr,"_pheno_age.txt",sep=""), sep="\t", col.names = FALSE, row.names = FALSE)
}

################################################
# do survival analysis on permutated dataset
################################################

data <- read.table("/users/mcvean/xilin/xilin/UK_biobank/covariates_for_asso_over_time.txt", sep="\t") %>%
  select(-V2,  -V3) %>%
  rename(eid = V1)
genetics <- read.table(paste0(ds_id,"/",ds_id,".raw"), header=T) %>%
  select(-IID, -PAT, -MAT, -SEX, -PHENOTYPE) 

SNPs <- names(genetics)[2:dim(genetics)[2]]
# get a PC list for covariates
PC_list <- "V4 "
for(pc in 5:43){
  PC_list <- paste0(PC_list , " + V", pc)
}
# creat a complete snp list
SNP_list <- SNPs[1]
for(snp in  SNPs[2:length(SNPs)]){
  SNP_list <- paste0(SNP_list, " + ", snp)
}

flip_snp <- read.csv(paste(ds_id,"/",ds_id,"_flipped_snp.csv", sep = ""))

# define a function to create the T*p matrix X
X_cb_spline <- function(P,t){
  X_base <- cbind(rep(1,t), 1:t, (1:t)^2, (1:t)^3)
  if(P <= 4){
    return(as.matrix(X_base[,1:P]))
  }else{
    X <- sapply(1+(1:(P-4))*(t-1)/(P-3), function(x) pmax(((1:t)-x)^3,0))
    return(cbind(X_base, X))
  }
}

# load("smoothness_d_f.Rdata")
# p <- best_p  %>% filter(coding == ds_id) %>% pull(2)
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

p <- 3
X <- X_cb_spline(p,8)

# function to compute the log likelihoood
comp_ll <- function(pi_saver, d_beta_hat){
  pi_d_beta_hat <- mapply(function(x,y) x * y, pi_saver, d_beta_hat, SIMPLIFY = FALSE)
  vector_for_sum <- sapply(1:S, function(j) Reduce("+", lapply(pi_d_beta_hat, function(x) x[j])))
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

# add prior sigma0 here 
comp_theta <- function(X, z_prob, betaj, sigmainverse, sigmabeta){
  p <- dim(X)[2]
  sigma0inv <- solve(matrix(.2, p, p) + (1-.2) * diag(p)) # add prior, super unefficient but good for now
  z_betaj_sum <- lapply(z_prob, function(z) rowSums(matrix(sapply(1:S, function(j) z[j] * t(X) %*% sigmabeta[[j]]), ncol = S)))
  z_sigma_sum <- lapply(z_prob, function(z) sigma0inv + Reduce("+", lapply(1:S, function(j) z[j] * t(X) %*% sigmainverse[[j]] %*% X)))
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

EM_K <- function(K, num_itr, flat_flag, frailty_flag){ # flat flag here refer to whether the first line need to be flat
  if(frailty_flag){
    load("adjustment.Rdata")
    adj <- adj[[as.character(ds_id)]]
  }
  # initialize parameter
  pi_saver <- list()
  theta <- list()
  BETA <- list()
  d_beta_hat <- list()
  for(j in 1:K){
    x <- runif(K)
    pi_saver[[j]] <- x[j]/sum(x)
    # pi_saver[[j]] <- 1/K
    theta[[j]] <-  matrix(rnorm(p,sd=.0001),p,1) # using a small sd value to avoid NA
    if(j == 1 && (flat_flag | frailty_flag)){
      # always force the first component to be a flat line
      BETA[[j]] <- matrix(0, M, 1)
    }else BETA[[j]] <- X %*% theta[[j]]
    d_beta_hat[[j]] <- sapply(1:S, function(x) 
      dmvnorm(t(betaj[[x]]), mean = BETA[[j]], sigma = sigmaSAMPLE[[x]]))
  }
  saver_ll <- matrix(0, num_itr, 1)
  for(itr in 1:num_itr){
    if(itr %% 50 ==0) print(itr)
    # compute ll
    ll <- comp_ll(pi_saver, d_beta_hat)
    if(is.na(ll)) break
    # assert_that(!is.na(ll),  msg = "NA values found!")
    saver_ll[itr,1] <- ll
    # M step
    z_prob <- comp_z(pi_saver, d_beta_hat)
    # E step
    theta <- comp_theta(X, z_prob, betaj, sigmainverse, sigmabeta)
    pi_saver <- comp_pi(z_prob)
    for(j in 1:K){
      if(j == 1 && flat_flag){
        # always force the first component to be a flat line
        z <- z_prob[[1]]
        beta0 <- sum(sapply(1:S, function(j) z[j] * sigmabeta[[j]]))/sum(sapply(1:S, function(j) z[j] * sigmainverse[[j]]))
        BETA[[j]] <- matrix(beta0, M, 1)
      }else if(j == 1 && frailty_flag){
        z <- z_prob[[1]]
        beta0 <- sum(sapply(1:S, function(j) z[j] *t(as.matrix(adj)) %*% sigmabeta[[j]]))/
          sum(sapply(1:S, function(j)  z[j] * t(as.matrix(adj)) %*% sigmainverse[[j]] %*%  as.matrix(adj)))
        BETA[[j]] <- beta0 * as.matrix(adj)
      }else BETA[[j]] <- X %*% theta[[j]]
      d_beta_hat[[j]] <- sapply(1:S, function(x) 
        dmvnorm(t(betaj[[x]]), mean = BETA[[j]], sigma = sigmaSAMPLE[[x]]))
    }
  }
  # plot(saver_ll, type = "l")
  return(list(pi_saver, z_prob, theta, BETA, saver_ll, ll))
}

K <- 6 # need a bit more replicate so I could actually identify multiple profiles 
permute_ll <- matrix(nrow = iter_rep, ncol = K)
for(itr in 1:iter_rep){
  DIR <- paste(ds_id,"/permutation",ID, sep = "")
  pheno_data <- read.table(paste(DIR, "/",ds_id,"_per_",ID, "_itr_",itr,"_pheno_age.txt",sep=""), sep="\t") %>%
    rename(eid = V1, time = V2, group = V3, censor = V4)
  data_coxph <- pheno_data %>%
    left_join(data, by = "eid") %>%
    left_join(genetics, by = c("eid" = "FID"))
  snp_for_clustering <- list()
  
  for(gp in 1:8){
    df <- data_coxph %>%
      filter(group == gp)
    # initiate with void values for OR and SE
    rslt <- data_frame(SNP = SNPs, OR = rep(exp(0),length(SNPs)), SE = rep(10,length(SNPs)), PHE = rep(gp,length(SNPs)))
    # need to handle numeric error
    try({
      mdl <- coxph(formula = as.formula(paste("Surv(time, censor) ~", SNP_list, " + ", PC_list)), data = df)
      rslt[, 2] <- exp(mdl$coefficients[1:length(SNPs)])
      rslt[, 3] <- sqrt(diag(mdl$var)[1:length(SNPs)])
        })
    # there might be different errors here, need  to use "try"
    snp_for_clustering[[gp]] <- rslt
  }
  snp_for_clustering <- bind_rows(snp_for_clustering)
  snp_for_clustering$SNP <- strsplit(snp_for_clustering$SNP, "_") %>% 
    sapply(function(x) x[1])
  
  useful_snp_data <- snp_for_clustering %>%
    # try to remove the numeric failure for being too close to no-effect
    mutate(SE = ifelse(abs(OR-1)< 10^(-5), 10, SE)) %>%
    mutate(OR = ifelse(SE > 10, 1, OR)) %>% 
    mutate(SE = ifelse(SE > 10, 10, SE)) %>% 
    # I need to remove the NAs
    mutate(SE = ifelse(is.na(SE), 10, SE)) %>%
    mutate(OR = ifelse(is.na(OR), 1, OR)) %>%
    mutate(OR = log(OR)) %>%
    # I will have to control the overflowing of coef (some -> inf)
    mutate(SE = ifelse(abs(OR) > 10, 10, SE)) %>%
    mutate(OR = ifelse(abs(OR) > 10, 0, OR)) %>% 
    mutate(OR = ifelse(SNP %in% flip_snp$SNP, -OR, OR))
  
  # save the permutation estiamte for debug purpose
  write.csv(useful_snp_data, paste(ds_id,"/permutation",ID, "/",ds_id, "_joint_per",ID,"_itr_",itr, "_useful_snp_data.csv", sep = ""), row.names = F)
  
  snp_lst <- useful_snp_data %>% 
    group_by(SNP)  %>%
    summarise(n()) %>%
    select(SNP) 
  sl <- snp_lst[["SNP"]]
  snp_lst <- as.list(sl)
  names(snp_lst) <- sl
  
  # not using sigmasnp
  sigmasnp <- 0.0004
  sigmaSAMPLE <- snp_lst %>%
    lapply(function(x) dplyr::filter(useful_snp_data, SNP == x)) %>%
    lapply(function(x) diag(x[["SE"]]^2) + sigmasnp)
  sigmainverse <- sigmaSAMPLE %>%
    lapply(function(x) solve(x))
  
  betaj <- snp_lst %>%
    lapply(function(x) dplyr::filter(useful_snp_data, SNP == x)) %>%
    lapply(function(x) matrix(x[["OR"]]))
  sigmabeta <- mapply(function(x,y) x %*% y, sigmainverse, betaj, SIMPLIFY = FALSE) 
  
  M <- length(betaj[[1]])
  S <- length(betaj)
  
  num_itr <- 3  #20
  ll <- matrix(0, K, num_itr)
  for(k in 1:K){
    for(rep in 1:num_itr){# 1241
      if(k == 1){
        rslt <- EM_K(k,100, FALSE, FALSE)
      }else if(k == 2){
        rslt <- EM_K(k,100, FALSE, FALSE)
      }else if(k == 3){
        rslt <- EM_K(k,100, FALSE, FALSE)
      }else if(k == 4){
        rslt <- EM_K(k,100, FALSE, FALSE)
      }else if(k==5){ # the single flat curve 
        rslt <- EM_K(1, 100, TRUE, FALSE)
      }else if(k==6){ # the frailty model
        rslt <- EM_K(1, 100, FALSE, TRUE)
      }
      ll[k, rep] <- rslt[[6]]
    }
  }
  permute_ll[itr,] <- apply(ll, 1, function(x) ifelse(anyNA(x), NA, x[which.max(x)]))
}
write.csv(permute_ll, paste(ds_id,"/permutation",ID, "/",ds_id,"_jointly_", ID, "_ll_matrix.csv", sep = ""), row.names = F)

####################################
# do the plottings 
####################################





