######################################
# survival analysis framework
######################################
require(survival)
require(dplyr)
library("mvtnorm")

args <- commandArgs(trailingOnly = TRUE)                                  
ds_id <- args[1]

data <- read.table("/users/mcvean/xilin/xilin/UK_biobank/covariates_for_asso_over_time.txt", sep="\t") %>%
  select(-V2,  -V3) %>%
  rename(eid = V1)
genetics <- read.table(paste0(ds_id,"/",ds_id,".raw"), header=T) %>%
  select(-IID, -PAT, -MAT, -SEX, -PHENOTYPE) 

data_coxph <- read.table(paste0(ds_id,"/",ds_id,"_pheno_over_age.txt"), header=F) %>%
  rename(eid = V1, time = V2, group = V3, censor = V4)

data_coxph <- data_coxph %>%
  left_join(data, by = "eid") %>%
  left_join(genetics, by = c("eid" = "FID"))

snp_for_clustering <- list()
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
for(gp in 1:8){
  df <- data_coxph %>%
    filter(group == gp)
  rslt <- data_frame(SNP = SNPs, OR = rep(0,length(SNPs)), SE = rep(0,length(SNPs)), PHE = rep(gp,length(SNPs)))
  mdl <- coxph(formula = as.formula(paste("Surv(time, censor) ~", SNP_list, " + ", PC_list)), data = df)
  rslt[, 2] <- exp(mdl$coefficients[1:length(SNPs)])
  rslt[, 3] <- sqrt(diag(mdl$var)[1:length(SNPs)])
  snp_for_clustering[[gp]] <- rslt
}
snp_for_clustering <- bind_rows(snp_for_clustering)
snp_for_clustering$SNP <- strsplit(snp_for_clustering$SNP, "_") %>% 
  sapply(function(x) x[1])

dir.create(paste(ds_id,"/joint_estimation", sep = ""))
snp_for_clustering <- snp_for_clustering %>%
  mutate(SE = ifelse(abs(OR-1)< 10^(-5), 10, SE)) %>%
  mutate(OR = ifelse(SE > 10, 1, OR)) %>% 
  mutate(SE = ifelse(SE > 10, 10, SE))
write.csv(snp_for_clustering, paste(ds_id,"/joint_estimation/",ds_id,"_multivariate_snp_for_clustering.csv", sep = ""), row.names = F)

###############################################
# doing EM
###############################################
flip_snp <- read.csv(paste(ds_id,"/",ds_id,"_flipped_snp.csv", sep = ""))
useful_snp_data <- snp_for_clustering %>% 
  # I need to remove the NAs
  mutate(SE = ifelse(is.na(SE), 10, SE)) %>%
  mutate(OR = ifelse(is.na(OR), 1, OR)) %>%
  mutate(OR = log(OR)) %>%
  # I will have to control the overflowing of coef (some -> inf)
  mutate(SE = ifelse(abs(OR) > 10, 10, SE)) %>%
  mutate(OR = ifelse(abs(OR) > 10, 0, OR)) %>% 
  mutate(OR = ifelse(SNP %in% flip_snp$SNP, -OR, OR))

source("Genetic_longitudinal_functions.R")

p <- 3
X <- X_cb_spline(p,8)
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
num_itr <- 20  #20
K <- 6
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
ll <- apply(ll, 1, function(x) x[which.max(x)])
write.csv(ll, paste(ds_id,"/joint_estimation/",ds_id, "_ll_matrix.csv", sep = ""), row.names = F)

list_of_rslt <- list()
for(k in 1:K){
  likelihood_thre <- ll[k] - .1
  likelihood <- -Inf
  while (likelihood < likelihood_thre) {
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
    likelihood <- rslt[[6]]
  }
  list_of_rslt[[k]] <- rslt
}

# Save an object to a file
save(list_of_rslt, file = paste(ds_id,"/joint_estimation/",ds_id, "_true_data.RData", sep = ""))

##########################
# plot profile inferred
##########################



