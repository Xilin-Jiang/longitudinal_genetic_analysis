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

########################
# plot joint estimation 
########################
pt <- "*_snp_for_clustering.csv"
temp <- list.files(paste("Joint_estimation/", sep=""), pattern=pt)

ds_list <- lapply(temp, function(x) strsplit(x, "_")[[1]][1])

id <- 2

ds_id <- as.character(ds_list[[id]])

ds_id <- "Z955"
load(paste("./Joint_estimation/",ds_id,"_true_data.RData", sep = ""))


# determine which profile to show

rslt <- list_of_rslt[[2]]

BETA <- rslt[[4]]

z_prob <- rslt[[2]]
lapply(z_prob, function(x) sum(x>.5))

flip_snp <- read.csv(paste("SNP_info/",ds_id,"_flipped_snp.csv", sep = ""))
useful_snp_data <- read.csv(paste0("Joint_estimation/", ds_id, "_multivariate_snp_for_clustering.csv")) %>%
  mutate(OR = log(OR)) %>%
  mutate(OR = ifelse(SNP %in% flip_snp$SNP, -OR, OR))


X_cb_spline <- function(P,t){
  X_base <- cbind(rep(1,t), 1:t, (1:t)^2, (1:t)^3)
  if(P <= 4){
    return(as.matrix(X_base[,1:P]))
  }else{
    X <- sapply(1+(1:(P-4))*(t-1)/(P-3), function(x) pmax(((1:t)-x)^3,0))
    return(cbind(X_base, X))
  }
}
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

P <- .2
sigma0 <- matrix(P, p, p) + (1-P) * diag(p)
# we need sigma0 if no snp is assigned to it
sigma0inv <- solve(sigma0)
# find the allocation of curves

z_asign <- sapply(1:S, function(s) which.max(lapply(z_prob, function(x) x[s])))

sigma_i <- list()
for(j in 1:length(z_prob)){
  ZSigma <- mapply(function(x,y) x*y, z_prob[[j]],  sigmainverse,SIMPLIFY = FALSE)
  sum_sigmainv <- Reduce('+', ZSigma)
  Aj <- sigma0inv + t(X) %*% sum_sigmainv %*% X
  sigma_i[[j]] <- X %*% solve(Aj) %*% t(X)
}

df_mean <- data.frame("PHE"= 1:8)
for(i in 1:2){
  df_mean[paste("beta_", i, sep = "")] = BETA[[i]]
  df_mean[paste("variance_", i, sep = "")] = sqrt(diag(sigma_i[[i]]))
}

ggplot(df_mean, aes(PHE)) + xlab("Age phase (per 5
                                 year from <45 to >75)") +
  # ylim(.06, .125)+
  ggtitle(paste0("Disease ID:",ds_id))+
  theme(legend.position = "none",panel.background=element_blank()) + 
  xlab("Age group") + ylab("Log Odds Ratio") + 
  geom_line(aes(y = beta_1), colour = "red") + 
  geom_ribbon(aes(ymin = beta_1 - 2*variance_1, ymax = beta_1 + 2*variance_1), fill = "red", alpha = 0.3) +
  geom_line(aes(y = beta_2), colour = "blue") +
  geom_ribbon(aes(ymin = beta_2 - 2*variance_2, ymax = beta_2 +2*variance_2), fill = "blue", alpha = 0.3) +
  geom_line(aes(y = beta_3), colour = "green") +
  geom_ribbon(aes(ymin = beta_3 - 2*variance_3, ymax = beta_3 +2*variance_3), fill = "green", alpha = 0.3) +
  geom_line(aes(y = beta_4), colour = "orange") +
  geom_ribbon(aes(ymin = beta_4 - 2*variance_4, ymax = beta_4 +2*variance_4), fill = "orange", alpha = 0.3)




