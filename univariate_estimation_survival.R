######################################
# survival analysis framework
######################################
require(survival)
require(dplyr)
library("mvtnorm")

args <- commandArgs(trailingOnly = TRUE)                                  
ds_id <- args[1]

##############################################################################################################
# be aware the section below need to be accomodated for the specific data source.
# The preprocessing for UK biobank dataset are not fully shown to comply with data sharing requirement of the UK Biobank 

# load the covariates for individuals
data <- read.table("/XXX/covariates_for_asso_over_time.txt", sep="\t") %>%
  select(-V2,  -V3) %>%
  rename(eid = V1)
# load getnetic data; the the genotype are encoded as in 0,1,2 for each locus;
# here we use plink/1.90b2  with flag --recode A to generate the genetic file for the subjects and SNPs of interest
genetics <- read.table(paste0(ds_id,"/",ds_id,".raw"), header=T) %>%
  select(-IID, -PAT, -MAT, -SEX, -PHENOTYPE) 

# load longitudinal data, which contain rows of events; the data frame has 
# has four columns: 1. id of the individual 2. event time, which is event point with respect to the beginning of interval; event 
# 3. group, an integer specify the age interval that cover the events 4. censor, which specify whether the event is disease (indicated by 1) or 
# a censoring event (indicated by 0). 
data_coxph <- read.table(paste0(ds_id,"/",ds_id,"_pheno_over_age.txt"), header=F) %>%
  rename(eid = V1, time = V2, group = V3, censor = V4)


# be aware the section above need to be accomodated for the specific data source.
# The section above should load longitudinal information, genotype information, and covaraites. 
##############################################################################################################

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
for(gp in 1:8){
  df <- data_coxph %>%
    filter(group == gp)
  rslt <- data_frame(SNP = SNPs, OR = rep(0,length(SNPs)), SE = rep(0,length(SNPs)), PHE = rep(gp,length(SNPs)))
  for(snp in SNPs){
    # I need to include all PCs as covariates
    mdl <- coxph(formula = as.formula(paste("Surv(time, censor) ~", snp, "+", PC_list)), data = df)
    # doing an adjustment for unobserved covariate here
    rslt[which(rslt$SNP == snp), 2] <- exp(mdl$coefficients[1])
    rslt[which(rslt$SNP == snp), 3] <- sqrt(mdl$var[1,1])
  }
  snp_for_clustering[[gp]] <- rslt
}
snp_for_clustering <- bind_rows(snp_for_clustering)
snp_for_clustering$SNP <- strsplit(snp_for_clustering$SNP, "_") %>% 
  sapply(function(x) x[1])

dir.create(paste(ds_id,"/univariate_estimation", sep = ""))
snp_for_clustering <- snp_for_clustering %>%
  mutate(SE = ifelse(abs(OR-1)< 10^(-5), 10, SE)) %>%
  mutate(OR = ifelse(SE > 10, 1, OR)) %>% 
  mutate(SE = ifelse(SE > 10, 10, SE))
write.csv(snp_for_clustering, paste(ds_id,"/result_of_true_data/",ds_id,"_snp_for_clustering.csv", sep = ""), row.names = F)

###############################################
# doing EM
###############################################
flip_snp <- read.csv(paste(ds_id,"/",ds_id,"_flipped_snp.csv", sep = ""))

# first get the convergence of EMs
para <- load_SNP_coef(paste("univariate_estimation_rm_LD/",ds_id,"_snp_for_clustering.csv", sep = ""), 
                      paste("SNP_info/",ds_id,"_flipped_snp.csv", sep = ""))

num_itr <- 2 #20
K <- 6
ll <- matrix(0, K, num_itr)
for(k in 1:K){
  for(rep in 1:num_itr){# 1241
    para$p <- 3
    para$X <- X_cb_spline(para$p,8)
    para$sigma0inv <- solve(diag(para$p) * 1) 
    if(k == 1){
      rslt <- EM_K(k,100, para)
    }else if(k == 2){
      rslt <- EM_K(k,100, para)
    }else if(k == 3){
      rslt <- EM_K(k,100, para)
    }else if(k == 4){
      rslt <- EM_K(k,100, para)
    }else if(k==5){ # the single flat curve 
      para$p <- 1
      para$X <- X_cb_spline(para$p,8)
      para$sigma0inv <- solve(diag(para$p) * 1) 
      rslt <- EM_K(1, 100, para)
    }else if(k==6){ # the frailty model
      para$p <- 1
      para$X <- X_cb_spline(para$p,8)
      load("adjustment.Rdata")
      adj <- adj[[as.character(ds_id)]]
      para$X[,1] <- adj
      para$sigma0inv <- solve(diag(para$p) * 1) 
      rslt <- EM_K(1, 100, para)
    }
    ll[k, rep] <- rslt[[6]]
  }
}
ll <- apply(ll, 1, function(x) x[which.max(x)])
write.csv(ll, paste(ds_id,"/result_of_true_data/",ds_id, "_ll_matrix.csv", sep = ""), row.names = F)

list_of_rslt <- list()
for(k in 1:K){
  likelihood_thre <- ll[k] - .1
  likelihood <- -Inf
  while (likelihood < likelihood_thre) {
    para$p <- 3
    para$X <- X_cb_spline(para$p,8)
    para$sigma0inv <- solve(diag(para$p) * 1) 
    if(k == 1){
      rslt <- EM_K(k,100, para)
    }else if(k == 2){
      rslt <- EM_K(k,100, para)
    }else if(k == 3){
      rslt <- EM_K(k,100, para)
    }else if(k == 4){
      rslt <- EM_K(k,100, para)
    }else if(k==5){ # the single flat curve 
      para$p <- 1
      para$X <- X_cb_spline(para$p,8)
      para$sigma0inv <- solve(diag(para$p) * 1) 
      rslt <- EM_K(1, 100, para)
    }else if(k==6){ # the frailty model
      para$p <- 1
      para$X <- X_cb_spline(para$p,8)
      load("adjustment.Rdata")
      adj <- adj[[as.character(ds_id)]]
      para$X[,1] <- adj
      para$sigma0inv <- solve(diag(para$p) * 1) 
      rslt <- EM_K(1, 100, para)
    }
    likelihood <- rslt[[6]]
  }
  list_of_rslt[[k]] <- rslt
}
# Save an object to a file
save(list_of_rslt, file = paste(ds_id,"/result_of_true_data/",ds_id, "_true_data.RData", sep = ""))
