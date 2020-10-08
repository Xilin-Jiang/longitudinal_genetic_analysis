library("mvtnorm")
library("dplyr")
library("ggplot2")
library("tidyr")
library("pdist")
library("assertthat")
library("stringr")
library("survival")
setwd("/Users/xilin/Desktop/genetics_longitudinal_data/longitudinal_data/")
# The palette with grey:
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

# The palette with black:
cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

# To use for fills, add
scale_fill_manual(values=cbPalette)
# To use for line and point colors, add  
scale_colour_manual(values=cbPalette)

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
X_cb_spline <- function(P,t){
  X_base <- cbind(rep(1,t), 1:t, (1:t)^2, (1:t)^3)
  if(P <= 4){
    return(as.matrix(X_base[,1:P]))
  }else{
    X <- sapply(1+(1:(P-4))*(t-1)/(P-3), function(x) pmax(((1:t)-x)^3,0))
    return(cbind(X_base, X))
  }
}
#####################################
# searching for the best smoothness
#####################################
# define a function to create the T*p matrix X of cubic spline
# when P == 1, it sould be a constant matrix, while P == 2,3 
# we will be adding the quadratic and cubic term
# the last term is going to be 
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
# get the likelihood gaining
setwd("~/Desktop/genetics_longitudinal_data/longitudinal_data/") # create a file for local testing:

pt <- "*_snp_for_clustering.csv"
temp <- list.files(paste("join_estimation_rm_LD/", sep=""), pattern=pt)

ds_list <- lapply(temp, function(x) strsplit(x, "_")[[1]][1])
icd_list <- data.frame(coding = unlist(ds_list)) %>%
  filter(!coding == "I258", !coding == "Z951",!coding == "Z955")
sv_ll <- matrix(nrow = dim(icd_list)[1], ncol = 10)
for(id in 1:dim(icd_list)[1]){
  ds_id <- icd_list$coding[id]
  flip_snp <- read.csv(paste("SNP_info/",ds_id,"_flipped_snp.csv", sep = ""))
  snp_for_clustering <- read.csv(paste0("join_estimation_rm_LD/",ds_id,"_multivariate_snp_for_clustering.csv"))
  useful_snp_data <- snp_for_clustering %>% 
    # I need to remove the NAs
    mutate(SE = ifelse(is.na(SE), 10, SE)) %>%
    mutate(OR = ifelse(is.na(OR), 1, OR)) %>%
    mutate(OR = log(OR)) %>%
    # I will have to control the overflowing of coef (some -> inf)
    mutate(SE = ifelse(abs(OR) > 10, 10, SE)) %>%
    mutate(OR = ifelse(abs(OR) > 10, 0, OR)) %>% 
    mutate(OR = ifelse(SNP %in% flip_snp$SNP, -OR, OR))
  
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
  
  rslt <- list()
  ll <- rep(0,10)
  for(p in 1:6){
    X <- X_cb_spline(p,8)
    rslt[[p]] <- EM_K(1,100, FALSE, FALSE)
    ll[p] <- rslt[[p]][[6]]
  }
  # add the natural cubic spline
  for(p in 3:4){
    X <- X_nat_cb(p,8)
    rslt[[p+4]] <- EM_K(1,100, FALSE, FALSE)
    ll[p+4] <- rslt[[p+4]][[6]]
  }
  # do fitting with a fixed frailty
  p <- 3
  X <- X_nat_cb(3,8)
  rslt[[9]] <- EM_K(1,100, FALSE, TRUE)
  ll[9] <- rslt[[9]][[6]]
  # do the analysis, with a cubic spline but the first spline factor to be frailed
  load("adjustment.Rdata")
  adj <- adj[[ds_id]]
  X[,1] <- adj
  rslt[[10]] <- EM_K(1,100, FALSE, FALSE)
  ll[10] <- rslt[[10]][[6]]
  sv_ll[id,] <- ll
}
ll_gain <- sv_ll[,2:6] - sv_ll[,1]
ll_gain <- sapply(1:5, function(x) -log10(1-pchisq(2*ll_gain[,x], df=x)))
rowMeans(sapply(1:dim(icd_list)[1], function(x) rank(ll_gain[x,]))) # get the rank for different d.f.
# this part is considering the natrual cubic spline
ll_nat_cb_gain <- sv_ll[,7:8] - sv_ll[,1]
ll_nat_cb_gain <- sapply(2:3, function(x) -log10(1-pchisq(2*ll_nat_cb_gain[,x-1], df=x)))
ll_gain <- cbind(ll_gain, ll_nat_cb_gain)
best_p <- sapply(1:dim(icd_list)[1], function(x) which.max(ll_gain[x,])) + 1
best_p <- data.frame(icd_list, best_p)
save(best_p, file = "smoothness_d_f.Rdata")
load("smoothness_d_f.Rdata")
p <- best_p  %>% filter(coding == ds_id) %>% pull(2)
plt_ll_gain <- data.frame(icd_list, ll_gain)

plt <- ggplot(plt_ll_gain, aes(x = coding, y = X1)) + 
  geom_point(size = 3, aes(color = "linear")) + 
  theme(legend.position = "top",panel.background=element_blank()) + 
  ylab("-log10(P)") +
  # scale_y_continuous(trans='log10') + # using log scale
  xlab("Disease (ICD 10 code)") +
  geom_point(data = plt_ll_gain, aes(x = coding, y = X2, color = "quadratic"), shape = 19, size = 3) +
  geom_hline(yintercept=-log10(0.05), linetype="dashed", color = "red") + 
  geom_point(data = plt_ll_gain, aes(x = coding, y = X3, color = "cubic"), shape = 19, size = 3) +
  geom_point(data = plt_ll_gain, aes(x = coding, y = X4, color = "one knot"), shape = 19, size = 3) +
  geom_point(data = plt_ll_gain, aes(x = coding, y = X5, color = "two knot"), shape = 19, size = 3) + 
  scale_color_manual(name="Model Type", values=c("linear" = blue, "quadratic"  = red, "cubic" = orange, "one knot" = green, "two knot" = grey))
  
ggsave("~/Desktop/Writting_up_genetic_risk/Choose_smoothness.png",plt,width = 8, height =5)
# testing effect other than frailty
ll_test_frailty <- sv_ll[,c(10,7)] - sv_ll[,c(9,1)]
ll_test_frailty <- -log10(1-pchisq(2*ll_test_frailty, df=2))
plt_ll_frailty <- data.frame(icd_list, ll_test_frailty) 
ggplot(plt_ll_frailty, aes(x = coding, y = X1)) + 
  geom_point(size = 3, color = "#0072B2") + 
  theme(legend.position = "none",panel.background=element_blank()) + 
  ylab("-log_10(Likelihood gain)") +
  # scale_y_continuous(trans='log10') + # using log scale
  xlab("Disease (ICD 10 code)") +
  geom_point(data = plt_ll_frailty, aes(x = coding, y = X2), color = "#D55E00", shape = 8, size = 5) +
  geom_hline(yintercept=-log10(0.05), linetype="dashed", color = "#D55E00") 

#########################
# plot the curve
#########################
setwd("/Users/xilin/Desktop/genetics_longitudinal_data/longitudinal_data/")
# load the disease list
pt <- "*_snp_for_clustering.csv"
temp <- list.files(paste("Joint_estimation/", sep=""), pattern=pt)
icd_code <- lapply(temp, function(x) strsplit(x, "_")[[1]][1])
icd_list <- data.frame(coding = unlist(icd_code)) %>%
  filter(! coding == "I258", ! coding == "Z955", ! coding == "Z951")
plt <-list()

# local distribution
prevelance <- read.table("HES_prevalence_by_code.txt", 
                         sep = "\t", quote = "", header  = TRUE)
snp_assoc <- read.table("SNP_Assocs_per_code.txt", 
                        sep = "\t", quote = "", header  = TRUE)
snp_assoc <- snp_assoc %>%
  right_join(prevelance, by = c("coding" = "Code"))
snp_of_interest <- snp_assoc %>% 
  filter(NUMSNPASSOCS > 20, Prevalence > .5) %>% # using prevalence > .5 will get 6 more interesting diseases
  select(coding, NUMSNPASSOCS, meaning)
one_cluster <- T # FALSE

# also save the posterior and standard error
mean_profile <- matrix(nrow = dim(icd_list)[1], ncol = 8)
se_profile <- matrix(nrow = dim(icd_list)[1], ncol = 8)
for(id in 1:dim(icd_list)[1]){
  ds_id <- icd_list$coding[id]
  load(paste("./join_estimation_rm_LD/",ds_id,"_true_data.RData", sep = ""))
  
  # determine which profile to show
  load("smoothness_d_f.Rdata")
  # p <- best_p  %>% filter(coding == ds_id) %>% pull(2)
  p <- 3
  X <- X_cb_spline(p,8)
  if(one_cluster == TRUE){
    rslt <- list_of_rslt[[1]]
  }else{
    rslt <- list_of_rslt[[2]]
  }
  BETA <- rslt[[4]]
  z_prob <- rslt[[2]]
  flip_snp <- read.csv(paste("SNP_info/",ds_id,"_flipped_snp.csv", sep = ""))
  useful_snp_data <- read.csv(paste0("join_estimation_rm_LD/", ds_id, "_multivariate_snp_for_clustering.csv")) %>%
    mutate(OR = log(OR)) %>%
    mutate(OR = ifelse(SNP %in% flip_snp$SNP, -OR, OR))
  
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
  
  df_mean <- data.frame("PHE"= c(45,50,55,60,65,70,75,80))
  # get the meaning of the code
  title_meaning <- snp_of_interest %>% filter(coding == ds_id) %>% pull(3) %>% as.character()
  
  if(one_cluster == TRUE){
    df_mean[paste("beta_", 1, sep = "")] = BETA[[1]]
    df_mean[paste("variance_", 1, sep = "")] = sqrt(diag(sigma_i[[1]]))
    
    mean_profile[id,] <- BETA[[1]]
    se_profile[id,] <- sqrt(diag(sigma_i[[1]]))
    
    plt[[id]] <- ggplot(df_mean, aes(PHE)) + xlab("Age phase (per 5
                                                  year from <45 to >75)") +
      # ylim(.06, .125)+
      ggtitle(paste0("Disease ID:", title_meaning))+
      theme(legend.position = "none",panel.background=element_blank(),plot.title = element_text(size = 10, face = "bold")) + 
      xlab("Age (years)") + ylab("Log Odds Ratio") + 
      geom_line(aes(y = beta_1), colour = red) + 
      geom_ribbon(aes(ymin = beta_1 - 2*variance_1, ymax = beta_1 + 2*variance_1), fill = "#D55E00", alpha = 0.3) 
    ggsave(paste0("~/Desktop/Writting_up_genetic_risk/supplementary/multivariate_single_profile_",  # multivariate_single_profile_", 
                  as.character(ds_id), ".png"), width = 5, height = 4,plt[[id]])
  }else{
    print(ds_id)
    print(lapply(z_prob, function(x) sum(x>.5)))
    for(i in 1:2){
      df_mean[paste("beta_", i, sep = "")] = BETA[[i]]
      df_mean[paste("variance_", i, sep = "")] = sqrt(diag(sigma_i[[i]]))
    }
    
    plt[[id]] <- ggplot(df_mean, aes(PHE)) + xlab("Age phase (per 5
                                                  year from <45 to >75)") +
      # ylim(.06, .125)+
      ggtitle(paste0("Disease ID:", title_meaning))+
      theme(legend.position = "none",panel.background=element_blank(),plot.title = element_text(size = 10, face = "bold"))  + 
      xlab("Age (years)") + ylab("Log Odds Ratio") + 
      geom_line(aes(y = beta_1), colour = blue) + 
      geom_ribbon(aes(ymin = beta_1 - 2*variance_1, ymax = beta_1 + 2*variance_1),fill=blue, alpha = 0.3) +
      geom_line(aes(y = beta_2), colour = red) +
      geom_ribbon(aes(ymin = beta_2 - 2*variance_2, ymax = beta_2 +2*variance_2),fill=red, alpha = 0.3) 
      
    ggsave(paste0("~/Desktop/Writting_up_genetic_risk/supplementary/multivariate_two_profile_", 
                  as.character(ds_id), ".png"), width = 5, height = 4,plt[[id]])
  }
}
pasted_rslt <- matrix(mapply(function(x,y) paste0(as.character(x), " (SE = ", as.character(y), ")"), 
                             round(mean_profile,digits = 3), round(se_profile, digits = 3)), nrow = dim(icd_list)[1],ncol = 8)
profiles_df <- data.frame(icd_list, pasted_rslt)
write.csv(profiles_df,"~/Desktop/Writting_up_genetic_risk/Figures_supplementary/profiles_all_ds.csv")
#########################################
# testing for frailty, likelihood ratio
#########################################
setwd("/Users/xilin/Desktop/genetics_longitudinal_data/longitudinal_data/")
source("Genetic_longitudinal_functions.R")
# load the disease list
pt <- "*_snp_for_clustering.csv"
temp <- list.files(paste("Joint_estimation/", sep=""), pattern=pt)
icd_code <- lapply(temp, function(x) strsplit(x, "_")[[1]][1])
icd_list <- data.frame(coding = unlist(icd_code)) %>%
  filter(! coding == "I258", ! coding == "Z951", ! coding == "Z955")
all_ds_test <- matrix(nrow = dim(icd_list)[1], ncol = 2)
# fitting either a joint manner or univariate
joint_flag <- F
for(id in 1:dim(icd_list)[1]){
  ds_id <- icd_list$coding[id]
  # load the likelihood
  if(joint_flag){
    true_likelihood <- read.csv(paste("./join_estimation_rm_LD/",ds_id,"_ll_matrix.csv", sep = ""))
  }else{
    true_likelihood <- read.csv(paste("./univariate_estimation_rm_LD/",ds_id,"_ll_matrix.csv", sep = ""))
  }
  # first get the p for using one flexible line again one constant line
  true_ll_gain <- true_likelihood$x[1] - true_likelihood$x[5]
  ###############################
  ###############################
  flip_file <- paste("SNP_info/",ds_id,"_flipped_snp.csv", sep = "")
  # I will be using the univariate estimation which is more robust
  if(joint_flag){
    snp_file <- paste0("Joint_estimation/",ds_id,"_multivariate_snp_for_clustering.csv")
  }else{
    snp_file <- paste0("univariate_estimation_rm_LD/",ds_id,"_snp_for_clustering.csv")
  }
  para <- load_SNP_coef(SNP_FILE = snp_file, FLIP_FILE = flip_file)
  load("adjustment.Rdata")
  adj <- adj[[as.character(ds_id)]]
  para$p <- 1
  para$X <- X_cb_spline(para$p,8)
  para$sigma0inv <- solve(diag(para$p) * 1) 
  para$X[,1] <- adj
  rslt_1 <- EM_K(1,100, para)
  # do the analysis, with a cubic spline but the first spline factor to be frailty (intercept)
  ###############################################################################
  # here we change X to make sure it is properly nested for the likelihood ratio test
  ###############################################################################
  para$p <- 3
  para$X <- X_cb_spline(para$p,8)
  para$sigma0inv <- solve(diag(para$p) * 1) 
  para$X[,1] <- adj
  rslt_2 <- EM_K(1,100, para)
  true_ll_gain <- c(true_ll_gain, rslt_2[[6]] - rslt_1[[6]])
  ###############################
  ###############################
  p_value <- 1- c(pchisq(2*true_ll_gain, df=2))
  # then compute all the test against a single flat line
  all_ds_test[id,] <- p_value
}
fdr1 <- p.adjust(all_ds_test[,1],method = "fdr") 
fdr2 <- p.adjust(all_ds_test[,2],method = "fdr")
fdr_list <- icd_list
fdr_list$test_against_flat <- -log10(fdr1)
fdr_list$test_against_frailty <- -log10(fdr2)
plt <- ggplot(fdr_list, aes(x = coding, y = test_against_flat)) + 
  geom_point(size = 3, aes(color = "Deviation from uniformity", shape = "Deviation from uniformity")) + 
  theme(legend.position = "bottom",panel.background=element_blank()) + 
  ylab("-Log10Q (FDR adjusted)") +
  # scale_y_continuous() + # using log scale
  xlab("Disease (ICD 10 code)") +
  geom_hline(yintercept=-log10(0.05), linetype="dashed", color = red) +
  geom_point(data = fdr_list, aes(x = coding, y = test_against_frailty, color = "Deviation from frailty",shape = "Deviation from frailty"), size = 5)+ 
  scale_color_manual(name = "Test",values = c( "Deviation from uniformity" = blue, "Deviation from frailty" = red))+
  scale_shape_manual(name = "Test", values = c( "Deviation from uniformity" = 16, "Deviation from frailty" = 8))
ggsave("~/Desktop/Writting_up_genetic_risk/univariate_frailty_p_value_adjusted.png",plt,width = 8, height =5) 
#############################################################
# plotting 2: get the P-value plot!
#############################################################
setwd("/Users/xilin/Desktop/genetics_longitudinal_data/longitudinal_data/")
# load the disease list
pt <- "*_snp_for_clustering.csv"
temp <- list.files(paste("join_estimation_rm_LD/", sep=""), pattern=pt)
icd_code <- lapply(temp, function(x) strsplit(x, "_")[[1]][1])
icd_list <- data.frame(coding = unlist(icd_code)) %>%
  filter(! coding == "I258", ! coding == "Z955", ! coding == "Z951")

P_all_ds <- matrix(0,dim(icd_list)[1], 4) # the test is, 1) against flat 2) two against 1, 3) three agains two 4) four against three 5) against frailty
cnt <- matrix(nrow = dim(icd_list)[1])
for(id in 1:dim(icd_list)[1]){
  ds_id <- icd_list$coding[id]
  temp <- list.files("./joint_permutation_LD_rm/", pattern=as.character(ds_id))
  # need to deal with potential NAs 
  permute_ll <- list()
  cnt[id] <- length(temp) * 100
  for(i in 1:length(temp)){
    ll <- read.csv(paste("./joint_permutation_LD_rm/", temp[i], sep = ""))
    permute_ll[[i]] <- ll
  }
  permute_ll <- bind_rows(permute_ll)
  if(anyNA(permute_ll)){
    cnt[id] <- (cnt[id] - sum(is.na(permute_ll))/6)
    permute_ll[is.na(permute_ll)] <- -Inf
  }
  # print for debugging
  print(c(ds_id, cnt[id]))
  per_ll_gain <- permute_ll[,c(1,2,3,4,1)] - permute_ll[,c(5,1,2,3,6)]
  # remove the NA value from Inf - Inf
  per_ll_gain[is.na(per_ll_gain)] <- -Inf
  # load the likelihood
  true_likelihood <- read.csv(paste("./join_estimation_rm_LD/",ds_id,"_ll_matrix.csv", sep = ""))
  # first get the p for using two flexible line again one flexible line
  true_ll_gain <- true_likelihood$x[c(1,2,3,4)] - true_likelihood$x[c(5,1,2,3)]
  p_value <- sapply(1:length(true_ll_gain), function(s) sum(per_ll_gain[,s] > true_ll_gain[s])/cnt[id])
  P_all_ds[id,] <- p_value
}
# # this is likelihood ratio test
# for(id in 1:dim(icd_list)[1]){
#   ds_id <- icd_list$coding[id]
#   # load the likelihood
#   true_likelihood <- read.csv(paste("./Joint_estimation/",ds_id,"_ll_matrix.csv", sep = ""))
#   # first get the p for using two flexible line again one flexible line
#   true_ll_gain <- true_likelihood$x[c(1,2,3,4)] - true_likelihood$x[c(5)]
#   p_value <- 1- c(pchisq(2*true_ll_gain[1], df=2),pchisq(2*true_ll_gain[2], df=6), 
#                                                          pchisq(2*true_ll_gain[3], df=10),
#                                                         pchisq(2*true_ll_gain[4], df=14))
#   # then compute all the test against a single flat line
#   P_all_ds[id,] <- p_value
# }
# compute a FDR threshold, adjusted
fdr1 <- p.adjust(P_all_ds[,1],method = "fdr")  + 1/cnt
fdr2 <- p.adjust(P_all_ds[,2],method = "fdr")  + 1/cnt
fdr3 <- p.adjust(P_all_ds[,3],method = "fdr")  + 1/cnt
fdr4 <- p.adjust(P_all_ds[,4],method = "fdr")  + 1/cnt
# fdr5 <- p.adjust(P_all_ds[,5],method = "fdr") + 1/cnt
fdr_list <- icd_list
fdr_list$test_against_flat <- -log10(fdr1)
fdr_list$test_two_flexible_against_one_flexible <- -log10(fdr2)
fdr_list$test_three_against_two <- -log10(fdr3)
fdr_list$test_four_against_three <- -log10(fdr4) 
# fdr_list$test_against_frailty <- -log10(fdr5)

plt <- ggplot(fdr_list, aes(x = coding)) + 
  geom_point(size = 3, aes(y = test_against_flat,color = "Single Non-constant against constant", shape = "Single Non-constant against constant")) + 
  theme(panel.background=element_blank()) + 
  ylab("-Log10Q (FDR adjusted)") +
  # scale_y_continuous(trans='log10') + # using log scale
  xlab("Disease (ICD 10 code)") +
  geom_point(aes(x = coding, y = test_two_flexible_against_one_flexible, color = "Two clusters against one cluster", shape = "Two clusters against one cluster"), size = 3) +
  geom_hline(yintercept=-log10(0.1), linetype="dashed", color = red) +
  geom_point(aes(x = coding, y = test_three_against_two, color = "Three clusters against two clusters", shape = "Three clusters against two clusters"), size = 3) +
  geom_point(aes(x = coding, y = test_four_against_three, color = "Four clusters against three clusters", shape = "Four clusters against three clusters"), size = 3) +
  scale_color_manual(name="Test",values=c("Single Non-constant against constant" = blue, "Two clusters against one cluster" = red,"Three clusters against two clusters" = orange, "Four clusters against three clusters" = green))+
  scale_shape_manual(name = "Test", values = c("Single Non-constant against constant"= 19, "Two clusters against one cluster"= 18, "Three clusters against two clusters" = 19,"Four clusters against three clusters" = 18))
  # geom_point(data = fdr_list, aes(x = coding, y = test_against_frailty), color = "#D55E00", shape = 8, size = 5) 
ggsave("~/Desktop/Writting_up_genetic_risk/FDR_joint.png",plt,width = 10, height =5) 

# save a file for the table in paper
prevelance <- read.table("HES_prevalence_by_code.txt", 
                         sep = "\t", quote = "", header  = TRUE)
snp_assoc <- read.table("SNP_Assocs_per_code.txt", 
                        sep = "\t", quote = "", header  = TRUE)
snp_assoc <- snp_assoc %>%
  right_join(prevelance, by = c("coding" = "Code"))
snp_of_interest <- snp_assoc %>% 
  filter(NUMSNPASSOCS > 20, Prevalence > .5) %>% # using prevalence > .5 will get 6 more interesting diseases
  select(coding, NUMSNPASSOCS, meaning)
icd_list <- snp_of_interest %>% select(coding)
keep <- read.table(file="~/Desktop/genetics_longitudinal_data/longitudinal_data/keep.txt", 
                   header=FALSE, sep=" ") 
sz <- dim(keep)[1]
onset_age <- matrix(nrow = length(icd_list[[1]]), ncol = 2)
for(id in 1:length(icd_list[[1]])){
  ds_id <- icd_list[[1]][id]
  diag_age <- read.csv(paste0("gwas_25/",ds_id, "_diag_age.csv"))
  onset_age[id, 1] <- diag_age %>% 
    summarise(mean(age_diag)) %>%
    pull(1)
  onset_age[id, 2] <- diag_age %>% 
    summarise(n()) %>%
    pull(1)/sz
}
snp_of_interest <- data.frame(snp_of_interest,onset_age) %>%
  rename(onset_age = X1, prevalence = X2)

setwd("/Users/xilin/Desktop/genetics_longitudinal_data/longitudinal_data/")
# load the disease list
pt <- "*_snp_for_clustering.csv"
temp <- list.files(paste("Joint_estimation/", sep=""), pattern=pt)
icd_code <- lapply(temp, function(x) strsplit(x, "_")[[1]][1])
icd_list <- data.frame(coding = unlist(icd_code)) %>%
  filter(! coding == "I258", ! coding == "Z955", ! coding == "Z951")

q_list <- data.frame(icd_list, fdr1, fdr2, P_all_ds[,1]+ 1/cnt, P_all_ds[,2]+ 1/cnt) %>%
  rename(q1 = fdr1, q2 = fdr2, p1 = P_all_ds...1....1.cnt, p2 = P_all_ds...2....1.cnt)

tb_all <- q_list %>%
  left_join(snp_of_interest, by = "coding") %>%
  select(coding, prevalence, NUMSNPASSOCS, onset_age, meaning, q1, q2, p1, p2)
write.csv(tb_all, "~/Desktop/Writting_up_genetic_risk/Non_LD_indep_multivariate_q_value.csv")

# plot p value before adjustment
p_list <- icd_list
p_list$test_against_flat <- -log10(P_all_ds[,1]+ 1/cnt)
p_list$test_two_flexible_against_one_flexible <- -log10(P_all_ds[,2]+ 1/cnt)
p_list$test_three_against_two <- -log10(P_all_ds[,3]+ 1/cnt)
p_list$test_four_against_three <- -log10(P_all_ds[,4]+ 1/cnt) 
# fdr_list$test_against_frailty <- -log10(fdr5)

ggplot(p_list, aes(x = coding, y = test_against_flat)) + 
  geom_point(size = 3, color = "#0072B2") + 
  theme(legend.position = "none",panel.background=element_blank()) + 
  ylab("log10(P)") +
  # scale_y_continuous(trans='log10') + # using log scale
  xlab("Disease (ICD 10 code)") +
  geom_point(data = p_list, aes(x = coding, y = test_two_flexible_against_one_flexible), color = "#D55E00", shape = 18, size = 3) +
  geom_hline(yintercept=-log10(0.05), linetype="dashed", color = "#D55E00") +
  geom_point(data = p_list, aes(x = coding, y = test_three_against_two), color = "#E69F00", shape = 19, size = 3) +
  geom_point(data = p_list, aes(x = coding, y = test_four_against_three), color = "#009E73", shape = 18, size = 3) # +


####################################### 
# chooce the number of knots for spline 

#####################################
# plot frailty against true inference
#####################################
# compare how the decreasing pattern is explained by the unobserved coefficient model
pt <- "*_snp_for_clustering.csv"
temp <- list.files(paste("univariate_estimation_rm_LD/", sep=""), pattern=pt)

ds_list <- lapply(temp, function(x) strsplit(x, "_")[[1]][1])

# in the plot, also save the disease name for convenience
prevelance <- read.table("HES_prevalence_by_code.txt", 
                         sep = "\t", quote = "", header  = TRUE)
snp_assoc <- read.table("SNP_Assocs_per_code.txt", 
                        sep = "\t", quote = "", header  = TRUE)
snp_assoc <- snp_assoc %>%
  right_join(prevelance, by = c("coding" = "Code"))
snp_of_interest <- snp_assoc %>% 
  filter(NUMSNPASSOCS > 20, Prevalence > .5) %>% # using prevalence > .5 will get 6 more interesting diseases
  select(coding, NUMSNPASSOCS, meaning)

plt <- list()
for(id in 1:length(ds_list)){
  ds_id <- as.character(ds_list[[id]])
  load(paste("./univariate_estimation_rm_LD/",ds_id,"_true_data.RData", sep = ""))
  load("adjustment.Rdata")
  adj <- adj[[ds_id]]
  # determine which profile to show
  
  rslt <- list_of_rslt[[1]]
  
  BETA <- rslt[[4]]
  
  z_prob <- rslt[[2]]
  lapply(z_prob, function(x) sum(x>.5))
  
  flip_snp <- read.csv(paste("SNP_info/",ds_id,"_flipped_snp.csv", sep = ""))
  useful_snp_data <- read.csv(paste0("univariate_estimation_rm_LD/", ds_id, "_snp_for_clustering.csv")) %>%
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
  
  # compute the adjustment for the unobserved effect
  beta0 <- sum(sapply(1:S, function(j) t(as.matrix(adj)) %*% sigmabeta[[j]]))/sum(sapply(1:S, function(j)  t(as.matrix(adj)) %*% sigmainverse[[j]] %*%  as.matrix(adj)))
  
  df_mean <- data.frame("PHE"= 1:8)
  for(i in 1:1){
    df_mean[paste("beta_", i, sep = "")] = BETA[[i]]
    df_mean[paste("variance_", i, sep = "")] = sqrt(diag(sigma_i[[i]]))
  }
  df_mean["unobserve"] = beta0 * as.matrix(adj)
  
  # get the meaning of the code
  title_meaning <- snp_of_interest %>% filter(coding == ds_id) %>% pull(3) %>% as.character()
  
  plt[[ds_id]] <- ggplot(df_mean, aes(PHE)) + xlab("Age phase (per 5
                                                   year from <45 to >75)") +
    ggtitle(paste0("Disease ID:", title_meaning))+
    xlab("Age group") + ylab("Effect Size") + 
    geom_line(aes(y = beta_1), colour = "#D55E00") + 
    geom_ribbon(aes(ymin = beta_1 - 2*variance_1, ymax = beta_1 + 2*variance_1), fill = "#D55E00", alpha = 0.3) +
    geom_line(aes(y = unobserve), colour = "#0072B2", size = 2,linetype = "dashed" )+
    theme(legend.position = c(0.95, 0.95),
          legend.justification = c("right", "top"),panel.background=element_blank(),
          plot.title = element_text(size = 10))
    
}
for(id in 1:length(ds_list)){
  ggsave(paste0("~/Desktop/Writting_up_genetic_risk/supplementary/frailty_compare_to_risk_", 
                as.character(ds_list[id]), ".png"), width = 5, height = 4,plt[[id]])
}

########################################################################
# checking if there are interaction between genetics 2020-05-18
########################################################################
library(Hmisc)
library(reshape2)
df <- load("snp_data_291019.rdata")

pt <- "*.raw"
temp <- list.files(paste("Survival/", sep=""), pattern=pt)

ds_list <- lapply(temp, function(x) strsplit(x, "[.]")[[1]][1])
for(id in 1:length(ds_list)){
  ds_id <- ds_list[[id]]
  print(ds_id)
  genetics <- read.table(paste0("Survival/",ds_id,".raw"), header=T, check.names=FALSE) %>%
    select(-IID, -PAT, -MAT, -SEX, -PHENOTYPE)
  res <- rcorr(as.matrix(select(genetics, - FID)), type = "pearson")
  longData <- melt(res[[1]])
  longData<-longData[longData$value!=0,]
  longData<-longData[longData$value!=1,]
  cor_data <- longData %>% filter(value^2 > 0.04)
  snp_lst <- names(genetics)[2:length(genetics)]
  rm_snp <- c()
  for(snp in snp_lst){
    cor <- cor_data %>%
      filter(Var1 == snp| Var2 == snp)
    if(dim(cor)[1] > 0){
      rm_snp <- c(rm_snp,snp)
      cor_data <- cor_data %>%
        filter(!Var1 == snp & !Var2 == snp)
    }
  }
  rm_snp <- sapply(rm_snp, function(x) strsplit(x, "_")[[1]][1])
  pos <- which(t$coding == ds_id)
  t[pos,3] <- t[pos,3] - length(rm_snp)
  tables[[pos]] <- tables[[pos]] %>%
    filter(! SNP %in% rm_snp)
}
save(t,tables,file = "SNP_data_LD_rm_20200521.RData")

new_Assoc <- t$NUMSNPASSOCS

df <- load("snp_data_291019.rdata")
t$new_assoc <- new_Assoc
t %>% mutate(rm_assoc = NUMSNPASSOCS - new_assoc) 
# genetics <- read.table(paste0(ds_id,"/",ds_id,".raw"), header=T, check.names=FALSE) %>%
#   select(-IID, -PAT, -MAT, -SEX, -PHENOTYPE) 

res <- rcorr(as.matrix(select(genetics, - FID)), type = "pearson")
logP_matrix <- -log10(res[[3]])
longData <- melt(res[[1]])
longData<-longData[longData$value!=0,]
longData<-longData[longData$value!=1,]
ggplot(longData, aes(x = Var2, y = Var1)) +
  geom_raster(aes(fill=value)) +
  scale_fill_gradient(low="grey90", high="red") +
  labs(x="letters", y="LETTERS", title="Matrix") +
  theme_bw() + theme(axis.text.x=element_text(size=9, angle=0, vjust=0.3),
                     axis.text.y=element_text(size=9),
                     plot.title=element_text(size=11))

##################################
# create a SNP list to be removed 
##################################
cor_data <- longData %>% filter(abs(value) > 0.2)
snp_lst <- names(genetics)[2:length(genetics)]
rm_snp <- c()
for(snp in snp_lst){
  cor <- cor_data %>%
    filter(Var1 == snp| Var2 == snp)
  if(dim(cor)[1] > 0){
    rm_snp <- c(rm_snp,snp)
    cor_data <- cor_data %>%
      filter(!Var1 == snp & !Var2 == snp)
  }
}


fdr_cor_G <- -log10(p.adjust(res[[3]],method = "fdr"))

df_logP_G_adj <- matrix(fdr_cor_G, nrow = 42,ncol = 42) # , names(genetics)[2:length(genetics)])
rownames(df_logP_G_adj) <- names(genetics)[2:length(genetics)]
colnames(df_logP_G_adj) <- names(genetics)[2:length(genetics)]
longData <- melt(df_logP_G_adj)
longData<-longData[longData$value!=0,]
ggplot(longData, aes(x = Var2, y = Var1)) +
  geom_raster(aes(fill=value)) +
  scale_fill_gradient(low="grey90", high="red") +
  labs(x="SNPs", y="SNPs", title="-Log10P") +
  theme_bw() + theme(axis.text.x=element_text(size=9, angle=0, vjust=0.3),
                     axis.text.y=element_text(size=9),
                     plot.title=element_text(size=11))

SNP_info <- read.csv(paste0("Simple_gwas_25_codes/",ds_id,"_SNP_OR.csv")) 
snp_lst <- names(genetics)[2:length(genetics)] %>%
  strsplit( "_") %>% 
  sapply(function(x) x[1])
pos1<- c()
pos2<- c()
dst_BP <- c()
corr_sgf <- c() # save the significance level
for(i in 1:length(snp_lst)){
  if(length(which(df_logP_G_adj[,i] > 1))){
    BP1 <- SNP_info %>% 
      filter(SNP == snp_lst[i]) %>%
      pull(3)
    BP2 <- SNP_info %>%
      filter(SNP %in% snp_lst[which(df_logP_G_adj[,i] > 1)]) %>%
      pull(3)
    pos1<- c(pos1, rep(snp_lst[i], length(which(df_logP_G_adj[,i] > 1)))) 
    pos2<- c(pos2, snp_lst[which(df_logP_G_adj[,i] > 1)])  
    dst_BP <- c(dst_BP, abs(BP1 - BP2))
  }
}
SNP_dst <- data.frame(pos1, pos2, dst_BP)
longData <- melt(SNP_dst)
longData<-longData[longData$value!=0,]
ggplot(longData, aes(x = pos1, y = pos2)) +
  geom_raster(aes(fill=log10(value))) +
  scale_fill_gradient(low="grey90", high="red") +
  labs(x="SNPs", y="SNPs", title="-Log10 of BP distance, only show those which are correlated") +
  theme_bw() + theme(axis.text.x=element_text(size=9, angle=0, vjust=0.3),
                     axis.text.y=element_text(size=9),
                     plot.title=element_text(size=11))

########################################################
# compare univariate/multivariate SNP profile 2020-05-20
########################################################
ds_id <- "C443"
useful_snp_data <- read.csv(paste0("univariate_estimation_rm_LD/", ds_id, "_snp_for_clustering.csv")) %>%
  mutate(OR = log(OR)) %>%
  mutate(OR = ifelse(SNP %in% flip_snp$SNP, -OR, OR))

snp_list <- useful_snp_data %>% 
  group_by(SNP) %>%
  slice(1) %>%
  pull(1)

plot_single_snp <- function(ds_id, snp){
  data_univariate <- read.csv(paste0("univariate_estimation_rm_LD/", ds_id, "_snp_for_clustering.csv")) %>%
    mutate(OR = log(OR)) %>%
    mutate(OR = ifelse(SNP %in% flip_snp$SNP, -OR, OR)) %>%
    dplyr::filter(SNP == snp)
  data_multivariate <- read.csv(paste0("join_estimation_rm_LD/", ds_id, "_multivariate_snp_for_clustering.csv")) %>%
    mutate(OR = log(OR)) %>%
    mutate(OR = ifelse(SNP %in% flip_snp$SNP, -OR, OR)) %>%
    dplyr::filter(SNP == snp)

  plt <- ggplot(data_univariate, aes(PHE)) + 
    geom_line(aes(y = OR), color = red) +
    geom_ribbon(aes(ymin=OR-SE, ymax=OR+SE), fill= red, alpha=.3)+ 
    geom_line(data = data_multivariate, aes(y = OR),color = blue) +
    geom_ribbon(data = data_multivariate,aes(ymin=OR-SE, ymax=OR+SE), fill= blue,alpha=.3)+
    ggtitle(snp) + 
    theme(legend.position = "none",panel.background=element_blank()) + 
    xlab("Age group") + ylab("Log Odds Ratio")
  ggsave(paste0("~/Desktop/Writting_up_genetic_risk/supplementary/LD_single_snp_compare_",ds_id,"_",snp, ".png"), plt, width = 5, height = 5)
}
sapply(snp_list, function(x) plot_single_snp(ds_id, x))

#############################################
# small task, do simulation ploting online
# 2020-05-26
#############################################
# save the many objects needed
s_lst <- (-100:100)/10000
p_non_flat <- function(SIM){
  rslt_alt_hyp <- SIM[[3]]
  rslt_null_hyp <- SIM[[4]]
  ratio_ll <- 2*(rslt_alt_hyp[[6]] - rslt_null_hyp[[6]])
  return(1 - pchisq(ratio_ll, df=2))
}
# save coverage
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

print("Non-flat")
p_covrg <- matrix(nrow = length(s_lst), ncol = 8)
p_non_flat_test <- matrix(nrow = length(s_lst), ncol = 400)
for(i in 1:201){
  print(i)
  load(paste0("Non_flat_test_jointly/joint_non_flat_test",i,".Rdata"))
  p_non_flat_test[i,1:length(SIM)] <- sapply(1:length(SIM), function(x) ifelse(length(SIM[[x]]) == 4,p_non_flat(SIM[[x]]),NA))
  if(is.na(SIM[[1]][[1]])){
    cv <- SIM[[2]][[1]]$latent_curve
  }else cv <- SIM[[1]][[1]]$latent_curve
  p_covrg[i,] <- rowMeans(sapply(1:length(SIM), function(x) comp_covrg(SIM[[x]], cv) ), na.rm = TRUE)
}
p_thre <- 0.05
power_mean <- rowMeans(p_non_flat_test < p_thre, na.rm = T)
power_var <- sqrt(power_mean*(1-power_mean)/400)
df_non_flat_test <- data.frame("Slope" = s_lst, "Power" = power_mean, "Power_std" = power_var)

power_mean <- rowMeans(p_covrg, na.rm = T)
power_std <- sqrt(power_mean*(1-power_mean)/3200)
df_covrg <- data.frame("Slope" = s_lst, "Power" = power_mean, "Power_std" = power_std)

print("Multi-curve")
# multi curve
s_lst <- (-150:150)/4000
p_multi_cluster <- function(SIM){
  rslt_alt_hyp <- SIM[[3]]
  rslt_null_hyp <- SIM[[4]]
  ratio_ll <- 2*(rslt_alt_hyp[[6]] - rslt_null_hyp[[6]])
  # degree of freedom is 5
  return(1 - pchisq(ratio_ll, df=4)) 
}
p_multi_curve_test <- matrix(nrow = length(s_lst), ncol = 400)
for(i in 1:301){ # for(i in 1:201)
  load(paste0("multi_curve_test_jointly/multi_curve_test_jointly_",i,".Rdata"))
  # p_multi_curve_test[i,] <- sapply(1:length(SIM), function(x) ifelse(anyNA(SIM[[x]]),NA,p_multi_cluster(SIM[[x]])))
  p_multi_curve_test[i,] <- sapply(1:length(SIM), function(x) ifelse(length(SIM[[x]]) == 4,p_multi_cluster(SIM[[x]]),NA))
}
p_thre <- 0.05
power_mean <- rowMeans(p_multi_curve_test< p_thre, na.rm = TRUE)
power_var <- sqrt(power_mean*(1-power_mean)/400)
df_multi_curv_test <- data.frame("Slope" = s_lst, "Power" = power_mean, "Power_std" = power_var)

print("Frailty")
# frailty
u_lst <- exp(0.1*(1:100))
p_frailty <- matrix(nrow = 100, ncol = 400)
for(i in c(1:100)){
  load(paste0("frailty_test_jointly/jointly_frailty_test",i,".Rdata"))
  p_frailty[i,1:length(SIM)] <- sapply(1:length(SIM), function(x) ifelse(length(SIM[[x]]) == 4,p_non_flat(SIM[[x]]),NA))
}
p_thre <- 0.05
power_mean <- rowMeans(p_frailty < p_thre, na.rm = TRUE)
power_var <- sqrt(power_mean*(1-power_mean)/400)
df_frailty <- data.frame("frailty" = 1/u_lst, "Power" = power_mean, "Power_std" = power_var)

print("Old population selection")
# healthy population
s_lst <- (-100:100)/100
p_non_flat_test <- matrix(nrow = length(s_lst), ncol = 400)
for(i in 1:201){
  load(paste0("healthy_elder_effect/joint_healthy_elder_test",i,".Rdata"))
  p_non_flat_test[i,1:length(SIM)] <- sapply(1:length(SIM), function(x) ifelse(length(SIM[[x]]) == 4,p_non_flat(SIM[[x]]),NA))
}
p_thre <- 0.05
power_mean <- rowMeans(p_non_flat_test < p_thre, na.rm = T)
power_var <- sqrt(power_mean*(1-power_mean)/400)
df_age_effect <- data.frame("Slope" = s_lst, "Power" = power_mean, "Power_std" = power_var)

save(df_non_flat_test, df_covrg, df_multi_curv_test, df_frailty, df_age_effect, file = "simulation_power.Rdata")


load("simulation_power.Rdata")
plt <- ggplot(df_non_flat_test) + 
  geom_line(aes(x=Slope, y = Power), color = blue) +
  geom_ribbon(aes(x=Slope, ymax = Power + 1.96*Power_std, ymin = Power - 1.96*Power_std ), alpha = 0.3) +
  # geom_smooth(aes(x=Slope, y = P), span = 0.1) + 
  theme(legend.position = "none",panel.background=element_blank()) +
  geom_hline(yintercept= 0.05, linetype="dashed", color = red,size = 1) + 
  ylab("Power")

ggsave("~/Desktop/Writting_up_genetic_risk/Joint_estimation_non_flat.png", width = 6, height = 4)

plt <- ggplot(df_multi_curv_test) + 
  geom_line(aes(x=Slope, y = Power), color = red) +
  geom_ribbon(aes(x=Slope, ymax = Power + 1.96*Power_std, ymin = Power - 1.96*Power_std ), alpha = 0.3) +
  # geom_smooth(aes(x=Slope, y = P), span = 0.1) + 
  theme(legend.position = "none",panel.background=element_blank()) +
  geom_hline(yintercept= 0.05, linetype="dashed", color = red, size = 1) + 
  ylab("Power")
ggsave("~/Desktop/Writting_up_genetic_risk/Joint_estimation_multiple_curve.png", width = 6, height = 4)

plt <- ggplot(df_covrg) + 
  geom_line(aes(x=Slope, y = Power), color = blue) +
  geom_ribbon(aes(x=Slope, ymax = Power + 1.96*Power_std, ymin = Power - 1.96*Power_std ), alpha = 0.3) +
  # geom_smooth(aes(x=Slope, y = P), span = 0.1) + 
  geom_hline(yintercept= 0.95, linetype="dashed", color = red, size = 1) + 
  theme(legend.position = "none",panel.background=element_blank()) +
  ylab("Coverage")

ggsave("~/Desktop/Writting_up_genetic_risk/coverage_over_slope.png" ,plt, width = 6, height = 4)

plt <- ggplot(df_frailty) + 
  geom_line(aes(x=frailty, y = Power), color = blue) +
  geom_ribbon(aes(x=frailty, ymax = Power + 1.96*Power_std, ymin = Power - 1.96*Power_std ), alpha = 0.3) +
  # geom_smooth(aes(x=Slope, y = P), span = 0.1) + 
  theme(legend.position = "none",panel.background=element_blank()) +
  ylab("False positive rate") + 
  scale_x_continuous(trans='log10') +
  geom_hline(yintercept= 0.05, linetype="dashed", color = red,size = 1) + 
  xlab("Gamma scale parameter (variance of risk)")
ggsave("~/Desktop/Writting_up_genetic_risk/frailty_power.png" ,plt, width = 6, height = 4)


plt <- ggplot(df_age_effect) + 
  geom_line(aes(x=Slope, y = Power), color = blue) +
  geom_ribbon(aes(x=Slope, ymax = Power + 1.96*Power_std, ymin = Power - 1.96*Power_std ), alpha = 0.3) +
  # geom_smooth(aes(x=Slope, y = P), span = 0.1) + 
  theme(legend.position = "none",panel.background=element_blank()) +
  ylab("False positive rate")

ggsave("~/Desktop/Writting_up_genetic_risk/Age_effect.png" ,plt, width = 6, height = 4)


load("simulation_linear.Rdata")
plt <- ggplot(df_non_flat_test) + 
  geom_line(aes(x=Slope, y = Power), color = blue) +
  geom_ribbon(aes(x=Slope, ymax = Power + 1.96*Power_std, ymin = Power - 1.96*Power_std ), alpha = 0.3) +
  # geom_smooth(aes(x=Slope, y = P), span = 0.1) + 
  theme(legend.position = "none",panel.background=element_blank()) +
  geom_hline(yintercept= 0.05, linetype="dashed", color = red,size = 1) + 
  ylab("Power")

ggsave("~/Desktop/Writting_up_genetic_risk/linear_estimation_non_flat.png", width = 6, height = 4)

plt <- ggplot(df_multi_curv_test) + 
  geom_line(aes(x=Slope, y = Power), color = red) +
  geom_ribbon(aes(x=Slope, ymax = Power + 1.96*Power_std, ymin = Power - 1.96*Power_std ), alpha = 0.3) +
  # geom_smooth(aes(x=Slope, y = P), span = 0.1) + 
  theme(legend.position = "none",panel.background=element_blank()) +
  geom_hline(yintercept= 0.05, linetype="dashed", color = red, size = 1) + 
  ylab("Power")
ggsave("~/Desktop/Writting_up_genetic_risk/linear_estimation_multiple_curve.png", width = 6, height = 4)

# plot linear together with quadratic polynomial
load("simulation_linear.Rdata")
df_tmp <- df_non_flat_test
load("simulation_power.Rdata")
plt <- ggplot(df_non_flat_test) + 
  geom_line(aes(x=Slope, y = Power, color = "Quadratic polynomial")) +
  geom_ribbon(aes(x=Slope, ymax = Power + 1.96*Power_std, ymin = Power - 1.96*Power_std ), alpha = 0.3) +
  geom_line(data = df_tmp,aes(x=Slope, y = Power, color = "Linear model")) +
  geom_ribbon(data = df_tmp,aes(x=Slope, ymax = Power + 1.96*Power_std, ymin = Power - 1.96*Power_std ), alpha = 0.3, fill = blue) +
  # geom_smooth(aes(x=Slope, y = P), span = 0.1) + 
  theme(legend.position = "none",panel.background=element_blank()) +
  geom_hline(yintercept= 0.05, linetype="dashed", color = red,size = 1) + 
  ylab("Power") + 
  scale_color_manual(name="Model type",values=c("Linear model" = blue, "Quadratic polynomial" = red))

ggsave("~/Desktop/Writting_up_genetic_risk/quadratic_linear_estimation_non_flat.png", width = 6, height = 4)

plt <- ggplot(df_non_flat_test) + 
  geom_line(aes(x=Slope, y = Power, color = "Quadratic polynomial")) +
  geom_ribbon(aes(x=Slope, ymax = Power + 1.96*Power_std, ymin = Power - 1.96*Power_std ), alpha = 0.3) +
  geom_line(data = df_tmp,aes(x=Slope, y = Power, color = "Linear model")) +
  geom_ribbon(data = df_tmp,aes(x=Slope, ymax = Power + 1.96*Power_std, ymin = Power - 1.96*Power_std ), alpha = 0.3, fill = blue) +
  # geom_smooth(aes(x=Slope, y = P), span = 0.1) + 
  theme(legend.position = "left",panel.background=element_blank()) +
  geom_hline(yintercept= 0.05, linetype="dashed", color = red,size = 1) + 
  ylab("Power") + 
  scale_color_manual(name="Model type",values=c("Linear model" = blue, "Quadratic polynomial" = red))

ggsave("~/Desktop/Writting_up_genetic_risk/legend_estimation_non_flat.png", width = 6, height = 4)

load("simulation_linear.Rdata")
df_tmp <- df_multi_curv_test
load("simulation_power.Rdata")
plt <- ggplot(df_multi_curv_test) + 
  geom_line(aes(x=Slope, y = Power, color = "Quadratic polynomial")) +
  geom_ribbon(aes(x=Slope, ymax = Power + 1.96*Power_std, ymin = Power - 1.96*Power_std ), alpha = 0.3) +
  geom_line(data = df_tmp,aes(x=Slope, y = Power, color = "Linear model")) +
  geom_ribbon(data = df_tmp,aes(x=Slope, ymax = Power + 1.96*Power_std, ymin = Power - 1.96*Power_std ), alpha = 0.3, fill = blue) +
  # geom_smooth(aes(x=Slope, y = P), span = 0.1) + 
  theme(legend.position = "none",panel.background=element_blank()) +
  geom_hline(yintercept= 0.05, linetype="dashed", color = red,size = 1) + 
  ylab("Power") + 
  scale_color_manual(name="Model type",values=c("Linear model" = blue, "Quadratic polynomial" = red))

ggsave("~/Desktop/Writting_up_genetic_risk/quadratic_linear_estimation_multiple_curve.png", width = 6, height = 4)

#############################################################
# checking: compute log likelihood using BETA and check prior 
#############################################################
source("Genetic_longitudinal_functions.R")
ds_id <- "C443"

############################################################################# 
# following part could be put into univariate_estimation_survival.R directly
#############################################################################

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


BETA <- list_of_rslt[[1]][[4]]
pi_saver <- list_of_rslt[[1]][[1]]
# load the multivariate fitting data
para <- load_SNP_coef(SNP_FILE = paste("join_estimation_rm_LD/",ds_id,"_multivariate_snp_for_clustering.csv", sep = ""),
                      FLIP_FILE = paste("SNP_info/",ds_id,"_flipped_snp.csv", sep = ""))

d_beta_hat <- list()
# BETA is from univariate fitting, while the data is from joint
d_beta_hat[[1]] <- sapply(1:para$S, function(x) 
  dmvnorm(t(para$betaj[[x]]), mean = BETA[[1]], sigma = para$sigmaSAMPLE[[x]]))
print(comp_ll(pi_saver,d_beta_hat, para))

############################################################################# 
# Do a sign test of the decreasing pattern and frailty 
#############################################################################
library(BSDA)
# compare how the decreasing pattern is explained by the unobserved coefficient model
pt <- "*_snp_for_clustering.csv"
temp <- list.files(paste("univariate_estimation_rm_LD/", sep=""), pattern=pt)

ds_list <- lapply(temp, function(x) strsplit(x, "_")[[1]][1]) 
ds_list <- data.frame(coding = unlist(ds_list)) %>%
  filter(! coding == "I258", ! coding == "Z951", ! coding == "Z955")
t_st_lst <- rep(NA, length(ds_list[[1]]))
total_slop <- list()
for(id in 1:length(ds_list[[1]])){
  ds_id <- as.character(ds_list[[1]][id])
  load(paste("./univariate_estimation_rm_LD/",ds_id,"_true_data.RData", sep = ""))
  load("adjustment.Rdata")
  adj <- adj[[ds_id]]
  # determine which profile to show
  
  rslt <- list_of_rslt[[1]]
  
  BETA <- rslt[[4]]

  # # get the momentum
  mom_fit <- rep(0,length(BETA[[1]]) - 1)
  mom_adj <- rep(0,length(BETA[[1]]) - 1)
  for(i in 1:(length(BETA[[1]]) - 1)){
    mom_fit[i] <- (BETA[[1]][i+1] - BETA[[1]][i])/BETA[[1]][i]
    mom_adj[i] <- (adj[i+1] - adj[i])/adj[i]
  }
  df_mom <- data.frame("PHE"= (1:7)+.5)
  df_mom["fitting"] <- mom_fit
  df_mom["adjustment"] <- mom_adj
  
  # test_st <- SIGN.test(mom_adj - mom_fit)
  test_st <- t.test(mom_adj,mom_fit,paired=TRUE, alternative="greater")
  print(c(ds_id,test_st$p.value < 0.05))
  t_st_lst[id] <- test_st$p.value
  total_slop[[ds_id]] <- c(sum(mom_fit), sum(mom_adj))
}
# read the q-value for comparison
tb_all <- read.csv("~/Desktop/Writting_up_genetic_risk/univariate_q_value.csv") %>%
  mutate(log10q1 = -log10(q1))
tb_all$t_test <- -log10(p.adjust(t_st_lst,method = "fdr") )
# plot logP
plt <-ggplot(tb_all, aes(x = coding)) + 
  geom_point(aes(y = log10q1, color = "Permutation test of non-constant effects", shape = "Permutation test of non-constant effects"), size = 3) + 
  geom_point(aes(y = t_test, color = "Paired t-test of frailty and profile slope",shape = "Paired t-test of frailty and profile slope"), size = 3) +
  theme(legend.position="bottom",panel.background=element_blank()) + 
  ylab("-Log10Q (FDR adjusted)") +
  # scale_y_continuous(trans='log10') + # using log scale
  xlab("Disease (ICD 10 code)") +
  geom_hline(yintercept=-log10(0.1), linetype="dashed", color = red) +
  scale_color_manual(name="Test",values=c("Paired t-test of frailty and profile slope" = green, "Permutation test of non-constant effects" = blue))+
  scale_shape_manual(name = "Test", values = c("Permutation test of non-constant effects"= 16, "Paired t-test of frailty and profile slope"= 2))

ggsave("~/Desktop/Writting_up_genetic_risk/Paired-t-test.png",plt,width = 8, height =5) 
##########################################
# fit a linear model to the data 2020-5-28
##########################################
source("Genetic_longitudinal_functions.R")
pt <- "*_snp_for_clustering.csv"
temp <- list.files(paste("join_estimation_rm_LD/", sep=""), pattern=pt)

ds_list <- lapply(temp, function(x) strsplit(x, "_")[[1]][1]) 
ds_list <- data.frame(coding = unlist(ds_list)) %>%
  filter(! coding == "I258", ! coding == "Z951", ! coding == "Z955")
intercept_profile <- rep(NA, length(ds_list[[1]]))
slope_profile <- rep(NA, length(ds_list[[1]]))
g_risk_ratio_late_age <- rep(NA, length(ds_list[[1]]))
for(id in 1:length(ds_list[[1]])){
  ds_id <- as.character(ds_list[[1]][id])
  snp_file <- paste("./join_estimation_rm_LD/",ds_id,"_multivariate_snp_for_clustering.csv", sep = "")
  flip_file <- paste("SNP_info/",ds_id,"_flipped_snp.csv", sep = "")
  para <- load_SNP_coef(SNP_FILE = snp_file, FLIP_FILE = flip_file)
  para$p <- 2
  para$X <- X_cb_spline(para$p,8)
  para$sigma0inv <- solve(diag(para$p) * 1) 
  rslt <- EM_K(1,100, para)
  theta <- rslt[[3]]
  intercept_profile[id] <- theta[[1]][1]
  slope_profile[id] <- theta[[1]][2]
  g_risk_ratio_late_age[id] <- (theta[[1]][1] + theta[[1]][2]*8)/(theta[[1]][1] + theta[[1]][2]*1)
}
df_linear_fit <- data.frame(ds_list, intercept_profile, slope_profile, g_risk_ratio_late_age)
write.csv(df_linear_fit, file = "~/Desktop/Writting_up_genetic_risk/linear_fitting.csv")


###################################################
# compute risk for individuals with I251 2020-06-12
###################################################
new_data <- read.csv("all_data.csv") %>%
  mutate(diag_icd10 = as.character(diag_icd10)) %>%
  mutate(eid = as.character(eid))

# delete withdraw group
with_draw <- read.csv(file="~/Desktop/genetics_longitudinal_data/longitudinal_data/with_draw_2020Feb11.csv", 
                      header=FALSE) 
keep <- read.table(file="~/Desktop/genetics_longitudinal_data/longitudinal_data/keep.txt", 
                   header=FALSE, sep=" ") %>% 
  anti_join(with_draw, by="V1")
keep <- keep %>%
  mutate(V1 = as.character(V1))
survive_age <- read.csv("survive_age.csv") %>%
  mutate(eid = as.character(eid))

pt <- "*_snp_for_clustering.csv"
temp <- list.files(paste("join_estimation_rm_LD/", sep=""), pattern=pt)

ds_list <- lapply(temp, function(x) strsplit(x, "_")[[1]][1]) 
ds_list <- data.frame(coding = unlist(ds_list)) %>%
  filter(! coding == "I258", ! coding == "Z951", ! coding == "Z955")

early_baseline_hz <- rep(NA, length(ds_list[[1]]))
late_baseline_hz <- rep(NA, length(ds_list[[1]]))
genetic_risk_early <- rep(NA, length(ds_list[[1]]))
genetic_risk_late<- rep(NA, length(ds_list[[1]]))
absolute_early_hz <- rep(NA, length(ds_list[[1]]))
absolute_late_hz<- rep(NA, length(ds_list[[1]]))
genetic_hz_ratio<-absolute_late_hz<- rep(NA, length(ds_list[[1]]))
abs_hz_ratio<-absolute_late_hz<- rep(NA, length(ds_list[[1]]))
for(id in 1:length(ds_list[[1]])){
  ds_id <- as.character(ds_list[[1]][id])
  print(ds_id)
  case_data <- read.table(paste0("Survival/",ds_id,"_pheno_over_age.txt"), header=F) %>%
    rename(eid = V1, time = V2, group = V3, censor = V4) %>% 
    filter(censor == 1) 
  
  earlySurvived <- survive_age %>%
    filter(survive_year >= 45) %>% 
    summarise(n()) %>% pull
  
  ds_earlier_than70 <- case_data %>%
    mutate(eid = as.character(eid)) %>%
    filter(group < 7) 
  lateSurvived <- survive_age %>%
    filter(survive_year >= 70) %>% 
    anti_join(ds_earlier_than70, by = "eid") %>%
    summarise(n()) %>% pull
  
  early_rate <- case_data %>%
    filter(group == 2) %>%
    summarise(n()) %>% pull
  early_rate <- early_rate/earlySurvived
  
  late_rate <- case_data %>%
    filter(group == 7) %>%
    summarise(n()) %>% pull
  late_rate <- late_rate/lateSurvived
  
  early_baseline_hz[id] <- early_rate
  late_baseline_hz[id] <- late_rate
  #################################
  # compute the genetic risk
  #################################
  patients <- case_data %>%
    select(eid)
  
  load("SNP_data_LD_rm_20200521.RData")
  pos <- which(t$coding == ds_id) 
  SNP_lst <- tables[[pos]]$SNP %>%
    make.names()
  
  genetics <- read.table(paste0("Survival/",ds_id,".raw"), header=T) %>%
    select(-IID, -PAT, -MAT, -SEX, -PHENOTYPE) %>%
    semi_join(patients,by = c("FID" = "eid")) 
  genetics <- setNames(genetics, sapply(names(genetics), function(x) strsplit(x, "_")[[1]][1]))
  genetics[is.na(genetics)] <- 0
  
  flip_snp <- read.csv(paste("SNP_info/",ds_id,"_flipped_snp.csv", sep = ""))
  
  snp_for_clustering <- read.csv(paste0("join_estimation_rm_LD/", ds_id, "_multivariate_snp_for_clustering.csv"))
  
  useful_snp_data <- snp_for_clustering %>% 
    # I need to remove the NAs
    mutate(SE = ifelse(is.na(SE), 10, SE)) %>%
    mutate(OR = ifelse(is.na(OR), 1, OR)) %>%
    mutate(OR = log(OR)) %>%
    # I will have to control the overflowing of coef (some -> inf)
    mutate(SE = ifelse(abs(OR) > 10, 10, SE)) %>%
    mutate(OR = ifelse(abs(OR) > 10, 0, OR)) # %>%
  # mutate(OR = ifelse(SNP %in% flip_snp$SNP, -OR, OR))
  
  early_effect <- useful_snp_data %>%
    filter(PHE == 2)
  
  late_effect <- useful_snp_data %>%
    filter(PHE == 7)
  
  early_risk <- rep(0, dim(genetics)[1])
  late_risk <- rep(0, dim(genetics)[1])
  
  for(i in 1:length(SNP_lst)){
    alleles <- genetics %>% select(SNP_lst[i]) %>% pull
    # if(SNP_lst[i] %in% flip_snp$SNP){
    #   alleles <- 2 - alleles
    # }
    early_sz <- early_effect %>% filter(SNP == SNP_lst[i]) %>% pull(2)
    early_risk <- early_risk + alleles * early_sz  
    late_sz <- late_effect %>% filter(SNP == SNP_lst[i]) %>% pull(2)
    late_risk <- late_risk + alleles * late_sz 
  }
  idx <- order(early_risk,decreasing = T)[1:floor(length(early_risk)/10)]
  genetic_hz_ratio[id] <- exp(mean(early_risk[idx]) - mean(late_risk[idx]))
  genetic_risk_early[id] <- exp(mean(early_risk[idx]))
  genetic_risk_late[id] <- exp(mean(late_risk[idx]))
  absolute_early_hz[id] <- exp(mean(early_risk[idx])) * early_baseline_hz[id]
  absolute_late_hz[id] <- exp(mean(late_risk[idx])) * late_baseline_hz[id]
  abs_hz_ratio[id] <- absolute_early_hz[id]/absolute_late_hz[id]
}
risk_early_late <- data.frame(ds_list,early_baseline_hz, late_baseline_hz,genetic_risk_early, genetic_risk_late, absolute_early_hz,absolute_late_hz, genetic_hz_ratio, abs_hz_ratio)
write.csv(risk_early_late, file = "~/Desktop/Writting_up_genetic_risk/early_late_risk_compare_first_decile.csv")

#############################################################
# Linear plotting 2: get the P-value plot!
#############################################################
setwd("/Users/xilin/Desktop/genetics_longitudinal_data/longitudinal_data/")
# load the disease list
pt <- "*_snp_for_clustering.csv"
temp <- list.files(paste("linear_joint_estimation_rm_LD//", sep=""), pattern=pt)
icd_code <- lapply(temp, function(x) strsplit(x, "_")[[1]][1])
icd_list <- data.frame(coding = unlist(icd_code)) %>%
  filter(! coding == "I258", ! coding == "Z955", ! coding == "Z951")

P_all_ds <- matrix(0,dim(icd_list)[1], 4) # the test is, 1) against flat 2) two against 1, 3) three agains two 4) four against three 5) against frailty
cnt <- matrix(nrow = dim(icd_list)[1])
for(id in 1:dim(icd_list)[1]){
  try({
  ds_id <- icd_list$coding[id]
  temp <- list.files("./linear_joint_permutation_rm_LD/", pattern=as.character(ds_id))
  # need to deal with potential NAs 
  permute_ll <- list()
  cnt[id] <- length(temp) * 100
  for(i in 1:length(temp)){
    ll <- read.csv(paste("./linear_joint_permutation_rm_LD/", temp[i], sep = ""))
    permute_ll[[i]] <- ll
  }
  permute_ll <- bind_rows(permute_ll)
  if(anyNA(permute_ll)){
    cnt[id] <- (cnt[id] - sum(is.na(permute_ll))/6)
    permute_ll[is.na(permute_ll)] <- -Inf
  }
  # print for debugging
  print(c(ds_id, cnt[id]))
  per_ll_gain <- permute_ll[,c(1,2,3,4,1)] - permute_ll[,c(5,1,2,3,6)]
  # remove the NA value from Inf - Inf
  per_ll_gain[is.na(per_ll_gain)] <- -Inf
  # load the likelihood
  true_likelihood <- read.csv(paste("./linear_joint_estimation_rm_LD/",ds_id,"_ll_matrix.csv", sep = ""))
  # first get the p for using two flexible line again one flexible line
  true_ll_gain <- true_likelihood$x[c(1,2,3,4)] - true_likelihood$x[c(5,1,2,3)]
  p_value <- sapply(1:length(true_ll_gain), function(s) sum(per_ll_gain[,s] > true_ll_gain[s])/cnt[id])
  P_all_ds[id,] <- p_value
  })
}

fdr1 <- p.adjust(P_all_ds[,1],method = "fdr")  + 1/cnt
fdr2 <- p.adjust(P_all_ds[,2],method = "fdr")  + 1/cnt
fdr3 <- p.adjust(P_all_ds[,3],method = "fdr")  + 1/cnt
fdr4 <- p.adjust(P_all_ds[,4],method = "fdr")  + 1/cnt
# fdr5 <- p.adjust(P_all_ds[,5],method = "fdr") + 1/cnt
fdr_list <- icd_list
fdr_list$test_against_flat <- -log10(fdr1)
fdr_list$test_two_flexible_against_one_flexible <- -log10(fdr2)
fdr_list$test_three_against_two <- -log10(fdr3)
fdr_list$test_four_against_three <- -log10(fdr4) 
# fdr_list$test_against_frailty <- -log10(fdr5)

plt <- ggplot(fdr_list, aes(x = coding)) + 
  geom_point(size = 3, aes(y = test_against_flat,color = "Single Non-constant against constant", shape = "Single Non-constant against constant")) + 
  theme(panel.background=element_blank()) + 
  ylab("-Log10Q (FDR adjusted)") +
  # scale_y_continuous(trans='log10') + # using log scale
  xlab("Disease (ICD 10 code)") +
  geom_point(aes(x = coding, y = test_two_flexible_against_one_flexible, color = "Two clusters against one cluster", shape = "Two clusters against one cluster"), size = 3) +
  geom_hline(yintercept=-log10(0.1), linetype="dashed", color = red) +
  geom_point(aes(x = coding, y = test_three_against_two, color = "Three clusters against two clusters", shape = "Three clusters against two clusters"), size = 3) +
  geom_point(aes(x = coding, y = test_four_against_three, color = "Four clusters against three clusters", shape = "Four clusters against three clusters"), size = 3) +
  scale_color_manual(name="Test",values=c("Single Non-constant against constant" = blue, "Two clusters against one cluster" = red,"Three clusters against two clusters" = orange, "Four clusters against three clusters" = green))+
  scale_shape_manual(name = "Test", values = c("Single Non-constant against constant"= 19, "Two clusters against one cluster"= 18, "Three clusters against two clusters" = 19,"Four clusters against three clusters" = 18))
# geom_point(data = fdr_list, aes(x = coding, y = test_against_frailty), color = "#D55E00", shape = 8, size = 5) 
ggsave("~/Desktop/Writting_up_genetic_risk/Linear_fitting_FDR_joint.png",plt,width = 10, height =5) 

setwd("/Users/xilin/Desktop/genetics_longitudinal_data/longitudinal_data/")
# load the disease list
pt <- "*_snp_for_clustering.csv"
temp <- list.files(paste("Joint_estimation/", sep=""), pattern=pt)
icd_code <- lapply(temp, function(x) strsplit(x, "_")[[1]][1])
icd_list <- data.frame(coding = unlist(icd_code)) %>%
  filter(! coding == "I258", ! coding == "Z955", ! coding == "Z951")

q_list <- data.frame(icd_list, fdr1, fdr2, P_all_ds[,1]+ 1/cnt, P_all_ds[,2]+ 1/cnt) %>%
  rename(q1 = fdr1, q2 = fdr2, p1 = P_all_ds...1....1.cnt, p2 = P_all_ds...2....1.cnt) %>%
  select(coding, p1, q1, p2, q2)

write.csv(q_list, "~/Desktop/Writting_up_genetic_risk/Linear_fittingLD_indep_multivariate_q_value.csv")





#################################################################
# compute age-specific PRS power, where PRS are computed from a normal one 2020-08-18
#################################################################
require(survival)
require(dplyr)
library("mvtnorm")
set.seed(19940110)

data <- read.table("covariates_for_asso_over_time.txt", sep="\t") %>%
  select(-V2,  -V3) %>%
  rename(eid = V1)

prevelance <- read.table("HES_prevalence_by_code.txt", 
                         sep = "\t", quote = "", header  = TRUE)
snp_assoc <- read.table("SNP_Assocs_per_code.txt", 
                        sep = "\t", quote = "", header  = TRUE)
snp_assoc <- snp_assoc %>%
  right_join(prevelance, by = c("coding" = "Code"))
snp_of_interest <- snp_assoc %>% 
  filter(NUMSNPASSOCS > 20, Prevalence > .5) %>% # using prevalence > .5 will get 6 more interesting diseases
  select(coding, NUMSNPASSOCS, meaning)
one_cluster <- T # FALSE

setwd("/Users/xilin/Desktop/genetics_longitudinal_data/longitudinal_data/")
# load the disease list
pt <- "*_snp_for_clustering.csv"
temp <- list.files(paste("Joint_estimation/", sep=""), pattern=pt)
icd_code <- lapply(temp, function(x) strsplit(x, "_")[[1]][1])
icd_list <- data.frame(coding = unlist(icd_code)) %>%
  filter(! coding == "I258", ! coding == "Z955", ! coding == "Z951")
plt <-list()

for(id in 1:dim(icd_list)[1]){
  ds_id <- icd_list$coding[id]
  print(ds_id)
  genetics <- read.table(paste0("Survival/",ds_id,".raw"), header=T) %>%
    select(-IID, -PAT, -MAT, -SEX, -PHENOTYPE) 
  
  data_coxph <- read.table(paste0("Survival/",ds_id,"_pheno_over_age.txt"), header=F) %>%
    rename(eid = V1, time = V2, group = V3, censor = V4)
  
  data_coxph <- data_coxph %>%
    left_join(data, by = "eid") %>%
    left_join(genetics, by = c("eid" = "FID"))
  
  # simple imputation since all SNPs are LD independent
  for(col_id in 45:dim(data_coxph)[2]){
    data_coxph[is.na(data_coxph[,col_id]),col_id] <- mean(data_coxph[,col_id], na.rm = T)
  }
  
  num_rep <- 20
  precision_logit_top10 <- matrix(NA, num_rep, 8)
  precision_logit_top20 <- matrix(NA, num_rep, 8)
  for(rp in 1:num_rep){
    case_idx <- list()
    control_idx <- list()
    for(gp in 1:8){
      group_case_data <- data_coxph %>% 
        filter(group == gp, censor == 1) 
      
      group_control_data <- data_coxph %>% 
        filter(group == gp, censor == 0)
      
      case_idx[[gp]] <- sample(1:dim(group_case_data)[1])
      control_idx[[gp]] <-  sample(1:dim(group_control_data)[1])
    }
    # obtain bootstrap samples
    cross_validation_samples <- list()
    cross_validation_top20 <- list()
    for(bt in 1:5){
      print(paste0("Rep: ",rp, "; cross validation fold: ", bt))
      testing_case_data <- list() 
      testing_control_data <- list()
      for(gp in 1:8){
        group_case_data <- data_coxph %>% 
          filter(group == gp, censor == 1) 
        group_control_data <- data_coxph %>% 
          filter(group == gp, censor == 0)
        
        num_case <- length(case_idx[[gp]])
        if(bt == 5){
          bt_gp_idx_case <- case_idx[[gp]][(1 + (bt -1)* floor(num_case/5)):num_case]
        }else{
          bt_gp_idx_case <- case_idx[[gp]][(1 + (bt -1)* floor(num_case/5)):(bt * floor(num_case/5))]
        }
        num_control <- length(control_idx[[gp]])
        if(bt == 5){
          bt_gp_idx_control <- control_idx[[gp]][(1 + (bt -1)* floor(num_control/5)):num_control]
        }else{
          bt_gp_idx_control <- control_idx[[gp]][(1 + (bt -1)* floor(num_control/5)):(bt * floor(num_control/5))]
        }
        
        testing_case_data[[gp]] <- group_case_data %>%
          slice(bt_gp_idx_case)
        
        testing_case_data[[gp]] <- group_case_data %>%
          sample_n(floor(dim(group_case_data)[1]/5))
        
        testing_control_data[[gp]] <- group_control_data %>%
          slice(bt_gp_idx_control)
      }
      
      testing_case_data <- bind_rows(testing_case_data)
      testing_control_data <- bind_rows(testing_control_data)
      testing_data <- bind_rows(testing_case_data, testing_control_data)
      
      data_training <- data_coxph %>%
        anti_join(testing_data,by = c("eid", "time", "group", "censor"))
      
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
      # survival_rslt <- matrix(NA, 8, length(SNPs))
      # for(gp in 1:8){
      #   df <- data_training %>%
      #     filter(group == gp)
      #   mdl <- coxph(formula = as.formula(paste("Surv(time, censor) ~", SNP_list, " + ", PC_list)), data = df)
      #   survival_rslt[gp, ] <- mdl$coefficients[1:length(SNPs)]
      # }
      
      ##############################
      # perform logistic regression
      ##############################
      data_logit_case <- data_training %>% 
        filter(censor == 1)
      
      data_logit_control <- data_training %>% 
        filter(censor == 0) %>%
        group_by(eid) %>%
        slice(1) %>% 
        anti_join(data_logit_case, by = "eid")
      
      data_logit <- bind_rows(data_logit_case, data_logit_control)
      logit_mdl <- glm(as.formula(paste("censor ~", SNP_list, " + ", PC_list)), data = data_logit, family = binomial)
      logit_coef <- logit_mdl$coefficients[2:(1+length(SNPs))] 
      
      # comparison of risk for each age interval
      # accuracy_survival <- rep(NA, 8)
      accuracy_logit <- rep(NA, 8)
      accuracy_logit_top20 <- rep(NA, 8)
      for(gp in 1:8){
        test_set_interval <- testing_data %>% 
          filter(group == gp)
        test_genetics <- test_set_interval[,45:dim(test_set_interval)[2]]
        
        # survival_risk <- as.matrix(test_genetics) %*% as.matrix(survival_rslt[gp,])
        # idx_survival <- order(survival_risk,decreasing = T)[1:floor(length(survival_risk)/5)]
        # accuracy_survival[gp] <- mean(test_set_interval$censor[idx_survival])
        
        logit_risk <- as.matrix(test_genetics) %*% as.matrix(logit_coef)
        
        # compute the log odds ratio
        idx_logit <- order(logit_risk,decreasing = T)[1:floor(length(logit_risk)/10)]
        accuracy_logit[gp] <- mean(test_set_interval$censor[idx_logit])
        idx_logit_top20 <- order(logit_risk,decreasing = T)[1:floor(length(logit_risk)/5)]
        accuracy_logit_top20[gp] <- mean(test_set_interval$censor[idx_logit_top20])
        
        # num_case <- test_set_interval %>% filter(censor == 1) %>% summarise(n()) %>% pull
        # num_control <- test_set_interval %>% filter(censor == 0) %>% summarise(n()) %>% pull
        # print(num_case/num_control)
      }
      cross_validation_samples[[bt]] <- accuracy_logit
      cross_validation_top20[[bt]] <- accuracy_logit_top20
    }
    precision_logit_top10[rp,] <- 4 * (Reduce("+", cross_validation_samples)/5)/(1-(Reduce("+", cross_validation_samples)/5))
    precision_logit_top20[rp,]<- 4 * (Reduce("+", cross_validation_top20)/5)/(1-(Reduce("+", cross_validation_top20)/5))
  }
  
  # get the meaning of the code
  title_meaning <- snp_of_interest %>% filter(coding == ds_id) %>% pull(3) %>% as.character()
  df_mean <- data.frame("AGE"= as.numeric(c(45,50,55,60,65,70,75,80)), top10 = colMeans(precision_logit_top10), top20 = colMeans(precision_logit_top20), 
                        err_top10 = apply(precision_logit_top10, 2, function(x) sd(x)/sqrt(length(x))), err_top20 = apply(precision_logit_top20, 2, function(x) sd(x)/sqrt(length(x))))
  
  plt[[id]] <- ggplot(df_mean, aes(x = AGE)) +
    #scale_y_continuous(labels = scales::percent) +
    scale_x_continuous(breaks = c(45,50,55,60,65,70,75,80),labels=c("<45", "45-50","50-55","55-60","60-65","65-70","70-75",">75")) +
    ggtitle(paste0("Disease ID:", title_meaning))+
    theme(legend.position = "none",panel.background=element_blank(),plot.title = element_text(size = 10, face = "bold"))  + 
    labs(title = paste0("Disease ID:", title_meaning), x = "Age (years)", y = "Odds ratio") + 
    geom_line(aes(y = top10), linetype = "dashed", color = red) + 
    geom_pointrange(aes(y = top10, ymin = top10 - 1.96*err_top10, ymax = top10 + 1.96*err_top10), size = .5, color = red) + 
    geom_line(aes(y = top20), linetype = "dashed", color = blue) + 
    geom_pointrange(aes(y = top20, ymin = top20 - 1.96*err_top20, ymax = top20 + 1.96*err_top20), size = .5, color = blue) 
  
  ggsave(paste0("~/Desktop/Writting_up_genetic_risk/supplementary/pridiction_precision_", 
                as.character(ds_id), ".png"), width = 5, height = 4,plt[[id]])
}


