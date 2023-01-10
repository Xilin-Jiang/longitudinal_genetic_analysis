# Genetic analysis for longitudinal data 
I will keep updating this document to provide three key utilities in the paper [The impact of age on genetic risk for common diseases](https://journals.plos.org/plosgenetics/article?id=10.1371/journal.pgen.1009723). 

1. Estimating genetic effect using diagnosis with age information. This provides an unbiased estimation of single SNP effect sizes when time-to-event data (survival data) are provided; note the common case-control regression framework would introduce bias, which is explained in S1 Analytical Note of [Jiang et al. 2021 PLoS Genetics](https://journals.plos.org/plosgenetics/article?id=10.1371/journal.pgen.1009723).

2. Clustering of SNPs based on age profiles and testing heterogeneity of effect across age bins. We note in further studies that this test could be performed on any E variables. It is straightforward to divide the population into groups based on a E variable and produce summary statistics for each E-bin. 

3. Estimating population risk distribution. Previous study has proposed [using instantaneous risk across age in populations to estimate the distribution of frailty](https://link.springer.com/article/10.2307/2061224). Here we provide an implementation using survival data of a population to fit a frailty model.   

#### The description of code describes how to estimate genetic risk with the interval-censored proportional hazard model; and inferring latent curves  .

The file [Genetic_longitudinal_functions.R](https://github.com/Xilin-Jiang/longitudinal_genetic_analysis/blob/main/Genetic_longitudinal_functions.R "Genetic_longitudinal_functions.R") has all functions needed for simulating disease event over age, and performing proportional hazard model analysis over the simulated population. 

###### What you need to for the simulation

The file [simulation_for_paper.R](https://github.com/Xilin-Jiang/longitudinal_genetic_analysis/blob/main/simulation_for_paper.R "simulation_for_paper.R") contains the code for the simulation shown in the paper. Though we recommend saparating blocks of code for submission to clusters. 

* We recommend using a central computation facillity to perform the simulation instead of a personal laptop.

* To submit job, you need to create a separate simulation file from [simulation_for_paper.R](https://github.com/Xilin-Jiang/longitudinal_genetic_analysis/blob/main/simulation_for_paper.R "simulation_for_paper.R") , depending on the tasks. There are four simulation you could perform, which are commented with SIMULATION BLOCK1/2/3/4; the separate simulation file should be created by deleting all the figure creating code (all the code at the end of [simulation_for_paper.R](https://github.com/Xilin-Jiang/longitudinal_genetic_analysis/blob/main/simulation_for_paper.R "simulation_for_paper.R") files, should be trivial to identify with the comments), and keep one of the  SIMULATION BLOCK code.

* the paramenter `s_id` at the beginning of the file serve as a handle for submission jobs with different parametrisation. For SIMULATION BLOCK1 `s_id` should be `1-201`, for SIMULATION BLOCK2, `s_id` should be `1-301`; for SIMULATION BLOCK3 and SIMULATION BLOCK4, `s_id`should be`1-100`.

* Once you have finished the simulations, the power test results could be generated using the plotting code.

###### Praparation for running with data

* Unfortunately we could not provide examples with UK Biobank dataset due to confidentiality; However, to adapt it to a dataset that you are familiar with should be easy. Experienced computational statistician could directly look into [Genetic_longitudinal_functions.R](https://github.com/Xilin-Jiang/longitudinal_genetic_analysis/blob/main/Genetic_longitudinal_functions.R "Genetic_longitudinal_functions.R") as well as [Supplementary Information](https://www.biorxiv.org/content/10.1101/2020.07.17.208280v1.supplementary-material) and carry out your own analysis. 

* To perform analysis with our code, you need to get four data files ready: one covariate files containing the covariate you want to include for survival analysis; one genetics file containing the genotypes; one case-control file that contains case control data for each group; one file containing the loci from the genetic file that have protective MAF. All files are discribed in  [joint_estimation_survival.R](https://github.com/Xilin-Jiang/longitudinal_genetic_analysis/blob/main/joint_estimation_survival.R "joint_estimation_survival.R") . 

[local_testing_file.R](https://github.com/Xilin-Jiang/longitudinal_genetic_analysis/blob/main/local_testing_file.R "local_testing_file.R") contains unorganised code to generate most of the figures, which is left in the repository for those interested. 
