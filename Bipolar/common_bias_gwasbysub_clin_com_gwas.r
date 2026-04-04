library(dplyr)
library(readr)
library(tidyr)
library(stringr)
library(ggplot2)
library(microshades)
library(openxlsx)
library(GenomicSEM)
library(data.table)

setwd(dir = "/local1/home/pazweifel/src")
#load the required files
covstruct_gwasbysub <- dget("matrices/ambivalent_bias_factor_clin_com.R")
sumstats <- readRDS("bipolar_bias/clin_comm_combinedbias_sumstats.R")

#specify the lavaan model
model <- "
Liability_clin =~ NA*Clin
Liability_comm =~ NA*Comm
Bias =~ NA*Combined_Bias + Comm

Liability_clin ~~ 0*Bias
Liability_comm ~~ 0*Bias
Clin ~~ 0*Clin
Comm ~~ 0*Comm
Combined_Bias ~~ 0*Combined_Bias
Clin ~~ Combined_Bias
Combined_Bias ~~ 0*Comm
Clin ~~ 0*Comm

Bias ~ SNP
Liability_comm ~ SNP
Liability_clin ~ SNP
"
    


#create the path to the object you want the output to store
bias_factor_gwas_path <- "bipolar_bias/GWAS_combined_bias_comm_liability.R"
start_time <- Sys.time()

cat("Starting commonfactorGWAS run at", Sys.time(), "\n")

bias_factor_bip_gwas <- userGWAS(covstruc = covstruct_gwasbysub, 
                                              SNPs = sumstats,
                                              model = model,
                                              sub = c("Bias ~ SNP", "Liability_comm ~ SNP", "Liability_clin ~ SNP"),
                                              toler = 1e-70,
                                              parallel = TRUE, 
                                              cores = 30,
                                              Q_SNP = FALSE
                                     )
saveRDS(bias_factor_bip_gwas, file = bias_factor_gwas_path)


runtime_mins <- as.numeric(difftime(Sys.time(), start_time, units="mins"))
cat("GWAS finished in", round(runtime_mins, 2), "minutes.\n")
cat("Finished at", Sys.time(), "\n")