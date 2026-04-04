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
covstruct_gwasbysub <- dget("matrices/neu_ext_comm_liability_factor.R")
sumstats <- readRDS("bipolar_bias/clin_comm_neu_sumstats.R")

#specify the lavaan model
model <- "
Liability_clin =~ NA*Clin
Liability_comm =~ NA*Comm
Bias =~ NA*Neu + Comm

Liability_clin ~~ 0*Bias
Liability_comm ~~ 0*Bias
Clin ~~ 0*Clin
Comm ~~ 0*Comm
Neu ~~ 0*Neu
Clin ~~ Neu
Neu ~~ 0*Comm
Clin ~~ 0*Comm

Bias ~ SNP
Liability_comm ~ SNP
Liability_clin ~ SNP
"
    


#create the path to the object you want the output to store
bias_factor_gwas_path <- "bipolar_bias/GWAS_Comm_liability.R"
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
cat("Model finished in", round(runtime_mins, 2), "minutes.\n")
cat("Finished at", Sys.time(), "\n")