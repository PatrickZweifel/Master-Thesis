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
covstruct_ambivalent_correlations <- dget("matrices/Clinical_ambivalent_external.R")
sumstats_ambivalent_bias <- readRDS("bipolar_bias/clin_ambivalent_bias_combined_sumstats.R")

#specify the lavaan model
ambivalent_correlation_traits <- c("MDD", "PTSD", "Neu", "AN", "ASD", "ADHD", "Smoking")


model <- paste0(
"Bias =~ NA*",paste(ambivalent_correlation_traits, collapse = " + "),"
Liability =~ NA*Clin + ",paste(ambivalent_correlation_traits, collapse = " + "),"
Bias ~~ 0*Liability
Clin ~~ 0*Clin
Bias ~~ 1*Bias
Liability ~~ 1*Liability
MDD ~~ a*MDD
a> 0.0001
    
Bias ~ SNP
Liability ~ SNP")
    


#create the path to the object you want the output to store
bias_factor_gwas_path <- "bipolar_bias/GWAS_Bias_factor.R"
start_time <- Sys.time()

cat("Starting commonfactorGWAS run at", Sys.time(), "\n")

bias_factor_bip_gwas <- userGWAS(covstruc = covstruct_ambivalent_correlations, 
                                              SNPs = sumstats_ambivalent_bias,
                                              model = model,
                                              sub = c("Bias ~ SNP"),
                                              toler = 1e-70,
                                              parallel = TRUE, 
                                              cores = 10,
                                              Q_SNP = FALSE
                                     )
saveRDS(bias_factor_bip_gwas, file = bias_factor_gwas_path)


runtime_mins <- as.numeric(difftime(Sys.time(), start_time, units="mins"))
cat("Model finished in", round(runtime, 2), "minutes.\n")
cat("Finished at", Sys.time(), "\n")