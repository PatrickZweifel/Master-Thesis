library(dplyr)
library(readr)
library(tidyr)
library(stringr)
library(ggplot2)
library(microshades)
library(openxlsx)
library(GenomicSEM)
library(data.table)

setwd(dir = "/local1/home/pazweifel/src/blood_analysis")
#load the required files
covstruct_blood <- dget("matrices/FSC_SSC_factor4.R")
blood_sumstats <- readRDS("commonfactor/fsc_ssc_factor4_sumstats.R")

hapmap3_snplist <- fread("/local1/scratch/pazweifel/heritability_analysis/w_hm3_snplist.gz", data.table = FALSE)
hapmap3_snps <- hapmap3_snplist$SNP
#filter the original sumstats file to a hapmap3 subset
blood_sumstats_hapmap3_subset <- blood_sumstats %>% filter(SNP %in% hapmap3_snps)
#specify the lavaan model
model <- 'F1 =~NA*RBCSSC+RBCFSC+RETFSC+IRFFSC
F2 =~NA*NEFSC+EOSSC+EOFSC+BASOFSC
F3 =~NA*PLTSSC+PLTFSC
F4 =~NA*EOSSC+EOFSC+BASOFSC+LYSSC+LYFSC
F1~~F2
F1~~F3
F1~~F4
F2~~F3
F2~~F4
F3~~F4
F1 ~ SNP
F2 ~ SNP
F3 ~ SNP
F4 ~ SNP
EOFSC ~~ a*EOFSC
a > 0.0001'


#create the path to the object you want the output to store
mental_common_factor_path <- "commonfactor/gwas_sumstats_factor4_qsnp.R"
start_time <- Sys.time()

cat("Starting commonfactorGWAS run at", Sys.time(), "\n")

mental_common_factor_gwas <- userGWAS(covstruc = covstruct_blood, 
                                              SNPs = blood_sumstats_hapmap3_subset,
                                              model = model,
                                              sub = c("F1~SNP", "F2~SNP", "F3~SNP", "F4~SNP"),
                                              toler = 1e-70,
                                              parallel = TRUE, 
                                              cores = 20,
                                              Q_SNP = TRUE
                                     )
saveRDS(mental_common_factor_gwas, file = mental_common_factor_path, compress = "gzip")


runtime_mins <- as.numeric(difftime(Sys.time(), start_time, units="mins"))
cat("Model finished in", round(runtime_mins, 2), "minutes.\n")
cat("Finished at", Sys.time(), "\n")
