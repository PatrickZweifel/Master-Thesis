library(dplyr)
library(readr)
library(tidyr)
library(stringr)
library(ggplot2)
library(microshades)
library(openxlsx)
library(GenomicSEM)

#sample prevalences can be set to 0.5 as N for the munged files has been converted to effective N. It could also be odne by specifying the correct sample prevalence when convertint to the liability scale, but for compatibility purposes the neffe variant has been applied
bip_prevalences <- c("Clin" = 0.02, "Comm" = 0.02, "BDI" = 0.01, "BDII" = 0.01)
bip_names <- names(bip_prevalences)

other_prevalences <- c(
  "ADHD" = 0.03, "ALC" = 0.16, "AN" = 0.01,
  "ASD" = 0.01, "MDD" = 0.15, "BMI" = NA,
  "CAD" = 0.05, "Chrono" = NA, "CRP" = NA,
  "Height" = NA, "Neu" = NA, "PD" = 0.03,
  "PTSD" = 0.3, "SCZ" = 0.01, "Smoking" = NA,
  "T2D" = 0.06, "PartEM" = 0.6, "PartMH" = 0.3, "PartSF" = 0.5
)
other_names <- names(other_prevalences)

all_prevalences <- c(bip_prevalences, other_prevalences)
all_names <- names(all_prevalences)

#"BIP" = 0.01

#load the paths of necessary documents as not all of them are in this directory

samples <- c("/local1/scratch/pazweifel/heritability_analysis/ldsc_correlation_bipolar/sumstats/bip2024_eur_clinical_no23andMe_wbeta_use_neff.sumstats.gz",
             "/local1/scratch/pazweifel/heritability_analysis/ldsc_correlation_bipolar/sumstats/bip2024_eur_community_no23andMe_wbeta_use_neff.sumstats.gz",
             "/local1/scratch/pazweifel/heritability_analysis/ldsc_correlation_bipolar/sumstats/pgc-bip2021-BDI_cleaned_neff.sumstats.gz",
             "/local1/scratch/pazweifel/heritability_analysis/ldsc_correlation_bipolar/sumstats/pgc-bip2021-BDII_cleaned_neff.sumstats.gz")
sample_names <-c("Clin", "Comm", "BDI", "BDII")
other_basepath <- "/local1/scratch/pazweifel/sumstats_ambits/munged/" #phenotype.sumstats.gz

path_single_cell_annotations <- "/local1/scratch/pazweifel/single_cell_Duncan_og_combined/"
ld_complete <- "/local1/hdata/REF/baseline_v2.2/"
weights <- "/local1/hdata/REF/eur_w_ld_chr/"
frq <- "/local1/scratch/pazweifel/heritability_analysis/1000G_Phase3_frq/"
weights <- "/local1/hdata/REF/eur_w_ld_chr/"


ld_vector <- c(ld_complete, path_single_cell_annotations)
sldsc_traits <- c("Clin")

output_path <- "S_LDSC_clin.RData"
samples_clin <- c(samples[1])

sldsc_pop_prevalences <- all_prevalences[sldsc_traits]
sldsc_sample_prevalences <- as.vector(ifelse(is.na(sldsc_pop_prevalences), NA, 0.5))


sldsc_liab_factor <- s_ldsc(
                    traits = samples_clin,
                    sample.prev = sldsc_sample_prevalences,
                    population.prev = sldsc_pop_prevalences,
                    ld = ld_vector,
                    wld = weights,
                    frq = frq,
                    trait.names = sldsc_traits
)

saveRDS(sldsc_liab_factor, file = output_path)









