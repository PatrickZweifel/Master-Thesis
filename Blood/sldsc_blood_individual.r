library(dplyr)
library(readr)
library(tidyr)
library(stringr)
library(ggplot2)
library(microshades)
library(openxlsx)
library(GenomicSEM)
library(lavaanPlot)
library(lavaan)
library(semPlot)

#load the metadata file
setwd(dir = "/local1/home/pazweifel/src/blood_analysis")
metadata <- read_tsv("blood_names.tsv", col_names = TRUE)

sample_names <- metadata$`Standard abbreviation`

ld <- "/local1/hdata/REF/eur_ref_ld_chr"
weights <- "/local1/hdata/REF/eur_w_ld_chr"
basepath <- "/local1/scratch/pazweifel/sumstats_blood/munged_sumstats/" #phenotype.sumstats.gz
sample_prevalences <- rep(NA, length(sample_names))
population_prevalences <- rep(NA, length(sample_names))

samples <- c()

for (trait in sample_names){
    samples <- c(samples, paste0(basepath, trait, ".sumstats.gz"))
    }

names(samples) <- sample_names
metadata <- rename(metadata, ABBRV = `Standard abbreviation`, Trait = `Long name`, N = `Number of participants contributing to GWAS`)
samples_fsc_ssc <- samples[grepl("SSC$|FSC$", names(samples))]


annotation_path_general <- "/local1/scratch/pazweifel/single_cell_1k1k_respective/"
a <- list.dirs(annotation_path_general, full.names = FALSE)
directory_names <- a[2:length(a)]

ld_complete <- "/local1/hdata/REF/baseline_v2.2/"
weights <- "/local1/hdata/REF/eur_w_ld_chr/"
frq <- "/local1/scratch/pazweifel/heritability_analysis/1000G_Phase3_frq/"
dir.create("sldsc_individual_annotation_covstructs", showWarnings = FALSE)

for (annotation in directory_names){
ld_vector <- c(ld_complete, paste0(annotation_path_general, annotation, "/"))

output_path <- paste0("sldsc_individual_annotation_covstructs/S_LDSC_blood_factors_fsc_ssc_",annotation,".rds")


sldsc_blood <- s_ldsc(
                    traits = samples_fsc_ssc,
                    sample.prev = sample_prevalences[1:length(samples_fsc_ssc)],
                    population.prev = population_prevalences[1:length(samples_fsc_ssc)],
                    ld = ld_vector,
                    wld = weights,
                    frq = frq,
                    trait.names = gsub("-", "", names(samples_fsc_ssc))
)
saveRDS(sldsc_blood, file = output_path)
}
