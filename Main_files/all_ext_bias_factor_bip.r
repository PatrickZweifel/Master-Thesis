#load the non-munged summary stats for which the error are missing and the original summary stats
original_sumstats <- c("daner_adhd_meta_filtered_NA_iPSYCH23_PGC11_sigPCs_woSEX_2ell6sd_EUR_Neff_70.meta.gz",  "iPSYCH-PGC_ASD_Nov2017.gz",  "PGC3_SCZ_wave3.european.autosome.public.v3.vcf.tsv.gz",  "pgcAN2.2019-07.vcf.tsv.gz",  "pgc-panic2019.vcf.tsv.gz")
original_path <- "/local1/scratch/pazweifel/sumstats_ambits/original_sumstats/"
original_paths <- paste0(original_path, original_sumstats)

bad_sumstats <- c("ADHD.txt", "ASD.txt", "SCZ.txt", "AN.txt", "PD.txt")
bad_path <- "/local1/scratch/pazweifel/sumstats_ambits/non_munged/"
bad_paths <- paste0(bad_path, bad_sumstats)

names(original_paths) <- c("ADHD", "ASD", "SCZ", "AN", "PD")
names(bad_paths) <- names(original_paths)
traits <- names(bad_paths)
names(traits) <- traits

#sumstats with overhead
overhead_sumstats <- c("SCZ", "AN", "PD")

original_datasets <- lapply(traits, function(trait) {
    if (trait %in% overhead_sumstats){
        fread(original_paths[[trait]], skip = "CHROM", data.table = FALSE)
    }
    else{
        fread(original_paths[[trait]], data.table = FALSE)
    }
    }
)

bad_datasets <- lapply(traits, function(trait) {
    fread(bad_paths[[trait]], data.table = FALSE)
    })

#check how it is named for the ones not calling it SNP
head(original_datasets$SCZ)
# its "ID", so rename it
original_datasets_curated <- lapply(traits, function(trait){
    if (!("SNP" %in% colnames(original_datasets[[trait]]))){
        rename(original_datasets[[trait]], SNP = ID)
    }
    else{
        original_datasets[[trait]]
    }
    }
)

for (trait in original_datasets_curated){
    print("SNP" %in% colnames(trait))
    }

#now match on the snp column and combined the bad_one with the original one, so that the new df now also contains the error column
good_sumstats <- lapply(traits, function(trait) {
    inner_join(x = bad_datasets[[trait]], y = select(original_datasets_curated[[trait]], SNP, SE), by = join_by(SNP))
    }
)

for (trait in good_sumstats){
    print(sum(is.na(trait$SE)))
    }

for (trait in good_sumstats){
    print(nrow(trait))
    }

for (trait in bad_datasets){
    print(nrow(trait))
    }

head(good_sumstats$AN)

new_folder_path <- "/local1/scratch/pazweifel/sumstats_ambits/gsem_ready/"

for (trait in traits){
    write.table(
        good_sumstats[[trait]],
        file = paste0(new_folder_path, trait, ".txt"),
        sep = "\t",
        col.names = TRUE,
        row.names = FALSE,
        quote = FALSE)
    }

#cleanup of the clinical dataset..

clinical_bd_bad <- fread("/local1/scratch/pazweifel/heritability_analysis/ldsc_bipolar_comclin/bip2024_eur_clinical_no23andMe_wbeta_use_neff.txt", data.table = FALSE)

head(clinical_bd_bad)

clinical_bd_bad_selected <- clinical_bd_bad %>% select(-c(Nca, Nco, N, BETA)) %>% rename(N = NEFF)

head(clinical_bd_bad_selected)

clinical_bd_good <- clinical_bd_bad_selected %>% filter(!str_detect(SNP, "^[0-9]"))

head(clinical_bd_good)
nrow(clinical_bd_good)
nrow(clinical_bd_bad)

clinical_path <- "/local1/scratch/pazweifel/sumstats_ambits/gsem_ready/BD_Clin.txt"

if (!file.exists(clinical_path)){
    write.table(x = clinical_bd_good,
                file = clinical_path,
                sep = "\t",
                col.names = TRUE,
                row.names = FALSE,
                quote = FALSE)
}

require(data.table)
comm_bd_bad <- fread("/local1/scratch/pazweifel/heritability_analysis/ldsc_bipolar_comclin/bip2024_eur_community_no23andMe_wbeta_use_neff.txt", data.table = FALSE)

head(comm_bd_bad)

comm_bd_bad_selected <- comm_bd_bad %>% select(-c(Nca, Nco, N, BETA)) %>% rename(N = NEFF)

head(comm_bd_bad_selected)

comm_bd_godd<- comm_bd_bad_selected %>% filter(!str_detect(SNP, "^[0-9]"))

comm_path <- "/local1/scratch/pazweifel/sumstats_ambits/gsem_ready/Comm.txt"

if (!file.exists(comm_path)){
    write.table(x = comm_bd_godd,
                file = comm_path,
                sep = "\t",
                col.names = TRUE,
                row.names = FALSE,
                quote = FALSE)
}

#calcualte one big covariance matrix involving the clinical bipolar trait and all external ones
covstruct_bias_selected <- here::here("matrices/Clinical_selected_external.R")

selected_traits <- c("ADHD", "ALC", "AN", "ASD", "BMI", "CAD", "Chrono", "CRP", "Height", "Neu", "PartEM", "PartMH", "PartSF", "PD", "PTSD", "SCZ", "Smoking", "T2D")

pop_prev <- c(bip_prevalences["Clin"], other_prevalences[selected_traits])

samp_prev <- as.vector(ifelse(is.na(pop_prev), NA, 0.5))

paths <- c(samples[1])
other_basepath <- "/local1/scratch/pazweifel/sumstats_ambits/munged/"

for (name in selected_traits){
    paths <- c(paths, paste0(other_basepath, name, ".sumstats.gz"))
}

    
if (!file.exists(covstruct_bias_selected)){
    covstruct_bias_selected_ldsc <- ldsc(
        traits = paths,
        trait.names = names(pop_prev),
        sample.prev = samp_prev,
        population.prev = pop_prev,
        ld = ld, 
        wld = weights,
        stand = TRUE
        )
    dput(covstruct_bias_selected_ldsc, covstruct_bias_selected, control = c("all", "digits17"))
}

covstruct.bias_selected <- dget("matrices/Clinical_selected_external.R")

usermodel(covstruc = bias_all.covstruct, estimation = "DWLS", model = model, CFIcalc = TRUE, imp_cov = TRUE)

sumstats_path <- "/local1/scratch/pazweifel/sumstats_ambits/gsem_ready/"
ref <- "/local1/scratch/pazweifel/jolien_paper_sumstats/reference.1000G.maf.0.005.txt"

traits <- c("Clin", "ADHD", "ALC", "AN", "ASD", "BMI", "CAD", "Chrono", "CRP", "Height", "Neu", "PartEM", "PartMH", "PartSF", "PD", "PTSD", "SCZ", "Smoking", "T2D")
se.logit <- c(T,T,F,T,T,F,T,F,F,F,F,T,T,T,T,T,T,F,T)
linprob <- c(F,F,T,F,F,F,F,F,F,F,T,F,F,F,F,F,F,F,F)
ols <- c(F,F,F,F,F,T,F,T,T,T,F,F,F,F,F,F,F,T,F)

path_clin <- paste0(sumstats_path,"BD_Clin.txt")
path_rest <- paste0(sumstats_path, traits[2:length(traits)],".txt")
all_paths <- c(path_clin, path_rest)

#run the sumstats function
dir.create("bipolar_bias", showWarnings = FALSE)

bip_bias_sumstats_path <- "bipolar_bias/clin_bias_combined_sumstats.R"

if (!file.exists(bip_bias_sumstats_path)){
    bias_sumstats <- sumstats(files = all_paths,
             ref = ref,
             trait.names = traits,
             se.logit = se.logit,
             OLS = ols, 
             linprob = linprob,
             parallel = TRUE,
             cores = 10)
    
    saveRDS(bias_sumstats, file = bip_bias_sumstats_path) 
    }
     

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
ld <- "/local1/hdata/REF/eur_ref_ld_chr"
weights <- "/local1/hdata/REF/eur_w_ld_chr"
samples <- c("/local1/scratch/pazweifel/heritability_analysis/ldsc_correlation_bipolar/sumstats/bip2024_eur_clinical_no23andMe_wbeta_use_neff.sumstats.gz",
             "/local1/scratch/pazweifel/heritability_analysis/ldsc_correlation_bipolar/sumstats/bip2024_eur_community_no23andMe_wbeta_use_neff.sumstats.gz",
             "/local1/scratch/pazweifel/heritability_analysis/ldsc_correlation_bipolar/sumstats/pgc-bip2021-BDI_cleaned_neff.sumstats.gz",
             "/local1/scratch/pazweifel/heritability_analysis/ldsc_correlation_bipolar/sumstats/pgc-bip2021-BDII_cleaned_neff.sumstats.gz")
sample_names <-c("Clin", "Comm", "BDI", "BDII")
other_basepath <- "/local1/scratch/pazweifel/sumstats_ambits/munged/" #phenotype.sumstats.gz

### create the new covstruct including the bias factor sumstats
covstruct_liability_comm_factor <- "matrices/neu_ext_comm_liability_factor.R"
necessary_traits <- c("Clin","Comm","Neu")

combined_paths <- c(samples[1:2], "/local1/scratch/pazweifel/sumstats_ambits/munged/Neu.sumstats.gz")

pop.prev <- c(all_prevalences[necessary_traits])

samp.prev <- c(ifelse(is.na(pop.prev), NA, 0.5))

if (!file.exists(covstruct_liability_comm_factor)){
    covstruct_neu_ext_comm_liab <- ldsc(
        traits = combined_paths,
        trait.names = c(necessary_traits),
        sample.prev = samp.prev,
        population.prev = pop.prev,
        ld = ld,
        wld = weights,
        stand = TRUE
    )
    #store the file in the previously defined path, now the path is, contorl tell how to deparse
    dput(covstruct_neu_ext_comm_liab, covstruct_liability_comm_factor, control = c("all", "digits17"))
}

gwas_by_sub_model <- "
Liability_clin =~ NA*Clin
Liability_other =~ NA*{internal_trait}
Bias =~ NA*Bias_GWAS + {internal_trait}

Liability_clin ~~ 0*Bias
Liability_other ~~ 0*Bias
Clin ~~ 0*Clin
{internal_trait} ~~ 0*{internal_trait}
Bias_GWAS ~~ 0*Bias_GWAS
Clin ~~ Bias_GWAS
Bias_GWAS ~~ 0*{internal_trait}
Clin ~~ 0*{internal_trait}
"

sumstats_path <- "/local1/scratch/pazweifel/sumstats_ambits/gsem_ready/"
ref <- "/local1/scratch/pazweifel/jolien_paper_sumstats/reference.1000G.maf.0.005.txt"

necessary_traits <- c("Clin","Comm","Neu")
se.logit <- c(T,T,F)
linprob <- c(F,F,T)
ols <- c(F,F,T)

path_clin <- paste0(sumstats_path,"BD_Clin.txt")
path_rest <- paste0(sumstats_path, necessary_traits[2:length(necessary_traits)],".txt")
all_paths <- c(path_clin, path_rest)

#run the sumstats function
dir.create("bipolar_bias", showWarnings = FALSE)

bip_bias_sumstats_path <- "bipolar_bias/clin_comm_neu_sumstats.R"

if (!file.exists(bip_bias_sumstats_path)){
    bias_sumstats <- sumstats(files = all_paths,
             ref = ref,
             trait.names = necessary_traits,
             se.logit = se.logit,
             OLS = ols, 
             linprob = linprob,
             parallel = TRUE,
             cores = 10)
    
    saveRDS(bias_sumstats, file = bip_bias_sumstats_path) 
    }
     
