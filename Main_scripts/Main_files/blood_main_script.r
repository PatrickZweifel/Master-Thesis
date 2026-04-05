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
library(ggrepel)

#load the metadata file
setwd(dir = "/local1/home/pazweifel/src/blood_analysis")
metadata <- read_tsv("blood_names.tsv", col_names = TRUE)

metadata$source_DOI <- "10.1038/s41467-023-40679-y"

head(metadata)

write_csv(metadata, "/local1/home/pazweifel/plots/blood/sample_data.csv")

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

metadata$ABBRV <- metadata$ABBRV %>% gsub(c("-", "#", "%", "(",")"), c("", "hash", "prctg", "_" ,"_"),.) %>% gsub("#","hash",.) %>% gsub("%","prctg",.) %>% gsub("\\(", "_", .) %>% gsub("\\)","_",.)

sample_names_all <- metadata$ABBRV

samples_fsc_ssc <- samples[grepl("SSC$|FSC$", names(samples))]

dir.create("matrices", showWarnings = FALSE)

covstruct_blood <- here::here("matrices","FSC_SSC.R")

if (!file.exists(covstruct_blood)){
    ldsc_blood <- ldsc(
                    traits = samples_fsc_ssc,
                    sample.prev = sample_prevalences[1:length(samples_fsc_ssc)],
                    population.prev = population_prevalences[1:length(samples_fsc_ssc)],
                    ld = ld,
                    wld = wld,
                    trait.names = gsub("-", "", names(samples_fsc_ssc))
        )
    dput(ldsc_blood, covstruct_blood, control = c("all", "digits17"))
}

ssc_fsc_covstruct <- dget(here::here("matrices","FSC_SSCR"))

covstruct_full <- here::here("matrices","all.R")

if (!file.exists(covstruct_full)){
    ldsc_blood <- ldsc(
                    traits = samples,
                    sample.prev = sample_prevalences,
                    population.prev = population_prevalences,
                    ld = ld,
                    wld = wld,
                    trait.names = sample_names_all,
                    stand = TRUE)
    saveRDS(object = ldsc_blood, file = covstruct_full)
}

#smooth the S matrix
require(Matrix)

ssc_fsc_covmmatrix_smoothed <- as.matrix((nearPD(ssc_fsc_covstruct$S, corr = FALSE))$mat)

#run EFA for different amount of factors specified
#preload the approximated number of model structures for memory efficiency
efa_list <- list(0,0,0,0,0)

require(stats)
#run the factanal function fro different amounts of factors
for (i in 1:length(efa_list)){
    EFA <- factanal(covmat = ssc_fsc_covmmatrix_smoothed, factors = i, rotation = "promax")
    efa_list[[i]] <- EFA

}

lavaan_structures_from_efa <- function(efa_list){
    #initialize the lavaan structure list
    lavaan_structure_list <- rep(list(0),length(efa_list))

    #iterate through the efa list to build the lavaan structures for each model
    for (i in 1:length(efa_list)){
        factor_lavaan_list <- rep(list(0),i)
    
        #extract the relevant trait on factor loadings and define the lavaan structure
        
        for (factor_index in 1:i){
            #for each factor loading column of an i-factor model, extract the names of the relevant traits (loading > 0.2)
            #and paste them together usign lavaan syntax
            loadings <- efa_list[[i]]$loadings[,factor_index]
            relevant_loadings_traits <- names(loadings[loadings > 0.2])
            factor_lavaan_list[[factor_index]] <- paste0("F",as.character(factor_index),
                                                 " =~NA*", 
                                                 paste(relevant_loadings_traits, collapse = "+")
                                                 )
            }
    
        #define the correlations between the individual factors if i > 1
             
        if (i > 1){
            factor_numbers <- seq(1,i)
            #get all the pairwise unique combinations of factors
            combination_matrix <- t(combn(factor_numbers, 2))
            print(nrow(combination_matrix))
            factor_correlation_list <- rep(list(0), nrow(combination_matrix))
            for (combination in 1:nrow(combination_matrix)){
                f_left <- as.character(combination_matrix[combination, 1])
                f_right <- as.character(combination_matrix[combination, 2])
        
                factor_correlation_list[[combination]] <- paste0("F",f_left,"~~","F",f_right)
            }
            #collapse the lists and create one lavaan structure object
            collapsed_factor_loadings <- paste(factor_lavaan_list, collapse = "\n")
            collapsed_factor_correlations <- paste(factor_correlation_list, collapse = "\n")
        
            lavaan_structure_list[[i]] <- paste(collapsed_factor_loadings, collapsed_factor_correlations, sep = "\n")
        }
        else{
            lavaan_structure_list[[i]] <- paste(factor_lavaan_list, collapse = "\n")
        }
    }
    lavaan_structure_list
}

lavaan_structure_list <- lavaan_structures_from_efa(efa_list)

lavaan_structure_list[[3]]

#run the cfa for all efa models using the lavaan strucutres created before and the covariance matrix
cfa_list <- lapply(lavaan_structure_list, function(lavaan_structure){
    usermodel(ssc_fsc_covstruct, estimation = "DWLS", model = lavaan_structure, CFIcalc = TRUE, std.lv = TRUE, imp_cov = TRUE)
    }
)

modelfit_list <- bind_rows(lapply(cfa_list, function(model) model$modelfit))

modelfit_list

write_csv(modelfit_list, "/local1/home/pazweifel/plots/blood/efa_modelfits.csv")

semPlotModel_GSEM <- function(gsem.object=cfa_bipolar){ 
        object <- gsem.object$results
        object$free=0
        numb=1:length(which(object$op!="~~"))
        object$free[which(object$op!="~~")]=numb
        varNames <- lavaanNames(object, type = "ov")
        factNames <- lavaanNames(object, type = "lv")
        factNames <- factNames[!factNames %in% varNames]
        n <- length(varNames)
        k <- length(factNames)
        if (is.null(object$label)) 
          object$label <- rep("", nrow(object))
        semModel <- new("semPlotModel")
        object$std <- object[,"STD_Genotype"]
        object$est <- object[,"Unstand_Est"]
        if (is.null(object$group)) 
          object$group <- ""
        semModel@Pars <- data.frame(label = object$label, lhs = ifelse(object$op == "~" | object$op == "~1", object$rhs, object$lhs), edge = "--", 
                                    rhs = ifelse(object$op == "~" | object$op == "~1", object$lhs, object$rhs), est = object$est, std = object$std, std = NA, group = object$group, 
                                    fixed = object$free==0, par = object$free, stringsAsFactors = FALSE)
        semModel@Pars$edge[object$op == "~~"] <- "<->"
        semModel@Pars$edge[object$op == "~*~"] <- "<->"
        semModel@Pars$edge[object$op == "~"] <- "~>"
        semModel@Pars$edge[object$op == "=~"] <- "->"
        semModel@Pars$edge[object$op == "~1"] <- "int"
        semModel@Pars$edge[grepl("\\|", object$op)] <- "|"
        semModel@Thresholds <- semModel@Pars[grepl("\\|", semModel@Pars$edge), 
                                             -(3:4)]
        semModel@Pars <- semModel@Pars[!object$op %in% c(":=", "<", 
                                                         ">", "==", "|", "<", ">"), ]
        semModel@Vars <- data.frame(name = c(varNames, factNames), 
                                    manifest = c(varNames, factNames) %in% varNames, exogenous = NA, 
                                    stringsAsFactors = FALSE)
        semModel@ObsCovs <- list()
        semModel@ImpCovs <- list()
        semModel@Computed <- FALSE
        semModel@Original <- list(object)
        return(semModel)
 }

for (i in 1:length(cfa_list)){
    #define the two sem objects
    current_fit <- semPlotModel_GSEM(gsem.object = cfa_list[[i]])

    #extract the std estimate and the corresponding se to plot them both, for common and independent model
    se <- cfa_list[[i]]$results$STD_Genotype_SE
    est <- cfa_list[[i]]$results$STD_Genotype
    est.se <- paste0(round(as.numeric(est), 2),"\n (",round(as.numeric(se), 2),")")
    print(est.se)

    pdf(paste0("/local1/home/pazweifel/plots/blood/fsc_ssc/",as.character(i),"_factor_model.pdf"), width = 30, height = 20)
    semPaths(current_fit, whatLabels = "std", layout = "tree", edge.color = "black", sizeMan = 5, sizeLat = 5, nCharNodes = 0, residuals = TRUE, style = "mx",
        edge.label.cex = 0.2,
        esize = 0.5,
        edgeLabels = est.se)
    dev.off()
}  

lavaan_structure_list[[4]]

fsc_ssc_4_factor_model <- 'F1 =~NA*RBCSSC+RBCFSC+RETFSC+IRFFSC
F2 =~NA*NEFSC+EOSSC+EOFSC+BASOFSC
F3 =~NA*PLTSSC+PLTFSC
F4 =~NA*EOSSC+EOFSC+BASOFSC+LYSSC+LYFSC
F1~~F2
F1~~F3
F1~~F4
F2~~F3
F2~~F4
F3~~F4
EOFSC ~~ a*EOFSC
a > 0.0001'

factor4_adjusted.fit <- usermodel(ssc_fsc_covstruct, estimation = "DWLS", model = fsc_ssc_4_factor_model, CFIcalc = TRUE, std.lv = TRUE, imp_cov = TRUE)

factor4_adjusted.fit

write_csv(factor4_adjusted.fit$modelfit, "/local1/home/pazweifel/plots/blood/modelfit_4.csv")
write_csv(factor4_adjusted.fit$results, "/local1/home/pazweifel/plots/blood/results_4.csv")

#extract the names of the remaining traits and subset the sample vector further to only include the relevant traits
names(samples_fsc_ssc) <- gsub("-", "", names(samples_fsc_ssc))
traits_to_keep_factor_4 <- names(samples_fsc_ssc) %in% rownames(factor4_adjusted.fit$resid_cov[[1]])
samples_fsc_ssc_factor_4 <- samples_fsc_ssc[traits_to_keep_factor_4]

names(samples_fsc_ssc_factor_4)

covstruct_blood <- here::here("matrices","FSC_SSC_factor4.R")

if (!file.exists(covstruct_blood)){
    ldsc_blood <- ldsc(
                    traits = samples_fsc_ssc_factor_4,
                    sample.prev = sample_prevalences[1:length(samples_fsc_ssc_factor_4)],
                    population.prev = population_prevalences[1:length(samples_fsc_ssc_factor_4)],
                    ld = ld,
                    wld = wld,
                    trait.names = names(samples_fsc_ssc_factor_4)
                    
        )
    dput(ldsc_blood, covstruct_blood, control = c("all", "digits17"))
}

#as soon as gwas done delete the one above

covstruct_blood <- here::here("matrices","FSC_SSC_factor4_stand.R")

if (!file.exists(covstruct_blood)){
    ldsc_blood <- ldsc(
                    traits = samples_fsc_ssc_factor_4,
                    sample.prev = sample_prevalences[1:length(samples_fsc_ssc_factor_4)],
                    population.prev = population_prevalences[1:length(samples_fsc_ssc_factor_4)],
                    ld = ld,
                    wld = wld,
                    trait.names = names(samples_fsc_ssc_factor_4),
                    stand = TRUE
                    
        )
    dput(ldsc_blood, covstruct_blood, control = c("all", "digits17"))
}

covstruct_blood_fsc_ssc_factor4 <- dget(here::here("matrices","FSC_SSC_factor4.R"))

covstruct_blood_fsc_ssc_factor4$S

factor4_adjusted_newcovstruct.fit <- usermodel(covstruct_blood_fsc_ssc_factor4, estimation = "DWLS", model = fsc_ssc_4_factor_model, CFIcalc = TRUE, std.lv = TRUE, imp_cov = TRUE)                                   

factor4_adjusted_newcovstruct.fit$resid_cov[[1]] - factor4_adjusted.fit$resid_cov[[1]]

factor4_adjusted_newcovstruct.fit

#extract the names of the traits which are relevant from the initial sample
sample_names_woline <- gsub("-", "", sample_names)
samples_to_keep_commonfactor <- sample_names_woline %in%  rownames(factor4_adjusted.fit$resid_cov[[1]])
commonfactor_samples <- sample_names[samples_to_keep_commonfactor]
names(commonfactor_samples) <- commonfactor_samples

basepath_unmunged_sumstats <- "/local1/scratch/pazweifel/sumstats_blood/sumstats_gsemcompatible/"

sumstats_unmunged_commonfactor <- sapply(commonfactor_samples, function(trait) paste0(basepath_unmunged_sumstats, trait, ".tsv"))

#continuous traits, so OLS = TRUE, se.logit = FALSE 
OLS <- rep(TRUE, length(sumstats_unmunged_commonfactor))
se.logit <- rep(FALSE, length(sumstats_unmunged_commonfactor))                                        
#reference file for maf
ref <- "/local1/scratch/pazweifel/jolien_paper_sumstats/reference.1000G.maf.0.005.txt"

#run the sumstats functio
dir.create("commonfactor", showWarnings = FALSE)

sumstats_path <- here::here("commonfactor","fsc_ssc_factor4_sumstats.R")

if (!file.exists(sumstats_path)){
    combined_sumstats <- sumstats(files = sumstats_unmunged_commonfactor,
                                  ref = ref,
                                  trait.names = names(samples_fsc_ssc_factor_4),
                                  se.logit = se.logit,
                                  OLS = OLS,
                                  parallel = TRUE,
                                  cores = 30)
    saveRDS(combined_sumstats, sumstats_path)
}
                                  
                                  

library(data.table)
hapmap3_snplist <- fread("/local1/scratch/pazweifel/heritability_analysis/w_hm3_snplist.gz", data.table = FALSE)
hapmap3_snps <-hapmap3_snplist$SNP

head(hapmap3_snps)

library(qqman)

blood_gwas <- readRDS("commonfactor/gwas_sumstats_factor4.R")

head(blood_gwas[[1]])

blood_gwas_filtered <- list()

for (i in 1:length(blood_gwas)){
    blood_gwas_filtered_temp <- blood_gwas[[i]] %>% filter(error == 0 & warning == 0 & Pval_Estimate < 0.1)
    blood_gwas_filtered[[i]] <- blood_gwas_filtered_temp
}

nrow(blood_gwas[[1]])
nrow(blood_gwas_filtered[[1]])

#factor 1
manhattan(blood_gwas_filtered[[1]], chr = "CHR", bp = "BP", p = "Pval_Estimate", snp = "SNP")

qq(filter(blood_gwas[[1]], warning == 0 & error == 0)$Pval_Estimate)

plot_path <- "/local1/home/pazweifel/plots/blood/fsc_ssc/gwas_factor_"

for (i in 1:length(blood_gwas_filtered)){
  manhattan_path <- paste0(plot_path, i,"_manhattan.pdf")
  qq_path <- paste0(plot_path, i, "_qq.pdf")

  #manhattan
  pdf(file = manhattan_path, width = 10, height = 12)
  manhattan(blood_gwas_filtered[[i]], chr = "CHR", bp = "BP", p = "Pval_Estimate", snp = "SNP")
  dev.off()

  #qq
  pdf(file = qq_path, width = 10, height = 12)
  qq(filter(blood_gwas[[i]], warning == 0 & error == 0)$Pval_Estimate)
  dev.off()
}

nrow(covstruct_blood_fsc_ssc_factor4$V)

covstruct_blood_fsc_ssc_factor4_stand <- dget(here::here("matrices","FSC_SSC_factor4_stand.R"))

sstand <- covstruct_blood_fsc_ssc_factor4_stand$S_Stand
vstand <- covstruct_blood_fsc_ssc_factor4_stand$V_Stand

#extract the standard error
vstand_diag <- sapply(c(1:nrow(covstruct_blood_fsc_ssc_factor4$V)), function(index) vstand[index, index])
                      
#copy the structure of the sstand matrix
error_matrix <- sstand
                      
#overwrite lower triangle with vstand values
error_matrix[lower.tri(error_matrix, diag = TRUE)] <- vstand_diag
                      
#transpose the matrix and copy the values to the upper triangle, making it symmetric
error_matrix[upper.tri(error_matrix)] <- t(error_matrix)[upper.tri(error_matrix)]
error_matrix <- sqrt(error_matrix)
#error_matrix
                      
#vectorize the error matrix so that you can use it for the ggplot
error_vectorized <- c(error_matrix)

#initialize the dataframe for the matrix values of both matrices
tempnames <- dimnames(sstand)[[2]]
expanded_df <- expand.grid(x = tempnames, y = tempnames)

#vectorize the standard matrix and store matrix values in df columns
sstand_vectorized <- c(sstand)
expanded_df$values <- round(sstand_vectorized, 2)
expanded_df$se <- round(error_vectorized, 2)
#expanded_df
#create the ggplot
correlation_plot <- ggplot(expanded_df, aes(x = x, y = y, fill = values)) +
scale_fill_gradient(low = "darkorchid", high = "orange", limits = c(min(expanded_df$values), max(expanded_df$values))) +
theme_minimal(base_size = 25) +
geom_tile() +
theme(axis.text.x = element_text(angle = 45, hjust = 1), axis.title = element_blank()) +
labs(fill = expression(r[g])) +
geom_text(aes(label = values)) +
geom_text(aes(label = paste0("(",se,")"), vjust = 2))

correlation_plot

ggsave("/local1/home/pazweifel/plots/blood/fsc_ssc/correlation_plot.pdf", correlation_plot, device = pdf, width = 20, height = 17)

write_csv(expanded_df, "/local1/home/pazweifel/plots/blood/correlations.csv")

#extract the heritabilities and their se
temp <- covstruct_blood_fsc_ssc_factor4_stand
S <- temp$S
V <- temp$V
herit_base <- sapply(c(1:nrow(S)), function(index) S[index,index])
                  
index_se <- rep(1, nrow(S))

for (i in 2:nrow(S)){
    index_se[i] <- index_se[i-1] + nrow(S) - i + 2
}
                     
se_base <- sapply(index_se, function(index) sqrt(V[index, index]))
           
names <- dimnames(S)[[2]][1:nrow(S)]
                
heritabilities <- tibble(trait = names, heritability = herit_base, se = se_base)

heritabilities


heritabilities$trait <- factor(heritabilities$trait, levels = unique(heritabilities$trait))


heritabilities_blood_paper <- ggplot(heritabilities, aes(x = trait, y = heritability)) + 
geom_col(color = "black", fill = "orange") +
theme_bw() +
geom_errorbar(aes(ymin = heritability - se, ymax = heritability + se), width = 0.3) +
theme_bw(base_size = 20) +
theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
panel.background = element_blank(), axis.line = element_line(colour = "black"), 
axis.text.x = element_text(angle = 45, hjust = 1)) +
labs(x = "Trait", y = expression(italic(h^2))) +
scale_y_continuous(limits = c(0, 0.6), expand = c(0,0))

heritabilities_blood_paper

ggsave("/local1/home/pazweifel/plots//blood/fsc_ssc/heritabilities.pdf", heritabilities_blood_paper, device = pdf, height = 10, width = 10)

write_csv(heritabilities, "/local1/home/pazweifel/plots/blood/heritabilities.csv")

sldsc_covstruct <- readRDS("S_LDSC_blood_factors_fsc_ssc.rds")


#specify the model syntax for a common factor model
model<-'F1 =~NA*RBCSSC+RBCFSC+RETFSC+IRFFSC
F2 =~NA*NEFSC+EOSSC+EOFSC+BASOFSC
F3 =~NA*PLTSSC+PLTFSC
F4 =~NA*EOSSC+EOFSC+BASOFSC+LYSSC+LYFSC
F1~~F2
F1~~F3
F1~~F4
F2~~F3
F2~~F4
F3~~F4
EOFSC ~~ a*EOFSC
a > 0.0001'

params<-c("F1~~F1","F2~~F2","F3~~F3","F4~~F4")

#use unit variance identification
std.lv=TRUE

#estimate enrichment using the enrich function
enrich_gwasbysub<-enrich(s_covstruc=sldsc_covstruct,model=model,params=params,std.lv=std.lv, rm_flank = TRUE)

enrich_blood_tau <- enrich(s_covstruc=sldsc_covstruct,model=model,params=params,std.lv=std.lv, rm_flank = TRUE, tau = TRUE)

annotations <- list.files(path = "/local1/scratch/pazweifel/single_cell_1k1k_newannot/", pattern = "\\.annot\\.gz$")

tr_annotations <- sub("\\..*", "", annotations)
tr_annotations <- unique(tr_annotations)

length(tr_annotations)

annotations_old <- list.files(path = "/local1/scratch/pazweifel/single_cell_1k1k/", pattern = "\\.annot\\.gz$")

tr_annotations <- sub("\\..*", "", annotations_old)
tr_annotations <- unique(tr_annotations)
length(tr_annotations)

enrich_df <- bind_rows(
    lapply(enrich_gwasbysub[1:4], function(factor) filter(factor, Annotation %in% tr_annotations)), .id = "Factor")

enrich_df_tau <- bind_rows(
    lapply(enrich_blood_tau[1:4], function(factor) filter(factor, Annotation %in% tr_annotations)), .id = "Factor")

tr_annotations

cell_type_groups <- list("T" = c('CD4_positive_alpha_beta_cytotoxic_T_cell', 'CD4_positive_alpha_beta_T_cell', 'CD8_positive_alpha_beta_T_cell','double_negative_thymocyte', 'central_memory_CD4_positive_alpha_beta_T_cell', 'central_memory_CD8_positive_alpha_beta_T_cell', 'effector_memory_CD4_positive_alpha_beta_T_cell', 'effector_memory_CD8_positive_alpha_beta_T_cell', 'gamma_delta_T_cell',  'mucosal_invariant_T_cell', 'naive_thymus_derived_CD4_positive_alpha_beta_T_cell', 'naive_thymus_derived_CD8_positive_alpha_beta_T_cell', 'regulatory_T_cell'),
                         "B" = c('memory_B_cell', 'naive_B_cell', 'transitional_stage_B_cell', 'plasmablast'),
                         "Precursors/Unspecific" = c('hematopoietic_precursor_cell', 'innate_lymphoid_cell', 'peripheral_blood_mononuclear_cell'),
                         "Dendritic" = c('dendritic_cell', 'conventional_dendritic_cell', 'plasmacytoid_dendritic_cell'),
                         "NK" = c('CD16_negative_CD56_bright_natural_killer_cell_human', 'natural_killer_cell'),
                         "Monocytes" = c('CD14_low_CD16_positive_monocyte', 'CD14_positive_monocyte'),
                         "RBC/platelet" = c('erythrocyte', 'platelet'),
                         "Control" = c('Control'))
                      

lookup_groups <- unlist(lapply(names(cell_type_groups), function(g) {
  setNames(rep(g, length(cell_type_groups[[g]])),
           cell_type_groups[[g]])
}))

lookup <- c("1" = "RBC", "2" = "Granulocytes", "3" = "Platelets", "4" = "Lymphocytes/Granulocytes")

enrich_df$Factor <- lookup[enrich_df$Factor]
significance_level <- 0.05/length(unique(enrich_df$Annotation))

enrich_df_tau$Factor <- lookup[enrich_df_tau$Factor]
significance_level <- 0.05/length(unique(enrich_df_tau$Annotation))

factor_order <- c("RBC", "Granulocytes","Platelets","Lymphocytes/Granulocytes")

enrich_df <- enrich_df %>% 
mutate(Bonf = Enrichment_p_value*length(unique(Annotation)), sig = Bonf < 0.05, log10p = -log10(Enrichment_p_value), group = lookup_groups[Annotation], Factor = factor(Factor, levels = factor_order)) %>% 
arrange(Factor, group, Annotation) %>% 
mutate(group = factor(group, levels = unique(group)), Annotation = factor(Annotation, levels = unique(Annotation))) %>%
filter(group != "Control")

enrich_df_tau <- enrich_df_tau %>% 
mutate(Bonf = Enrichment_p_value*length(unique(Annotation)), sig = Bonf < 0.05, log10p = -log10(Enrichment_p_value), group = lookup_groups[Annotation], Factor = factor(Factor, levels = factor_order)) %>% 
arrange(Factor, group, Annotation) %>% 
mutate(group = factor(group, levels = unique(group)), Annotation = factor(Annotation, levels = unique(Annotation))) %>%
filter(group != "Control")

group_colors <- c(
  "T" = "#0072B2",                   # blue
  "B" = "#E69F00",                   # orange
  "Precursors/Unspecific" = "#009E73",  # bluish green
  "Dendritic" = "#CC79A7",           # reddish purple
  "NK" = "#56B4E9",                  # sky blue
  "Monocytes" = "#D55E00",           # vermillion
  "RBC/platelet" = "#F0E442"         # yellow
)

enrich_df$Factor <- recode(enrich_df$Factor, Granulocytes = "F[Gran]", `Lymphocytes/Granulocytes` = "F[GranLy]", Platelets = "F[Plat]", RBC = "F[RBC]")
enrich_df_tau$Factor <- recode(enrich_df_tau$Factor, Granulocytes = "F[Gran]", `Lymphocytes/Granulocytes` = "F[GranLy]", Platelets = "F[Plat]", RBC = "F[RBC]")

enrich_plot <- ggplot(enrich_df, aes(x = Annotation, y = log10p)) +
geom_point(aes(color = group), size = 6, shape = 17) +
theme_bw(base_size = 25) +
geom_hline(yintercept = -log10(significance_level), color = "red", linetype = "dashed") +
scale_color_viridis_d(option = "D") +
geom_label_repel(data = filter(enrich_df, sig == TRUE), aes(label = Annotation)) +
labs(y = expression(-log[10](P)), color = "Group") +
facet_wrap(vars(Factor), nrow = 4, labeller = label_parsed) + 
theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
panel.background = element_blank(), axis.line = element_line(colour = "black"), axis.text.x = element_blank(), axis.ticks.x = element_blank(),
strip.background = element_blank(), strip.placement = "outside", axis.title.x = element_blank(), legend.position = "bottom") +
scale_y_continuous(limits = c(0, 7), expand = c(0, 0))



enrich_plot

ggsave(filename = "/local1/home/pazweifel/plots/blood/enrichment_4_factors.pdf", plot = enrich_plot, device = "pdf", width = 20, height = 16)

write_csv(enrich_df, "/local1/home/pazweifel/plots/blood/stratifiedgsem_zerocov_enrichment.csv")

enrich_plot_tau <- ggplot(enrich_df_tau, aes(x = Annotation, y = log10p)) +
geom_point(aes(color = group), size = 6, shape = 17) +
theme_bw(base_size = 25) +
geom_hline(yintercept = -log10(significance_level), color = "red", linetype = "dashed") +
scale_color_viridis_d(option = "D") +
geom_label_repel(data = filter(enrich_df_tau, sig == TRUE), aes(label = Annotation)) +
labs(y = expression(-log[10](P)), color = "Group") +
facet_wrap(vars(Factor), nrow = 4, labeller = label_parsed) + 
theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
panel.background = element_blank(), axis.line = element_line(colour = "black"), axis.text.x = element_blank(), axis.ticks.x = element_blank(),
strip.background = element_blank(), strip.placement = "outside", axis.title.x = element_blank(), legend.position = "bottom") +
scale_y_continuous(limits = c(0, 7), expand = c(0, 0))



enrich_plot_tau

ggsave(filename = "/local1/home/pazweifel/plots/blood/enrichment_4_factors_taucorrected.pdf", plot = enrich_plot_tau, device = "pdf", width = 20, height = 16)

write_csv(enrich_df_tau, "/local1/home/pazweifel/plots/blood/stratifiedgsem_tau_enrichment.csv")

blood_gwas <- readRDS("commonfactor/gwas_sumstats_factor4.R")

##Calculate Implied Sample Size for Factor
#restrict to MAF >= 10%
blood_gwas_N <- lapply(blood_gwas, function(gwas){
        subby<-subset(gwas, gwas$MAF >= .1)
        
        N_hat_F1<-round(mean(1/((2*subby$MAF*(1-subby$MAF))*subby$SE^2)), 0)
        mutate(tibble(gwas), N = N_hat_F1) %>% rename(P = Pval_Estimate, Z = Z_Estimate, Beta = est) %>% select(SNP, CHR, BP, MAF, A1, A2, Beta, SE, Z, P, N)
    })

dir.create("/local1/scratch/pazweifel/sumstats_blood/factor_sumstats", showWarnings = FALSE)
basepath <- "/local1/scratch/pazweifel/sumstats_blood/factor_sumstats/"
factors <- c("RBC_fsc_ssc_4.tsv", "Gran_fsc_ssc_4.tsv", "Plat_fsc_ssc_4.tsv", "GraLy_fsc_ssc_4.tsv")

blood_filepaths <- paste0(basepath, factors)

for (i in c(1:length(blood_gwas_N))){
    if (!file.exists(blood_filepaths[i])){
        write_tsv(blood_gwas_N[[i]], file = blood_filepaths[i])
        }}

cell_type_groups <- list("T" = c('CD4_positive_alpha_beta_cytotoxic_T_cell', 'CD4_positive_alpha_beta_T_cell', 'CD8_positive_alpha_beta_T_cell','double_negative_thymocyte', 'central_memory_CD4_positive_alpha_beta_T_cell', 'central_memory_CD8_positive_alpha_beta_T_cell', 'effector_memory_CD4_positive_alpha_beta_T_cell', 'effector_memory_CD8_positive_alpha_beta_T_cell', 'gamma_delta_T_cell',  'mucosal_invariant_T_cell', 'naive_thymus_derived_CD4_positive_alpha_beta_T_cell', 'naive_thymus_derived_CD8_positive_alpha_beta_T_cell', 'regulatory_T_cell'),
                         "B" = c('memory_B_cell', 'naive_B_cell', 'transitional_stage_B_cell', 'plasmablast'),
                         "Precursors/Unspecific" = c('hematopoietic_precursor_cell', 'innate_lymphoid_cell', 'peripheral_blood_mononuclear_cell'),
                         "Dendritic" = c('dendritic_cell', 'conventional_dendritic_cell', 'plasmacytoid_dendritic_cell'),
                         "NK" = c('CD16_negative_CD56_bright_natural_killer_cell_human', 'natural_killer_cell'),
                         "Monocytes" = c('CD14_low_CD16_positive_monocyte', 'CD14_positive_monocyte'),
                         "RBC/platelet" = c('erythrocyte', 'platelet'),
                         "Control" = c('Control'))

annotation_order_data_labels <- c(
    "B",
    "Dendritic",
    "Monocytes",
    "NK",
    "Precursors/Unspecific",
    "RBC/platelet",
    "T")

lookup_groups <- unlist(lapply(names(cell_type_groups), function(g) {
  setNames(rep(g, length(cell_type_groups[[g]])),
           cell_type_groups[[g]])
}))

#load the data
samples <- c("RBC", "Granulocytes", "Platelets", "Granulocytes/Lymphocytes")
names(samples) <- samples

pathlist <- c("RBC_fsc_ssc_4.h2.cell_type_results.txt", "Gran_fsc_ssc_4.h2.cell_type_results.txt", "Plat_fsc_ssc_4.h2.cell_type_results.txt", "GraLy_fsc_ssc_4.h2.cell_type_results.txt")
names(pathlist) <- samples

enrich_table_common <- bind_rows(lapply(samples, function(sample){
    a <- read_table(pathlist[sample])
    a %>% mutate(Phenotype = sample, Bonfer = Coefficient_P_value*nrow(a), sign = Bonfer < 0.05, BH = p.adjust(Coefficient_P_value, method = "fdr"), Supercluster= factor(lookup_groups[Name], levels = annotation_order_data_labels)) %>% arrange(Supercluster, Name)
    })) %>% mutate(Phenotype = factor(Phenotype, levels = unique(Phenotype)), Name = factor(Name, levels = unique(Name))) 

head(enrich_table_common)

p_cut <- max(enrich_table_common$Coefficient_P_value[enrich_table_common$BH <= 0.05], na.rm = TRUE)

enrich_table_common$Phenotype <- recode(enrich_table_common$Phenotype, RBC = "F[RBC]", Granulocytes = "F[Gran]", Platelets = "F[Plat]", `Granulocytes/Lymphocytes` = "F[GranLy]")

sldsc_plot <- ggplot(enrich_table_common, aes(x = Name, y = -log10(Coefficient_P_value))) +
geom_point(aes(color = Supercluster), size = 6, shape = 17) +
scale_color_viridis_d(option = "D") +
geom_label_repel(data = filter(enrich_table_common, sign == TRUE), aes(label = Name)) +
geom_hline(yintercept = -log10(0.05/length(unique(enrich_table_common$Name))), color = "red", linetype = "dashed") +
labs(y = expression(-log[10](P)), color = "Category") + 
facet_wrap(vars(Phenotype), nrow = 5, labeller = label_parsed) +
theme_bw(base_size = 25) +
theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
panel.background = element_blank(), axis.line = element_line(colour = "black"), axis.text.x = element_blank(), axis.ticks.x = element_blank(),
strip.background = element_blank(), strip.placement = "outside", axis.title.x = element_blank(), legend.position = "bottom") +
scale_y_continuous(limits = c(0, 7), expand = c(0, 0))

sldsc_plot

ggsave("/local1/home/pazweifel/plots/blood/sldsc_4_factors.pdf", plot = sldsc_plot, device = "pdf", width = 20, height = 16)

write_csv(enrich_table_common, "/local1/home/pazweifel/plots/blood/sldsc_enrichment.csv")

sldsc_results_paths <- list.files("sldsc_original_sumstats", pattern = "cell_type_results\\.txt$", full.names = TRUE)

sldsc_results_names <- list.files("sldsc_original_sumstats", pattern = "cell_type_results\\.txt$", full.names = FALSE)
sldsc_results_names <- str_remove(sldsc_results_names, pattern = ".h2.cell_type_results.txt")

names(sldsc_results_names) <- sldsc_results_names
names(sldsc_results_paths) <- sldsc_results_names

enrich_table_original_sumstats <- bind_rows(lapply(sldsc_results_names, function(sample){
    a <- read_table(sldsc_results_paths[sample])
    a %>% mutate(Phenotype = sample, Bonfer = Coefficient_P_value*nrow(a), sign = Bonfer < 0.05, BH = p.adjust(Coefficient_P_value, method = "fdr"), Supercluster= factor(lookup_groups[Name], levels = annotation_order_data_labels)) %>% arrange(Supercluster, Name)
    })) %>% mutate(Phenotype = factor(Phenotype, levels = unique(Phenotype)), Name = factor(Name, levels = unique(Name))) 

p_cut <- max(enrich_table_original_sumstats$Coefficient_P_value[enrich_table_original_sumstats$BH <= 0.05], na.rm = TRUE)

#enrich_table_original_sumstats$Phenotype <- recode(enrich_table_original_sumstats$Phenotype, RBC = "F[RBC]", Granulocytes = "F[Gran]", Platelets = "F[Plat]", `Granulocytes/Lymphocytes` = "F[GranLy]")

head(enrich_table_original_sumstats)

sldsc_plot_originals <- ggplot(enrich_table_original_sumstats, aes(x = Name, y = -log10(Coefficient_P_value))) +
geom_point(aes(color = Supercluster), size = 4, shape = 17) +
scale_color_viridis_d(option = "D") +
geom_label_repel(data = filter(enrich_table_original_sumstats, sign == TRUE), aes(label = Name)) +
geom_hline(yintercept = -log10(0.05/length(unique(enrich_table_original_sumstats$Name))), color = "red", linetype = "dashed") +
labs(y = expression(-log[10](P)), color = "Category") + 
facet_wrap(vars(Phenotype), nrow = 8, ncol = 8) +
theme_bw(base_size = 25) +
theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
panel.background = element_blank(), axis.line = element_line(colour = "black"), axis.text.x = element_blank(), axis.ticks.x = element_blank(),
strip.background = element_blank(), strip.placement = "outside", axis.title.x = element_blank()) 

sldsc_plot_originals

ggsave("/local1/home/pazweifel/plots/blood/sldsc_original_datasets.pdf", plot = sldsc_plot_originals, device = "pdf", width = 40, height = 40)

write_csv(enrich_table_original_sumstats, "/local1/home/pazweifel/plots/blood/sldsc_enrichment_original_sumstats.csv")

reticulocyte_measures <- c("LFR","MFR","RPI")

sign_table <- filter(enrich_table_original_sumstats, sign == TRUE) %>% mutate(cell_type = factor(ifelse(Phenotype %in% reticulocyte_measures, "Reticulocyte", "Platelet"), levels = c("Platelet", "Reticulocyte"))) %>% arrange(cell_type, Bonfer) %>% mutate(Phenotype = factor(Phenotype, levels = unique(Phenotype)))

sign_table

group_colors <- c(
    "Platelet" = "brown",
    "Reticulocyte" = "yellow")

sign_plot <- ggplot(sign_table, aes(x = Phenotype, y = -log10(Coefficient_P_value))) +
geom_point(aes(color = cell_type, shape = Name), size = 5) +
#geom_label_repel(aes(label = Name), ) +
theme_bw(base_size = 25) +
labs(color = "Cell type", shape = "Annotation", y = expression(-log[10](P)), x = "Phenotype") +
theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
panel.background = element_blank(), axis.line = element_line(colour = "black"), axis.text.x = element_text(angle = 45, hjust = 1)) +
scale_color_manual(values = group_colors) 

sign_plot

ggsave("/local1/home/pazweifel/plots/blood/sign_enrichments_original_sumstats.pdf", device = "pdf", plot = sign_plot, width = 18, height = 10)

knitr::kable(metadata)

full_covstruct <- readRDS("matrices/all.R")

sstand <- full_covstruct$S_Stand
vstand <- full_covstruct$V_Stand

#extract the standard error
vstand_diag <- sapply(c(1:dim(vstand)[1]), function(index) vstand[index, index])
                      
#copy the structure of the sstand matrix
error_matrix <- sstand
                      
#overwrite lower triangle with vstand values
error_matrix[lower.tri(error_matrix, diag = TRUE)] <- vstand_diag
                      
#transpose the matrix and copy the values to the upper triangle, making it symmetric
error_matrix[upper.tri(error_matrix)] <- t(error_matrix)[upper.tri(error_matrix)]
error_matrix <- sqrt(error_matrix)
error_matrix
                      
#vectorize the error matrix so that you can use it for the ggplot
error_vectorized <- c(error_matrix)

#initialize the dataframe for the matrix values of both matrices
tempnames <- dimnames(sstand)[[2]]
expanded_df <- expand.grid(x = tempnames, y = tempnames)

#vectorize the standard matrix and store matrix values in df columns
sstand_vectorized <- c(sstand)
expanded_df$values <- round(sstand_vectorized, 2)
expanded_df$se <- round(error_vectorized, 2)
expanded_df

write_csv(expanded_df, "/local1/home/pazweifel/plots/blood/all_correlations.csv")

metadata_expanded <- metadata %>%
  mutate(
    category = case_when(
      # --- MEP: megakaryocyte–erythroid progenitors (platelets / erythrocytes / reticulocytes)
      str_detect(ABBRV, "^(PLT|HIPF|IPF|RBC|RET|HFR|MFR|LFR|IRF|RDW|MacroR|MicroR|HyperHe|DeltaHe|DeltaHGB|RPI)") ~ "MEP",

      # --- GMP: granulocyte–monocyte progenitors (granulocytes + monocytes + dendritic)
      str_detect(ABBRV, "^(NE|IG|EO|BASO|MO)") |
        str_detect(Trait, regex("granulocyte|eosinophil|basophil|neutrophil|monocyte|dendritic", ignore_case = TRUE)) ~ "GMP",

      # --- CLP: common lymphoid progenitors (B/T/NK/lymphocyte)
      str_detect(ABBRV, "^(LY|RELYMP)") |
        str_detect(Trait, regex("lymphocyte|B cell|T cell|natural killer|NK|dendritic", ignore_case = TRUE)) ~ "CLP",

      TRUE ~ "unspecific"
    )
  )

metadata_expanded[metadata_expanded["Trait"] == "Platelet large cell ratio", ]$category <- "MEP"

expanded_df_annotated <- expanded_df %>% left_join(metadata_expanded, by = join_by("x" == "ABBRV"))
head(expanded_df_annotated)

# 1) set your desired category order (optional but recommended)
expanded_df_annotated$category <- factor(expanded_df_annotated$category, levels = c("MEP", "GMP", "CLP"))

# 2) build a master order for phenotypes (based on category grouping)
ord <- expanded_df_annotated %>%
  distinct(x, category) %>%
  arrange(category, x) %>%
  pull(x)

# 3) apply the same order to BOTH axes
expanded_df_annotated <- expanded_df_annotated %>%
  mutate(
    x = factor(x, levels = ord),
    y = factor(y, levels = ord)   # <-- key line for the diagonal
  )

head(expanded_df_annotated)

#create the ggplot
correlation_plot <- ggplot(expanded_df_annotated, aes(x = x, y = y, fill = values)) +
scale_fill_gradient(low = "darkorchid", high = "orange", limits = c(0, max(expanded_df$values))) +
theme_minimal() +
geom_tile() +
theme(axis.text = element_text(face = "bold", size = 12), axis.title = element_blank(), legend.text = element_text(size = 12), legend.title = element_text(size = 12), axis.text.x = element_text(angle = 45, hjust = 1)) +
labs(fill = expression(r[g])) +
geom_text(aes(label = values)) +
geom_text(aes(label = paste0("(",se,")"), vjust = 2))

correlation_plot

ggsave("/local1/home/pazweifel/plots/blood/correlations_all.pdf", plot = correlation_plot, width = 40, height = 40, device = "pdf")

expanded_df_annotated %>% group_by(category, x) %>% summarise(mean_corr = mean(values)) %>% arrange(category, desc(mean_corr)) %>% slice_head(n=3) %>% ungroup()

library(tibble)
mat <- expanded_df_annotated %>%
  select(x, y, values) %>%
  pivot_wider(names_from = y, values_from = values) %>%
  column_to_rownames("x") %>%
  as.matrix()

hc <- hclust(as.dist(1 - mat), method = "ward.D2")

ord <- hc$order
mat <- mat[ord, ord]

pdf("/local1/home/pazweifel/plots/blood/heatmap_all.pdf")
heatmap(mat,
        labRow = rownames(mat),
        labCol = colnames(mat),
        cexRow = 0.6,     # shrink row labels
        cexCol = 0.6) 
dev.off()

list.files("sldsc_individual_annotation_covstructs/")


model<-'F1 =~NA*RBCSSC+RBCFSC+RETFSC+IRFFSC
F2 =~NA*NEFSC+EOSSC+EOFSC+BASOFSC
F3 =~NA*PLTSSC+PLTFSC
F4 =~NA*EOSSC+EOFSC+BASOFSC+LYSSC+LYFSC
F1~~F2
F1~~F3
F1~~F4
F2~~F3
F2~~F4
F3~~F4
EOFSC ~~ a*EOFSC
a > 0.0001'

params<-c("F1~~F1","F2~~F2","F3~~F3","F4~~F4")

#use unit variance identification
std.lv=TRUE

#estimate enrichment using the enrich function
enrich_gwasbysub<-enrich(s_covstruc=sldsc_covstruct,model=model,params=params,std.lv=std.lv, rm_flank = TRUE)

library(dplyr)
library(stringr)
library(purrr)

covstruct_dir <- "sldsc_individual_annotation_covstructs"

covstruct_files <- list.files(
  path = covstruct_dir,
  pattern = "\\.rds$",
  full.names = TRUE
)

annotation_names <- str_extract(
  basename(covstruct_files),
  "(?<=ssc_).*(?=\\.rds$)"
)

model <- '
F1 =~ NA*RBCSSC + RBCFSC + RETFSC + IRFFSC
F2 =~ NA*NEFSC + EOSSC + EOFSC + BASOFSC
F3 =~ NA*PLTSSC + PLTFSC
F4 =~ NA*EOSSC + EOFSC + BASOFSC + LYSSC + LYFSC
F1~~F2
F1~~F3
F1~~F4
F2~~F3
F2~~F4
F3~~F4
EOFSC ~~ a*EOFSC
a > 0.0001
'

params <- c("F1~~F1", "F2~~F2", "F3~~F3", "F4~~F4")
std.lv <- TRUE

# turn full annotation path into plain annotation name
clean_annotation <- function(x) {
  x <- as.character(x)
  x <- str_replace(x, "/+$", "")   # remove trailing slash(es)
  x <- ifelse(str_detect(x, "/"), basename(x), x)
  x
}

all_enrichments <- map2_dfr(covstruct_files, annotation_names, function(file_path, annotation_name) {
  
  message("Processing: ", annotation_name)
  
  sldsc_covstruct <- readRDS(file_path)
  
  enrich_out <- enrich(
    s_covstruc = sldsc_covstruct,
    model = model,
    params = params,
    std.lv = std.lv,
    rm_flank = TRUE,
    tau = TRUE
  )
  
  enrich_df <- bind_rows(enrich_out, .id = "Factor") %>%
    mutate(
      Factor = as.integer(Factor),
      Annotation_clean = clean_annotation(Annotation)
    )
  
  out <- enrich_df %>%
    filter(Annotation_clean == annotation_name) %>%
    mutate(Source_annotation = annotation_name)
  
  if (nrow(out) == 0) {
    message("No rows matched for: ", annotation_name)
    message("Available cleaned annotations:")
    print(unique(enrich_df$Annotation_clean))
  }
  
  out
})

write.table(
  all_enrichments,
  file = "combined_factor_annotation_enrichments.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

library(dplyr)
library(stringr)

clean_annotation <- function(x) {
  x <- as.character(x)
  x <- str_replace(x, "/+$", "")   # remove trailing slashes
  x <- ifelse(str_detect(x, "/"), basename(x), x)
  x
}

all_enrichments <- all_enrichments %>%
  mutate(Annotation = clean_annotation(Annotation))

lookup <- c("1" = "RBC", "2" = "Granulocytes", "3" = "Platelets", "4" = "Lymphocytes/Granulocytes")

all_enrichments$Factor <- lookup[all_enrichments$Factor]
significance_level <- 0.05/length(unique(all_enrichments$Annotation))


factor_order <- c("RBC", "Granulocytes","Platelets","Lymphocytes/Granulocytes")

all_enrichments <- all_enrichments %>% 
mutate(Bonf = Enrichment_p_value*length(unique(Annotation)), sig = Bonf < 0.05, log10p = -log10(Enrichment_p_value), group = lookup_groups[Annotation], Factor = factor(Factor, levels = factor_order)) %>% 
arrange(Factor, group, Annotation) %>% 
mutate(group = factor(group, levels = unique(group)), Annotation = factor(Annotation, levels = unique(Annotation))) %>%
filter(group != "Control")

group_colors <- c(
  "T" = "#0072B2",                   # blue
  "B" = "#E69F00",                   # orange
  "Precursors/Unspecific" = "#009E73",  # bluish green
  "Dendritic" = "#CC79A7",           # reddish purple
  "NK" = "#56B4E9",                  # sky blue
  "Monocytes" = "#D55E00",           # vermillion
  "RBC/platelet" = "#F0E442"         # yellow
)

all_enrichments$Factor <- recode(all_enrichments$Factor, Granulocytes = "F[Gran]", `Lymphocytes/Granulocytes` = "F[GranLy]", Platelets = "F[Plat]", RBC = "F[RBC]")

enrich_plot_oneatatime <- ggplot(all_enrichments, aes(x = Annotation, y = log10p)) +
geom_point(aes(color = group), size = 6, shape = 17) +
theme_bw(base_size = 25) +
geom_hline(yintercept = -log10(significance_level), color = "red", linetype = "dashed") +
scale_color_viridis_d(option = "D") +
geom_label_repel(data = filter(all_enrichments, sig == TRUE), aes(label = Annotation)) +
labs(y = expression(-log[10](P)), color = "Group") +
facet_wrap(vars(Factor), nrow = 4, labeller = label_parsed) + 
theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
panel.background = element_blank(), axis.line = element_line(colour = "black"), axis.text.x = element_blank(), axis.ticks.x = element_blank(),
strip.background = element_blank(), strip.placement = "outside", axis.title.x = element_blank(), legend.position = "bottom") +
scale_y_continuous(limits = c(0, 7), expand = c(0, 0))



enrich_plot_oneatatime

ggsave(filename = "/local1/home/pazweifel/plots/blood/enrichment_4_factors_ldsc_onatatime_nocontrol.pdf", plot = enrich_plot_oneatatime, device = "pdf", width = 20, height = 16)

write_csv(all_enrichments, "/local1/home/pazweifel/plots/blood/stratifiedgsem_oneatatime.csv")

###figure 2: plots --> enrich_plot, sldsc_plot

library(cowplot)

common_legend <- get_legend(enrich_plot)

enrich_plot_wol <- enrich_plot + theme(legend.position = "none")

sldsc_plot_wol <- sldsc_plot + theme(legend.position = "none")

combined_plot_fig2 <- plot_grid(plotlist = list(enrich_plot_wol, sldsc_plot_wol), ncol = 2)
combined_plot_fig2_complete <- plot_grid(
    plotlist = list(combined_plot_fig2, common_legend),
    nrow = 2,
    rel_heights = c(1, 0.1)
    )

combined_plot_fig2_complete

ggsave("/local1/home/pazweifel/plots/blood/figure2.pdf", plot = combined_plot_fig2_complete, device = "pdf", width = 14, height = 14)

### Figure 3: plots --> enrich_plot_tau, enrich_plot_oneatatime

common_legend <- get_legend(enrich_plot_tau)

enrich_plot_tau_wol <- enrich_plot_tau + theme(legend.position = "none")

enrich_plot_oneatatime_wol <- enrich_plot_oneatatime + theme(legend.position = "none")

combined_plot_fig3 <- plot_grid(plotlist = list(enrich_plot_tau_wol, enrich_plot_oneatatime_wol), ncol = 2)
combined_plot_fig3_complete <- plot_grid(
    plotlist = list(combined_plot_fig3, common_legend),
    nrow = 2,
    rel_heights = c(1, 0.1)
    )

combined_plot_fig3_complete

ggsave("/local1/home/pazweifel/plots/blood/figure3.pdf", plot = combined_plot_fig3_complete, device = "pdf", width = 14, height = 14)


