library(dplyr)
library(readr)
library(tidyr)
library(stringr)
library(ggplot2)
library(microshades)
library(openxlsx)
library(GenomicSEM)



sample_sizes <- read_csv(
  "
Pheno,Nca,Nco,Neff
MDD,357636,1281936,1118502.79
ADHD,38691,186843,128213.795
AN,16992,55525,46321.9
ALC,NA,NA,26853.43
ASD,18381,27969,44366.62
BIP,41917,371549,150669.892
BMI,NA,NA,681275
CAD,95830,385536,285937.77
Chrono,NA,NA,449734
CRP,NA,NA,204402
Height,NA,NA,693529
Neu,NA,NA,390278
PD,2147,7760,6665.06
PTSD,23212,151447,140475
SCZ,53386,77258,126281.975077309
Smoking,NA,NA,263954
T2D,80154,853816,251739.50
"
)

#sample prevalences can be set to 0.5 as N for the munged files has been converted to effective N. It could also be odne by specifying the correct sample prevalence when convertint to the liability scale, but for compatibility purposes the neffe variant has been applied
mental_prevalences <- c("MDD" = 0.08, "ADHD" = 0.05, "ASD" = 0.012, "BIP" = 0.02, "SCZ" = 0.005)
mental_names <- names(mental_prevalences)

other_prevalences <- c(
  "ALC" = 0.16, "AN" = 0.01, "BMI" = NA,
  "CAD" = 0.05, "Chrono" = NA, "CRP" = NA,
  "Height" = NA, "Neu" = NA, "PD" = 0.03,
  "PTSD" = 0.3, "Smoking" = NA,
  "T2D" = 0.06
)
other_names <- names(other_prevalences)

all_prevalences <- c(mental_prevalences, other_prevalences)
all_names <- names(all_prevalences)

#load the paths of necessary documents as not all of them are in this directory
ld <- "/local1/hdata/REF/eur_ref_ld_chr"
weights <- "/local1/hdata/REF/eur_w_ld_chr"
samples <- c("/local1/scratch/pazweifel/jolien_paper_sumstats/munged_sumstats/mdd_munged.sumstats.gz",
             "/local1/scratch/pazweifel/jolien_paper_sumstats/munged_sumstats/adhd_munged.sumstats.gz",
             "/local1/scratch/pazweifel/jolien_paper_sumstats/munged_sumstats/asd_munged.sumstats.gz",
             "/local1/scratch/pazweifel/jolien_paper_sumstats/munged_sumstats/bip_munged.sumstats.gz",
            "/local1/scratch/pazweifel/jolien_paper_sumstats/munged_sumstats/scz_munged.sumstats.gz")
other_basepath <- "/local1/hdata/sumstats/munged/" #phenotype.sumstats.gz

#create a matrices directory
#dir.create("matrices_jolien", showWarnings = FALSE)
#create a path for a covstruct object where you store the covariance structure of the bipolar disorders and store it in the matrices directory
covstruct_mental <- here::here("matrices_jolien",paste("mental", "R", sep = "."))
#extract population prevalence 
pop_prev <- as.vector(all_prevalences[mental_names])
#fill the samp_prev vector depending on the pop_prev vector, if NA use NA and if not use 0.5
samp_prev <- as.vector(ifelse(is.na(pop_prev), NA, 0.5))
#check if the covstruct file already exists to not do the analysis unnecessarily, first time no file there
if (!file.exists(covstruct_mental)){
    covstruct_mental_ldsc <- ldsc(
        traits = samples,
        trait.names = mental_names,
        sample.prev = samp_prev,
        population.prev = pop_prev,
        ld = ld,
        wld = weights,
        stand = TRUE
    )
    #store the file in the previously defined path, now the path is, contorl tell how to deparse
    dput(covstruct_mental_ldsc, covstruct_mental, control = c("all", "digits17"))
}

samples_test <- samples <- c("/local1/scratch/pazweifel/jolien_paper_sumstats/munged_sumstats/mdd_munged_ndiv6.sumstats.gz",
             "/local1/scratch/pazweifel/jolien_paper_sumstats/munged_sumstats/adhd_munged.sumstats.gz",
             "/local1/scratch/pazweifel/jolien_paper_sumstats/munged_sumstats/asd_munged.sumstats.gz",
             "/local1/scratch/pazweifel/jolien_paper_sumstats/munged_sumstats/bip_munged.sumstats.gz",
            "/local1/scratch/pazweifel/jolien_paper_sumstats/munged_sumstats/scz_munged.sumstats.gz")

covstruct_mental <- here::here("matrices_jolien",paste("mental_ndiv6", "R", sep = "."))

if (!file.exists(covstruct_mental)){
    covstruct_mental_ldsc <- ldsc(
        traits = samples_test,
        trait.names = mental_names,
        sample.prev = samp_prev,
        population.prev = pop_prev,
        ld = ld,
        wld = weights,
        stand = TRUE
    )
    #store the file in the previously defined path, now the path is, contorl tell how to deparse
    dput(covstruct_mental_ldsc, covstruct_mental, control = c("all", "digits17"))
}

#First define the base model
#get the base covariance matrix from the matrices directory
covstruct_base_test <- dget(here::here("matrices_jolien/mental_ndiv6.R"))
#define the model using NA*Clinical to say that this parameter loading has to be freely estimated, we fix the factor variance to 1 and ensure residual genetic variance of phenotypes isnt negative
base_model <- "F1=~NA*MDD+ASD+ADHD+BIP+SCZ
F1~~1*F1
MDD~~a*MDD
a > 0.0001
"

base.fit_test <- usermodel(covstruct_base_test,
                      estimation = "DWLS",
                      model = base_model,
                      imp_cov = TRUE
)

#se <- base.fit_test$results$STD_Genotype_SE
#est <- base.fit_test$results$STD_Genotype
#est.se <- paste0(round(as.numeric(est), 2),"\n (",round(as.numeric(se), 2),")")

#fit_sem_base_test <- semPlotModel_GSEM(base.fit_test)

#semPaths(fit_sem_base_test, whatLabels = "std", layout = "tree", edge.color = "black", sizeMan = 12, sizeLat = 12, nCharNodes = 0, residuals = TRUE, style = "mx",
   #     edge.label.cex = 1.2,
  #      esize = 2,
  #      edgeLabels = est.se,
  #      edge.label.position = 0.5
  #      )


original_sumstats <- c("/local1/scratch/pazweifel/jolien_paper_sumstats/cleaned_sumstats/mdd_eur_neff_2.gz",
                       "/local1/scratch/pazweifel/jolien_paper_sumstats/cleaned_sumstats/adhd_eur_neff.gz",
                       "/local1/scratch/pazweifel/jolien_paper_sumstats/cleaned_sumstats/asd_eur_neff.gz",
                       "/local1/scratch/pazweifel/jolien_paper_sumstats/cleaned_sumstats/bip_eur_neff_2.gz",
                       "/local1/scratch/pazweifel/jolien_paper_sumstats/cleaned_sumstats/scz_eur_neff_2.gz")

ref <- "/local1/scratch/pazweifel/jolien_paper_sumstats/reference.1000G.maf.0.005.txt"
se.logit <- c(T,T,T,T,T)

mental_sumstats_path <- here::here("sumstats_jolien",paste("mental_sumstats","R", sep = "."))

if (!file.exists(mental_sumstats_path)){
    mental_sumstats <- sumstats(files = original_sumstats, trait.names = mental_names, ref = ref, se.logit = se.logit)
    dput(mental_sumstats, mental_sumstats_path, control = c("all", "digits17"))
}

#load the required files
#######covstruct_mental <- dget(here::here("matrices_jolien/mental.R"))
#mental_sumstats <- dget(here::here("sumstats_jolien/mental_sumstats.R"))
#create the path to the object you want the output to store
#####mental_common_factor_path <- here::here("common_factor_jolien","mental_common_factor.R")
####if (!file.exists(mental_common_factor_path)){
   ### mental_common_factor_gwas <- commonfactorGWAS(covstruc = covstruct_mental, SNPs = mental_sumstats, parallel = TRUE, cores = 50)
    ##dput(mental_common_factor_gwas, mental_common_factor_path, control = c("all", "digits17"))
#}

# load the resulsts to investigate them
#commonfactorgwas_jolien <- dget(here::here("common_factor_jolien/mental_common_factor.R"))

#library(qqman)

#head(commonfactorgwas_jolien)

#manhattan(filter(commonfactorgwas_jolien, warning == "0"), chr = "CHR", bp = "BP", p = "Pval_Estimate", snp = "SNP", suggestiveline = F)

#commonfactorgwas_jolien_filtered <- filter(commonfactorgwas_jolien, warning == "0")

#write_delim(commonfactorgwas_jolien_filtered, here::here("common_factor_jolien/mental_common_factor_filtered.tsv"))

##library(R.utils)
#gzip(here::here("common_factor_jolien/mental_common_factor_filtered.tsv"), remove = FALSE)

#filter(commonfactorgwas_jolien, Pval_Estimate < 1e-150)

#ggplot(commonfactorgwas_jolien, aes(x = as.factor(CHR), y = se_c)) + geom_boxplot(color = "darkorchid") + theme_bw() + scale_y_log10()

#ggplot(commonfactorgwas_jolien, aes(x = as.factor(CHR), y = est)) + geom_boxplot(color = "darkorchid") + theme_bw()

####summary(log10(commonfactorgwas_jolien$se_c))
###range(co####monfactorgwas_jolien$se_c, na.rm = TRUE)
##sum(commonfactorgwas_jolien$se_c < 1e-10, na.rm = TRUE)
#sum(!is.finite(commonfactorgwas_jolien$se_c))


###head(commonfactorgwas_jolien)
##commonfactorgwas_jolien_filter_est <- filter(commonfactorgwas_jolien, abs(est) < 0.05 & Pval_Estimate != 0)
#nrow(commonfactorgwas_jolien)
#nrow(commonfactorgwas_jolien_filter_est)
#commonfactorgwas_jolien_filter_est_se <- filter(commonfactorgwas_jolien_filter_est, abs(se_c) < 0.2 & abs(se_c) > 1e-5)
#nrow(commonfactorgwas_jolien_filter_est)
#nrow(commonfactorgwas_jolien_filter_est_se)

#ggplot(commonfactorgwas_jolien_filter_est_se, aes(x=as.factor(CHR), y = -log10(Pval_Estimate))) + geom_boxplot() + theme_bw()

#qqnorm(-log10(commonfactorgwas_jolien_filter_est_se$Pval_Estimate))
#qqline(-log10(commonfactorgwas_jolien_filter_est_se$Pval_Estimate))

#
#qq(filter(commonfactorgwas_jolien, Pval_Estimate != 0)$Pval_Estimate)

#head(mental_sumstats)

####extract the snp id column from the filtered dataset
##filtered_snps <- commonfactorgwas_jolien_filter_est_se$SNP
#head(filtered_snps)
##sort the sumstats file with the extracted snpvector
#mental_sumstats_filtered2 <- filter(mental_sumstats, SNP %in% filtered_snps)

##nrow(mental_sumstats_filtered2)
#nrow(mental_sumstats)

##save the filtered sumstats file in the correct folder
#dput(mental_sumstats_filtered2, here::here("sumstats_jolien","mental_sumstats_filtered2.R"), control = c("all", "digits17"))

######checking the distribution of pvalues in the entry sumstats
#load all the sumstats
#####sumstats_list <- list(0,0,0,0,0)
####i <- 1
###names(sumstats_list) <- mental_names
##mental_list <- mental_names
#names(mental_list) <- mental_names

#####
#####for (path in original_sumstats){
####    temp <- read.table(path, fill = TRUE, blank.lines.skip = TRUE, header = TRUE)
###    sumstats_list[[i]] <- temp
##    i <- i + 1
#}


###clean them
##sumstats_cleaned <- lapply(mental_list, function(trait){
##    filter(sumstats_list[[trait]], !is.na(N_EFF))
##    }
##                           )
    

##plot them
#library(qqman)

#for (trait in mental_list){
#    if (trait %in% c("MDD", "SCZ")){
#        print(trait)
#        qq(sumstats_cleaned[[trait]]$PVAL)
#    }
#    else {
#        print(trait)
#        qq(sumstats_cleaned[[trait]]$P)
#        }
#}

# do the same thing but include the other phenotypes one by one and name them after the current phenotype used
for (phenotype in other_names){
    covstruct_pheno <- here::here("matrices_jolien",paste(phenotype,"R", sep = "."))
    #extract population prevalences
    pop_prev <- c(all_prevalences[mental_names], all_prevalences[phenotype])
    #fill the sample prevalences depending on pop vector
    samp_prev <- as.vector(ifelse(is.na(pop_prev), NA, 0.5))
    #concatenate the sample vector with the bip filepaths and the specific single phenotype filepath
    samples_plus1 <- c(samples, paste0(other_basepath,phenotype,".sumstats.gz"))
    #check if file already exists
    if (!file.exists(covstruct_pheno)){
        covstruct_ldsc <- ldsc(
            traits = samples_plus1,
            trait.names = c(mental_names, phenotype),
            sample.prev = samp_prev,
            population.prev = pop_prev,
            ld = ld,
            wld = weights,
            stand = TRUE
        )
        dput(covstruct_ldsc, covstruct_pheno, control = c("all", "digits17"))
    }
}
        

#First define the base model
#get the base covariance matrix from the matrices directory
covstruct_base <- dget(here::here("matrices_jolien/mental.R"))
#define the model using NA*Clinical to say that this parameter loading has to be freely estimated, we fix the factor variance to 1 and ensure residual genetic variance of phenotypes isnt negative
base_model <- "F1=~NA*MDD+ASD+ADHD+BIP+SCZ
F1~~1*F1
MDD~~a*MDD
a > 0.0001
"

base.fit <- usermodel(covstruct_base,
                      estimation = "DWLS",
                      model = base_model,
                      imp_cov = TRUE
)

#same thing but with fixed loading of 1 on the clinical, factor variance freely estimated
anchor_model <- "F1=~1*MDD+ASD+ADHD+BIP+SCZ
F1~~F1
"
anchor.fit <- usermodel(covstruct_base,
                        estimation = "DWLS",
                        model = anchor_model,
                        imp_cov = TRUE
)

anchor.fit$results

sstand <- covstruct_base$S_Stand
vstand <- covstruct_base$V_Stand

#extract the standard error
vstand_diag <- sapply(c(1:15), function(index) vstand[index, index])
#copy the structure of the sstand matrix
error_matrix <- sstand
#overwrite lower triangle with vstand values
error_matrix[lower.tri(error_matrix, diag = TRUE)] <- vstand_diag
#transpose the matrix and copy the values to the upper triangle, making it symmetric
error_matrix[upper.tri(error_matrix)] <- t(error_matrix)[upper.tri(error_matrix)]
error_matrix <- sqrt(error_matrix)
error_matrix
#vectorize the matrix so that you can use it for the ggplot
error_vectorized <- c(error_matrix)
tempnames <- dimnames(sstand)[[2]]
expanded_df <- expand.grid(x = tempnames, y = tempnames)

sstand_vectorized <- c(sstand)
expanded_df$values <- round(sstand_vectorized, 2)
expanded_df$se <- round(error_vectorized, 2)
expanded_df
#create the ggplot
correlation_plot <- ggplot(expanded_df, aes(x = x, y = y, fill = values)) +
scale_fill_gradient(low = "darkorchid", high = "orange", limits = c(0, max(expanded_df$values))) +
theme_minimal() +
geom_tile() +
theme(axis.text = element_text(face = "bold", size = 12), axis.title = element_blank(), legend.text = element_text(size = 12), legend.title = element_text(size = 12)) +
labs(fill = expression(r[g])) +
geom_text(aes(label = values)) +
geom_text(aes(label = paste0("(",se,")"), vjust = 2))

correlation_plot

#define the common and the independent pathway models and introduce a placeholder which can filled by the repsective phenotype corresponding to the used covariance matrix
#set the factor variance to 1
common.model <- "F1=~NA*MDD+ASD+ADHD+BIP+SCZ
F1~~{trait}
F1~~1*F1
MDD~~a*MDD
a > 0.0001
"

independent.model <- "F1=~NA*MDD+ASD+ADHD+BIP+SCZ
MDD~~{trait}
ASD~~{trait}
ADHD~~{trait}
BIP~~{trait}
SCZ~~{trait}
F1~~1*F1
MDD~~a*MDD
a > 0.0001
"

#define a list with trait names named after themselves, so that that the naming remains when using this list with lapply and therefore specific elements can be accessed by name
traits <- other_names
names(traits) <- traits
traits

#make a named vector with all the common models and all the applied models, take advantage of the palceholder introduced
traits_common.models <- lapply(traits, 
                               function(trait) str_glue(common.model)
)

traits_independent.models <- lapply(traits,
                                    function(trait) str_glue(independent.model)
)

traits_independent.models

#load the genetic covariance matrices
trait.covstructs <- lapply(traits, 
                           function(trait) dget(here::here("matrices_jolien",paste(trait,"R", sep = ".")))
)

#now fit the common and the independent models using the two vectors
trait_common.fit <- lapply(traits, 
                           function(trait) {
                               usermodel(trait.covstructs[[trait]],
                               estimation = "DWLS",
                               model = traits_common.models[[trait]],
                               imp_cov = TRUE
    )
  }
)

trait_independent.fit <- lapply(traits,
                                function(trait) {
                                    usermodel(
                                        trait.covstructs[[trait]],
                                        estimation = "DWLS",
                                        model = traits_independent.models[[trait]],
                                        imp_cov = TRUE
    )
  }
)

#define the chisquared function: a smaller chi squared value means that the estimated models covariance structure fits the observed covariance structure well
#if you subtract the chi squared of the independent model from the common one, you can see how much better it is
#the difference follows a chisquared distribution and you can calculate the p value of this observation

model_diff_chisq <- function(common_model, independent_model) {
    q_chisq <- common_model$modelfit$chisq - independent_model$modelfit$chisq
    q_df <- common_model$modelfit$df - independent_model$modelfit$df
    q_p <- pchisq(q = q_chisq, df = q_df, lower.tail =FALSE)
    data.frame(q_chisq, q_df, q_p)
}

model_diffs <- bind_rows(
    lapply(traits, function(trait) {
        model_diff_chisq(trait_common.fit[[trait]], trait_independent.fit[[trait]])
        }
    ),
    .id = "Trait"
    ) %>% 
    mutate(FDR = p.adjust(q_p, method = "fdr"))

#merge with samplesize df first
model_diff_neff_df <- left_join(model_diffs, sample_sizes, by = c("Trait" = "Pheno"))
model_diff_neff_df$Significant <- ifelse(model_diff_neff_df$FDR < 0.05, "FDR < 0.05", "No")

require(ggrepel)
model_diff_neff_df

trait_common.fit[["Neu"]]$results

#create a dotplot with -log10(p) on the y axis and neff on the x axis

chisq_diff_plot <- ggplot(model_diff_neff_df, aes(x = Neff, y = -log10(FDR), label = Trait)) +
geom_point(aes(color = Significant), size = 2.5) +
scale_color_manual(values = c("red","grey")) +
theme_bw() +
geom_label_repel()

chisq_diff_plot

ggsave("../plots/jolien_paper/chisq_diff_bipolar.svg", chisq_diff_plot, width = 10, height = 10, create.dir = TRUE)

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

library(lavaanPlot)
library(lavaan)
library(semPlot)

for (trait in other_names){
    #define the two sem objects
    current_fit_common <- semPlotModel_GSEM(gsem.object = trait_common.fit[[trait]])
    current_fit_independent <- semPlotModel_GSEM(gsem.object = trait_independent.fit[[trait]])
    #extract the std estimate and the corresponding se to plot them both, for common and independent model
    se <- trait_common.fit[[trait]]$results$STD_Genotype_SE
    est <- trait_common.fit[[trait]]$results$STD_Genotype
    est.se <- paste0(round(as.numeric(est), 2),"\n (",round(as.numeric(se), 2),")")

    se_in <- trait_independent.fit[[trait]]$results$STD_Genotype_SE
    est_in <- trait_independent.fit[[trait]]$results$STD_Genotype
    est.se_in <- paste0(round(as.numeric(est_in), 2),"\n (",round(as.numeric(se_in), 2),")")
    #save the common model as an svg object, remove the last label as overlapping and known as factor artificially set to variance of 1
    svg(paste0("/local1/home/pazweifel/plots/jolien_paper/",trait,"_common_model.svg"), width = 11, height = 12)
    semPaths(current_fit_common, whatLabels = "std", layout = "tree", edge.color = "black", sizeMan = 12, sizeLat = 12, nCharNodes = 0, residuals = TRUE, style = "mx",
        edge.label.cex = 1,
        esize = 2,
        edgeLabels = est.se[1:12])
    dev.off()
    #save the indpendent mode as an svg object
    svg(paste0("/local1/home/pazweifel/plots/jolien_paper/",trait,"_independent_model.svg"), width = 12, height = 11)
    semPaths(current_fit_independent, whatLabels = "std", layout = "circle2", edge.color = "black", sizeMan = 8, sizeLat = 8, nCharNodes = 0, residuals = TRUE, style = "mx",
        edge.label.cex = 1,
        rotation = 1,
        esize = 2,
        edge.label.position = 0.45,
        edgeLabels = est.se_in)
    dev.off()
}

se <- base.fit$results$STD_Genotype_SE
est <- base.fit$results$STD_Genotype
est.se <- paste0(round(as.numeric(est), 2),"\n (",round(as.numeric(se), 2),")")

se_unstand <- base.fit$results$Unstand_SE
est_unstand <- base.fit$results$Unstand_Est
est.se_unstand <- paste0(round(as.numeric(est_unstand), 2),"\n (",round(as.numeric(se_unstand), 2),")")

se_in <- anchor.fit$results$STD_Genotype_SE
est_in <- anchor.fit$results$STD_Genotype
est.se_in <- paste0(round(as.numeric(est_in), 2),"\n (",round(as.numeric(se_in), 2),")")

fit_sem_base <- semPlotModel_GSEM(base.fit)
fit_sem_anchor <- semPlotModel_GSEM(anchor.fit)

sem_neu_independent <- semPlotModel_GSEM(trait_independent.fit[["T2D"]])
#svg(paste0("/local1/home/pazweifel/plots/jolien_paper/","Neu","_independent_model.svg"), width = 12, height = 11)
semPaths(sem_neu_independent,  whatLabels = "std", layout = "circle2", edge.color = "black", sizeMan = 8, sizeLat = 8, nCharNodes = 0, residuals = TRUE, style = "mx",
        edge.label.cex = 1,
        rotation = 1,
        edge.label.position = 0.45,
        esize = 2,
        edge.label.position = 0.5
        )
#dev.off()

positions <- sapply(c(1:length(est.se)), function(index) {
    if (index%%2==0){
        0.7
    }else {
        0.3}
}
)
semPaths(fit_sem_base, whatLabels = "std", layout = "tree", edge.color = "black", sizeMan = 12, sizeLat = 12, nCharNodes = 0, residuals = TRUE, style = "mx",
        edge.label.cex = 1.2,
        esize = 2,
        edgeLabels = est.se,
        edge.label.position = 0.5
        )

pdf("/local1/home/pazweifel/plots/jolien_paper/base_model_wse.pdf", width = 14, height = 14)
semPaths(fit_sem_base, whatLabels = "std", layout = "tree", edge.color = "black", sizeMan = 12, sizeLat = 12, nCharNodes = 0, residuals = TRUE, style = "mx",
        edge.label.cex = 1.2,
        esize = 2,
        edgeLabels = est.se,
        edge.label.position = 0.5
        )
dev.off()

semPaths(fit_sem_base, whatLabels = "std", layout = "tree", edge.color = "black", sizeMan = 12, sizeLat = 12, nCharNodes = 0, residuals = TRUE, style = "mx",
        edge.label.cex = 1.2,
        esize = 2,
        edgeLabels = est.se_unstand,
        edge.label.position = 0.5
        )

pdf("/local1/home/pazweifel/plots/jolien_paper/base_model_wse_unstand.pdf", width = 14, height = 14)
semPaths(fit_sem_base, whatLabels = "std", layout = "tree", edge.color = "black", sizeMan = 12, sizeLat = 12, nCharNodes = 0, residuals = TRUE, style = "mx",
        edge.label.cex = 1.2,
        esize = 2,
        edgeLabels = est.se_unstand,
        edge.label.position = 0.5
        )
dev.off()

base.fit$modelfit

anchor.fit$modelfit





semPaths(fit_sem_anchor, whatLabels = "std", layout = "tree", edge.color = "black", sizeMan = 12, sizeLat = 12, nCharNodes = 0, residuals = TRUE, style = "mx",
        edge.label.cex = 1.2,
        esize = 2,
        edgeLabels = est.se_in)

svg("/local1/home/pazweifel/plots/jolien_paper/anchor_model_wse.svg", width = 10, height = 12)
semPaths(fit_sem_anchor, whatLabels = "std", layout = "tree", edge.color = "black", sizeMan = 12, sizeLat = 12, nCharNodes = 0, residuals = TRUE, style = "mx",
        edge.label.cex = 1.2,
        esize = 2,
        edgeLabels = est.se_in)
dev.off()

trait.qtrait <- bind_rows(lapply(traits, function(trait) {
    QTrait(
    LDSCoutput = trait.covstructs[[trait]],
    indicators = mental_names,
    traits = trait)
    }),
    .id = "Trait")

trait.qtrait

help(QTrait)

cov2cor(trait_common.fit[["Neu"]]$resid_cov$'Model Implied Covariance Matrix')[1:5, 6]

#first extract the model implied covariances from the common model fit and convert to correlations
target_common_rg <- bind_rows(
    lapply(trait_common.fit, function(fit) {
        cov2cor(V = fit$resid_cov$'Model Implied Covariance Matrix')[1:5, 6]
        }
    ),
    .id = "trait")
    

#add the common pathway label and make the df longer so that you have one correlation per line to plot afterwards
target_common_rg <- target_common_rg %>% mutate(Pathway = "Common") %>% pivot_longer(ADHD:SCZ, names_to = "Mental", values_to = "rg")

target_common_rg
trait_independent.fit[[1]]$results

#now extract correlation of external phenotype with common factor form results table (F1~~E)
target_common_factor_rg <- bind_rows(
    lapply(trait_common.fit, function(fit){
        fit$results %>%
        filter(
            lhs == "F1",
            rhs != "F1",
            !rhs %in% mental_names
            ) %>%
        transmute(Mental = "F1", rg = STD_Genotype, se = as.numeric(STD_Genotype_SE)
                  )
        }
    ),
    .id = "trait"
    ) %>%
    mutate(Pathway = "Common")

#now the same for the independent model, extract from model-implied covariance matrix
target_independent_rg <- bind_rows(
    lapply(
        trait_independent.fit, function(fit){
            cov2cor(fit$resid_cov$'Model Implied Covariance Matrix')[1:5, 6]
            }),
    .id = "trait"
    ) %>% mutate(Pathway = "Independent") %>% pivot_longer(ADHD:SCZ, names_to = "Mental", values_to = "rg")

#for the independent model extract the correlation between factor and external trait form results
target_independent_rg_se <- bind_rows(
    lapply(trait_independent.fit, function(fit){
        fit$results %>% filter(
            lhs %in% mental_names,
            lhs != rhs
            ) %>% transmute(
            Mental = lhs, rg = STD_Genotype, se = as.numeric(STD_Genotype_SE))
        }),
        .id = "trait"
    ) %>% mutate(
    Pathway = "Independent")
        


target_sample_sizes <-
  bind_rows(
    lapply(trait.covstructs, function(covs) {
      # get sample sizes
      N <- covs$N
      # create a matrix the same size as the cov matrix
      Nmat <- covs$S
      # fill in lower triangle
      Nmat[TRUE] <- NA
      Nmat[lower.tri(Nmat, diag = TRUE)] <- N
      # extract sample size and names from matrix
      Ns <- Nmat[6, 1:5]
      N_names <- dimnames(Nmat)[[2]][1:5]
      # add total N for factor
      N_total <- sum(Ns)
      tibble(N = c(Ns, N_total), Mental = c(N_names, "F1"))
    }),
    .id = "trait"
  )

target_rg <- bind_rows(
    target_common_rg, 
    target_common_factor_rg, 
    target_independent_rg_se) %>% left_join(
    target_sample_sizes,
    by = c("trait", "Mental")) %>% left_join(
    model_diff_neff_df, by = c("trait" = "Trait"))

correlations_table <- target_rg |> arrange(desc(Pathway)) |> mutate(trait = reorder(trait, 1 - rg)) %>% select(-c(Nca, Nco, Neff, Significant))


correlation_plot <- ggplot(
  target_rg |> arrange(desc(Pathway)) |> mutate(trait = reorder(trait, 1 - rg)),
  aes(
    x = Mental, y = rg,
    ymin = rg + se * qnorm(0.025), ymax = rg + se * qnorm(0.975),
    shape = Pathway, colour = Pathway, size = N
  )
) +
  geom_linerange(aes(ymin = rg + se * qnorm(0.1), ymax = rg + se * qnorm(0.9)), linewidth = 0.75) +
  geom_pointrange() +
  facet_grid(rows = vars(trait)) +
  scale_x_discrete("Mental phenotype", limits = c("F1", rev(mental_names))) +
  scale_y_continuous(expression(r[g]), breaks = c(-0.25, 0, 0.25, 0.5, 0.75, 1)) +
  scale_size_area(breaks = c(2.5e4, 5e4, 1e5, 2.5e5, 5e5), max_size = 1) +
  scale_colour_discrete(
    type = c(
      microshades_cvd_palettes$micro_cvd_blue[4],
      microshades_cvd_palettes$micro_cvd_orange[2]
    )
  ) +
  coord_flip(ylim = c(-0.25, 1)) +
  theme_minimal() +
  theme(panel.grid.minor = element_blank(), panel.grid.major = element_line(linewidth = 0.1))

correlation_plot

ggsave("/local1/home/pazweifel/plots/jolien_paper/Correlation_comparison.svg", plot = correlation_plot, width = 10, height = 25, device = "svg") 

target_rg |>
    pivot_wider(names_from = c("Pathway"), values_from = c("rg", "se")) |>
    mutate(rgD = rg_Independent - rg_Common) |>
    select(-se_Common) |>
    arrange(desc(abs(rgD))) |>
    filter(abs(rgD) > se_Independent) |> select(-c(Nca, Nco, Neff, Significant))


rg_ordered_grouped_df <- target_rg %>% filter(Pathway == "Independent") %>% group_by(trait) %>% arrange(-rg, .by_group = TRUE) %>% subset(select = c(trait, Mental, rg))

rg_ordered_grouped_df

vec <- filter(rg_ordered_grouped_df, trait == "ADHD") %>% .$BD %>% .[1:(4-1)]
vec[5]

dynamic_model_1 <- "F1=~NA*Clinical+Community+BDI+BDII
{external_trait}~~{trait1}
{external_trait}~~{trait2}
{external_trait}~~{trait3}
"

dynamic_model_2 <- "F1~~1*F1
Clinical~~a*Clinical
a>0.0001
"

paste0(dynamic_model_1,dynamic_model_2)
dynamic_model_1
substr(dynamic_model_1, 1, nchar(dynamic_model_1))

sub <- nchar("{external_trait}~~{traitx}\n")
a <- substr(dynamic_model_1, 1, nchar(dynamic_model_1) - 2 * sub)

paste0(a,dynamic_model_2)
#trait 1 the one with the highest loading 

#lsit where i want to store the different follow up models with stepwise retracted traits
fitlist <- list(follow_up_1 = 0, follow_up_2 = 0, follow_up_3 = 0)

#list where i want to store the different follow up models with stepwise retracted traits
fitlist <- list(follow_up_1 = 0, follow_up_2 = 0, follow_up_3 = 0)

dynamic_model_1 <- "F1=~NA*Clinical+Community+BDI+BDII
F1~~1*F1
{trait1}~~{external_trait}
{trait2}~~{external_trait}
{trait3}~~{external_trait}
"

dynamic_model_2 <- "F1~b*{external_trait}
Clinical~~a*Clinical
a>0.0001
"

sub <- nchar("{external_trait}~~{traitx}\n")


for (index in 1:3){
    temp.fit <- lapply(traits, function(external_trait){
        #extract the internal traits ordered by rg; depending on index remove the one with the lowest correlation and store them in a vector
        trait_vector <- filter(rg_ordered_grouped_df, trait == external_trait) %>% .$BD %>% .[1:(4-index)]
        #assign values form vector to trait variables which match the placeholders in the dynamic model, if out of range jsut gets NA so no problem
        trait1 <- trait_vector[1]
        trait2 <- trait_vector[2]
        trait3 <- trait_vector[3]
        #adjust the model
        temp_model <- substr(dynamic_model_1, 1, nchar(dynamic_model_1) - (index -1) * sub)
        print(temp_model)
        #paste the traits into the model structure
        temp_model <- str_glue(temp_model)
        temp_model_2 <- str_glue(dynamic_model_2)
        print(temp_model)
        temp_model <- paste(temp_model, temp_model_2, sep = "\n")
        print(temp_model)
        usermodel(
                trait.covstructs[[external_trait]],
                estimation = "DWLS",
                model = temp_model,
                imp_cov = TRUE
                )
            }
        )
    fitlist[[index]] <- temp.fit
}


fitlist$follow_up_3[["ADHD"]]

for (i in fitlist){
    tempsem <- semPlotModel_GSEM(i[["Neu"]])
      semPaths(tempsem, whatLabels = "std", layout = "circle", edge.color = "black", sizeMan = 12, sizeLat = 12, nCharNodes = 0, residuals = TRUE, style = "mx",
        edge.label.cex = 1.2,
        rotation = 1,
        esize = 2,
        edge.label.position = 0.7)
    print(i[["Neu"]]$modelfit$chisq)
    }
print(trait_common.fit[["Neu"]]$modelfit$chisq)    
print(trait_independent.fit[["Neu"]]$modelfit$chisq)  


