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
Clin,27196,43792,67108.0876767904
Comm,32091,737230,123009.500221624
BDI,25060,449978,94951.971673845
BDII,6781,364075,26628.0451172423
ADHD,19099,34194,45684.86
AN,16992,55525,46321.9
ALC,NA,NA,26853.43
ASD,18381,27969,44366.62
MDD,357636,1281936,1118502.79
BMI,NA,NA,681275
CAD,95830,385536,285937.77
Chrono,NA,NA,449734
CRP,NA,NA,204402
Height,NA,NA,693529
Neu,NA,NA,390278
PD,2147,7760,6665.06
PTSD,23212,151447,140475
SCZ,53386,77258,117498.3
Smoking,NA,NA,263954
T2D,80154,853816,251739.50
PartEM,NA,NA,356560.32
PartMH,NA,NA,311999.52
PartSF,NA,NA,449041.88
"
)

#BIP,59287,781022,163367

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

#create a matrices directory
#dir.create("matrices", showWarnings = FALSE)
#create a path for a covstruct object where you store the covariance structure of the bipolar disorders and store it in the matrices directory
covstruct_bip <- here::here("matrices",paste("BD", "R", sep = "."))
#extract population prevalence 
pop_prev <- as.vector(all_prevalences[bip_names])
#fill the samp_prev vector depending on the pop_prev vector, if NA use NA and if not use 0.5
samp_prev <- as.vector(ifelse(is.na(pop_prev), NA, 0.5))
#check if the covstruct file already exists to not do the analysis unnecessarily, first time no file there
if (!file.exists(covstruct_bip)){
    covstruct_bip_ldsc <- ldsc(
        traits = samples,
        trait.names = sample_names,
        sample.prev = samp_prev,
        population.prev = pop_prev,
        ld = ld,
        wld = weights,
        stand = TRUE
    )
    #store the file in the previously defined path, now the path is, contorl tell how to deparse
    dput(covstruct_bip_ldsc, covstruct_bip, control = c("all", "digits17"))
}

# do the same thing but include the other phenotypes one by one and name them after the current phenotype used
for (phenotype in other_names){
    covstruct_pheno <- here::here("matrices",paste(phenotype,"R", sep = "."))
    #extract population prevalences
    pop_prev <- c(all_prevalences[bip_names], all_prevalences[phenotype])
    #fill the sample prevalences depending on pop vector
    samp_prev <- as.vector(ifelse(is.na(pop_prev), NA, 0.5))
    #concatenate the sample vector with the bip filepaths and the specific single phenotype filepath
    if (phenotype == "MDD"){
        samples_plus1 <- c(samples, "/local1/scratch/pazweifel/jolien_paper_sumstats/munged_sumstats/mdd_munged.sumstats.gz")
    } else {
        samples_plus1 <- c(samples, paste0(other_basepath,phenotype,".sumstats.gz"))
    }
    #check if file already exists
    if (!file.exists(covstruct_pheno)){
        covstruct_ldsc <- ldsc(
            traits = samples_plus1,
            trait.names = c(sample_names, phenotype),
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
covstruct_base <- dget(here::here("matrices/BD.R"))
#define the model using NA*Clin to say that this parameter loading has to be freely estimated, we fix the factor variance to 1 and ensure residual genetic variance of phenotypes isnt negative
base_model <- "BD=~NA*Clin+Comm+BDI+BDII
BD~~1*BD
Clin~~a*Clin
a > .0001"

base.fit <- usermodel(covstruct_base,
                      estimation = "DWLS",
                      model = base_model,
                      imp_cov = TRUE
)

base.fit

base.fit$results

sstand <- covstruct_base$S_Stand
vstand <- covstruct_base$V_Stand

#extract the standard error
vstand_diag <- sapply(c(1:10), function(index) vstand[index, index])
                      
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
#create the ggplot
correlation_plot <- ggplot(expanded_df, aes(x = x, y = y, fill = values)) +
scale_fill_gradient(low = "darkorchid", high = "orange", limits = c(0, max(expanded_df$values))) +
theme_minimal(base_size = 20) +
geom_tile() +
theme(axis.title = element_blank()) +
labs(fill = expression(r[g])) +
geom_text(aes(label = values)) +
geom_text(aes(label = paste0("(",se,")"), vjust = 2))

correlation_plot

ggsave("/local1/home/pazweifel/plots/bip_paper/correlation_matrix.pdf", correlation_plot, width = 7, height = 5, device = "pdf")

write_csv(expanded_df, "/local1/home/pazweifel/plots/bip_paper/correlations_base.csv")

#same thing but with fixed loading of 1 on the clinical, factor variance freely estimated
anchor_model <- "BD=~1*Clin+Comm+BDI+BDII
BD~~BD
Clin~~a*Clin
a > 0.0001
"
anchor.fit <- usermodel(covstruct_base,
                        estimation = "DWLS",
                        model = anchor_model,
                        imp_cov = TRUE
)

anchor.fit$modelfit

#define the common and the independent pathway models and introduce a placeholder which can filled by the repsective phenotype corresponding to the used covariance matrix
#set the factor variance to 1
common.model <- "BD=~NA*Clin+Comm+BDI+BDII
BD~~{trait}
BD~~1*BD
Clin~~a*Clin
a > 0.0001
"

independent.model <- "BD=~NA*Clin+Comm+BDI+BDII
Clin~~{trait}
Comm~~{trait}
BDI~~{trait}
BDII~~{trait}
BD~~1*BD
Clin~~a*Clin
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
                           function(trait) dget(here::here("matrices",paste(trait,"R", sep = ".")))
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
model_diff_neff_df %>% subset(select = c(Trait, q_chisq, q_df, q_p, FDR)) %>% mutate(FDR_stars = case_when(
    FDR < 0.001 ~ "***",
    FDR < 0.01 ~ "**",
    FDR < 0.05 ~ "*",
    TRUE ~ "")) 

#create a dotplot with -log10(p) on the y axis and neff on the x axis
library(scales)

chisq_diff_plot <- ggplot(model_diff_neff_df, aes(x = Neff, y = -log10(FDR), label = Trait)) +
geom_point(aes(color = Significant), size = 2.5) +
scale_color_manual(values = c("red","grey"), labels = c("FDR < 0.05" = expression(P[FDR]*"< 0.05"))) +
theme_bw(base_size = 20) +
geom_label_repel() +
scale_x_continuous(labels = label_comma()) +
labs(x = expression(N[eff]), y = expression(-log[10](P[FDR]))) + 
theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
panel.background = element_blank(), axis.line = element_line(colour = "black"), legend.position = "bottom") 



chisq_diff_plot

ggsave("../plots/bip_paper/chisq_diff_bipolar.pdf", chisq_diff_plot, width = 10, height = 10)

write_csv(model_diff_neff_df, "/local1/home/pazweifel/plots/bip_paper/chisquared_diffs.csv")

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

current_fit_common <- semPlotModel_GSEM(gsem.object = trait_common.fit[["ADHD"]])

se <- trait_common.fit[["ADHD"]]$results$STD_Genotype_SE
est <- trait_common.fit[["ADHD"]]$results$STD_Genotype
est.se <- paste0(round(as.numeric(est), 2),"\n (",round(as.numeric(se), 2),")")

semPaths(current_fit_common, layout = "tree", edge.color = "black", sizeMan = 12, sizeLat = 12, nCharNodes = 0, residuals = TRUE, style = "mx",
        edge.label.cex = 1.2,
        esize = 2,
        edgeLabels = est.se[1:(length(est.se)-1)])

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
    svg(paste0("/local1/home/pazweifel/plots/bip_paper/",trait,"_common_model.svg"), width = 11, height = 12)
    semPaths(current_fit_common, whatLabels = "std", layout = "tree", edge.color = "black", sizeMan = 12, sizeLat = 12, nCharNodes = 0, residuals = TRUE, style = "mx",
        edge.label.cex = 1.2,
        esize = 2,
        edgeLabels = est.se[1:10])
    dev.off()
    #save the indpendent mode as an svg object
    svg(paste0("/local1/home/pazweifel/plots/bip_paper/",trait,"_independent_model.svg"), width = 12, height = 11)
    semPaths(current_fit_independent, whatLabels = "std", layout = "circle", edge.color = "black", sizeMan = 12, sizeLat = 12, nCharNodes = 0, residuals = TRUE, style = "mx",
        edge.label.cex = 1.2,
        rotation = 1,
        esize = 2,
        edgeLabels = est.se_in)
    dev.off()
}

se <- base.fit$results$STD_Genotype_SE
est <- base.fit$results$STD_Genotype
est.se <- paste0(round(as.numeric(est), 2),"\n (",round(as.numeric(se), 2),")")

se_in <- anchor.fit$results$STD_Genotype_SE
est_in <- anchor.fit$results$STD_Genotype
est.se_in <- paste0(round(as.numeric(est_in), 2),"\n (",round(as.numeric(se_in), 2),")")

fit_sem_base <- semPlotModel_GSEM(base.fit)
fit_sem_anchor <- semPlotModel_GSEM(anchor.fit)

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

svg("/local1/home/pazweifel/plots/bip_paper/base_model_wse.svg", width = 10, height = 12)
semPaths(fit_sem_base, whatLabels = "std", layout = "tree", edge.color = "black", sizeMan = 12, sizeLat = 12, nCharNodes = 0, residuals = TRUE, style = "mx",
        edge.label.cex = 1.2,
        esize = 2,
        edgeLabels = est.se,
        edge.label.position = 0.5
        )
dev.off()

semPaths(fit_sem_anchor, whatLabels = "std", layout = "tree", edge.color = "black", sizeMan = 12, sizeLat = 12, nCharNodes = 0, residuals = TRUE, style = "mx",
        edge.label.cex = 1.2,
        esize = 2,
        edgeLabels = est.se_in)

svg("/local1/home/pazweifel/plots/bip_paper/anchor_model_wse.svg", width = 10, height = 12)
semPaths(fit_sem_anchor, whatLabels = "std", layout = "tree", edge.color = "black", sizeMan = 12, sizeLat = 12, nCharNodes = 0, residuals = TRUE, style = "mx",
        edge.label.cex = 1.2,
        esize = 2,
        edgeLabels = est.se_in)
dev.off()

trait.qtrait <- bind_rows(lapply(traits, function(trait) {
    QTrait(
    LDSCoutput = trait.covstructs[[trait]],
    indicators = sample_names,
    traits = trait)
    }),
    .id = "Trait")

trait.qtrait

cov2cor(trait_common.fit[["ADHD"]]$resid_cov$'Model Implied Covariance Matrix')[1:4, 5]

#first extract the model implied covariances from the common model fit and convert to correlations
target_common_rg <- bind_rows(
    lapply(trait_common.fit, function(fit) {
        cov2cor(V = fit$resid_cov$'Model Implied Covariance Matrix')[1:4, 5]
        }
    ),
    .id = "trait")
    

#add the common pathway label and make the df longer so that you have one correlation per line to plot afterwards
target_common_rg <- target_common_rg %>% mutate(Pathway = "Common") %>% pivot_longer(Clin:BDII, names_to = "BD", values_to = "rg")

trait_independent.fit[[1]]$results

#now extract correlation of external phenotype with common factor form results table (BD~~E)
target_common_factor_rg <- bind_rows(
    lapply(trait_common.fit, function(fit){
        fit$results %>%
        filter(
            lhs == "BD",
            rhs != "BD",
            !rhs %in% bip_names
            ) %>%
        transmute(BD = "BD", rg = STD_Genotype, se = as.numeric(STD_Genotype_SE)
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
            cov2cor(fit$resid_cov$'Model Implied Covariance Matrix')[1:4, 5]
            }),
    .id = "trait"
    ) %>% mutate(Pathway = "Independent") %>% pivot_longer(Clin:BDII, names_to = "BD", values_to = "rg")

#for the independent model extract the correlation between factor and external trait form results
target_independent_rg_se <- bind_rows(
    lapply(trait_independent.fit, function(fit){
        fit$results %>% filter(
            lhs %in% bip_names,
            lhs != rhs
            ) %>% transmute(
            BD = lhs, rg = STD_Genotype, se = as.numeric(STD_Genotype_SE))
        }),
        .id = "trait"
    ) %>% mutate(
    Pathway = "Independent")
        

#extract the observed correlation and se
test <- trait.covstructs[["MDD"]]
S <- test$S_Stand
V <- test$V_Stand
rgs <- S[5, 1:4]
print(rgs)
names <- dimnames(S)[[2]][1:4]
print(names)
se_index <- c(5,9,12,14)
ses <- sapply(se_index, function(index) V[index, index])
print(ses)
tibble(BD = names, rg = rgs, se = ses)
#basically the same as the ones form the independent model, so not necessary

#extract the heritabilities and their se
#for the base phenotypes
temp <- covstruct_base
S <- temp$S
V <- temp$V
herit_base <- sapply(c(1:4), function(index) S[index,index])
                  
index_se <- c(1, 5, 8, 10)
se_base <- sapply(index_se, function(index) sqrt(V[index, index]))
           
names <- dimnames(S)[[2]][1:4]
                
base_heritabilities <- tibble(trait = names, heritability = herit_base, se = se_base, group = "internal")

ext_heritabilities <- bind_rows(
                            lapply(traits, function(trait) {
                                temp <- trait.covstructs[[trait]]
                                heritability <- temp$S[5,5]
                                se <- sqrt(temp$V[15,15])
                                tibble(heritability = heritability, se = se, group = "external")
                                }),
                                .id = "trait")
                                

all_heritabilities <- bind_rows(base_heritabilities, ext_heritabilities)
all_heritabilities$trait <- factor(all_heritabilities$trait, levels = unique(all_heritabilities$trait))
all_heritabilities$group <- factor(all_heritabilities$group, levels = c("internal", "external"))

heritabilities_bip_paper <- ggplot(all_heritabilities, aes(x = trait, y = heritability, fill = group)) + 
geom_col(color = "black") +
theme_bw(base_size = 20) +
geom_errorbar(aes(ymin = heritability - se, ymax = heritability + se), width = 0.6) +
theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.title = element_blank()) +
labs(x = "Trait", y = expression(italic(h^2))) +
scale_y_continuous(limits = c(0, 0.6), expand = c(0,0)) +
scale_fill_manual(values = c("darkorchid", "orange"))

heritabilities_bip_paper

ggsave("/local1/home/pazweifel/plots/bip_paper/heritabilities_all.pdf", heritabilities_bip_paper, device = pdf, width = 22, height = 10)

write_csv(all_heritabilities, "/local1/home/pazweifel/plots/bip_paper/heritabilities.csv")


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
      Ns <- Nmat[5, 1:4]
      N_names <- dimnames(Nmat)[[2]][1:4]
      # add total N for factor
      N_total <- sum(Ns)
      tibble(N = c(Ns, N_total), BD = c(N_names, "BD"))
    }),
    .id = "trait"
  )

target_sample_sizes

target_rg <- bind_rows(
    target_common_rg, 
    target_common_factor_rg, 
    target_independent_rg_se) %>% left_join(
    target_sample_sizes,
    by = c("trait", "BD")) %>% left_join(
    model_diff_neff_df, by = c("trait" = "Trait"))

correlation_table <- target_rg |> arrange(desc(Pathway)) |> mutate(trait = reorder(trait, 1 - rg)) %>% select(-c(Nca, Nco, Neff, Significant))
write_csv(x = correlation_table, file = "/local1/home/pazweifel/plots/bip_paper/correlations_table.csv")

target_rg_noheight <- filter(target_rg, trait!= "Height")

mdd_correlations <- read_delim("/local1/scratch/pazweifel/sumstats_ambits/correlations_md.csv", delim = ";")

head(mdd_correlations)
head(target_rg)

mdd_correlations_adjusted <- mdd_correlations %>% mutate(Variable = MDD, Analysis = "MDD") %>% select(-MDD)
bd_correlations_adjusted <- target_rg %>% mutate(Variable = BD, Analysis = "BD") %>% select(-c(BD, Nca, Nco, Neff, Significant))

combined_rgs <- bind_rows(mdd_correlations_adjusted, bd_correlations_adjusted)
head(combined_rgs)

correlation_plot_thes <- ggplot(
  target_rg |> arrange(desc(Pathway)) |> mutate(trait = reorder(trait, 1 - rg)),
  aes(
    x = BD, y = rg,
    ymin = rg + se * qnorm(0.025), ymax = rg + se * qnorm(0.975),
    shape = Pathway, colour = Pathway, size = N
  )
) +
  geom_linerange(aes(ymin = rg + se * qnorm(0.1), ymax = rg + se * qnorm(0.9)), linewidth = 0.75) +
  geom_pointrange() +
  facet_grid(rows = vars(trait)) +
  scale_x_discrete("BD phenotype", limits = c("BD", rev(bip_names))) +
  scale_y_continuous(expression(r[g]), breaks = c(-0.25, 0, 0.25, 0.5, 0.75, 1)) +
  scale_size_area(breaks = c(2.5e4, 5e4, 1e5, 2.5e5, 5e5), max_size = 1) +
  scale_colour_discrete(
    type = c(
      microshades_cvd_palettes$micro_cvd_blue[4],
      microshades_cvd_palettes$micro_cvd_orange[2]
    )
  ) +
  coord_flip(ylim = c(-0.25, 1)) +
  theme_minimal(base_size = 30) +
  theme(panel.grid.minor = element_blank(), panel.grid.major = element_line(linewidth = 0.1), legend.position = "bottom")
correlation_plot_thes

ggsave("/local1/home/pazweifel/plots/bip_paper/correlation_plot_use.pdf", plot = correlation_plot_thes, device = "pdf", height = 30, width = 25)

correlation_plot_bd <- ggplot(
  target_rg_noheight |> arrange(desc(Pathway)) |> mutate(trait = reorder(trait, 1 - rg)),
  aes(
    x = BD, y = rg,
    ymin = rg + se * qnorm(0.025), ymax = rg + se * qnorm(0.975),
    shape = Pathway, colour = Pathway, size = N
  )
) +
  geom_linerange(aes(ymin = rg + se * qnorm(0.1), ymax = rg + se * qnorm(0.9)), linewidth = 0.75) +
  geom_pointrange() +
  facet_grid(rows = vars(trait)) +
  scale_x_discrete("BD phenotype", limits = c("BD", rev(bip_names))) +
  scale_y_continuous(expression(r[g]), breaks = c(-0.25, 0, 0.25, 0.5, 0.75, 1)) +
  scale_size_area(breaks = c(2.5e4, 5e4, 1e5, 2.5e5, 5e5), max_size = 1) +
  scale_colour_discrete(
    type = c(
      microshades_cvd_palettes$micro_cvd_blue[4],
      microshades_cvd_palettes$micro_cvd_orange[2]
    )
  ) +
  coord_flip(ylim = c(-0.25, 1)) +
  theme_minimal(base_size = 30) +
  theme(panel.grid.minor = element_blank(), panel.grid.major = element_line(linewidth = 0.1))



correlation_plot_md <- ggplot(
  mdd_correlations |> arrange(desc(Pathway)) |> mutate(trait = reorder(trait, 1 - rg)),
  aes(
    x = MDD, y = rg,
    ymin = rg + se * qnorm(0.025), ymax = rg + se * qnorm(0.975),
    shape = Pathway, colour = Pathway, size = N
  )
) +
  geom_linerange(aes(ymin = rg + se * qnorm(0.1), ymax = rg + se * qnorm(0.9)), linewidth = 0.75) +
  geom_pointrange() +
  facet_grid(rows = vars(trait)) +
  scale_x_discrete("MD phenotype", limits = c("MD", "Help", "Quest", "EHR", "Clin")) +
  scale_y_continuous(expression(r[g]), breaks = c(-0.25, 0, 0.25, 0.5, 0.75, 1)) +
  scale_size_area(breaks = c(2.5e4, 5e4, 1e5, 2.5e5, 5e5), max_size = 1) +
  scale_colour_discrete(
    type = c(
      microshades_cvd_palettes$micro_cvd_blue[4],
      microshades_cvd_palettes$micro_cvd_orange[2]
    )
  ) +
  coord_flip(ylim = c(-0.25, 1)) +
  theme_minimal(base_size = 30) +
  theme(panel.grid.minor = element_blank(), panel.grid.major = element_line(linewidth = 0.1))

correlation_plot_md_legendbottom <- correlation_plot_md + theme(legend.position = "bottom")

unique(mdd_correlations$MDD)

library(cowplot)

shared_legend <- get_legend(correlation_plot_md)
shared_legend_bottom <- get_legend(correlation_plot_md_legendbottom)

md_wo_legend <- correlation_plot_md + theme(legend.position = "None")
bd_wo_legend <- correlation_plot_bd + theme(legend.position = "None")

correlation_plots <- plot_grid(plotlist = c(md_wo_legend, bd_wo_legend), ncol = 2)

correlation_plots_wlegend <- plot_grid(correlation_plots, shared_legend, ncol = 2, rel_widths = c(1, 0.25))
correlation_plots_wlegendbottom <- plot_grid(correlation_plots, shared_legend_bottom, nrow = 2, rel_heights = c(1, 0.05))

ggsave("/local1/home/pazweifel/plots/bip_paper/combined_correlation_plot.pdf", plot = correlation_plots_wlegend, device = "pdf", height = 30, width = 25)

ggsave("/local1/home/pazweifel/plots/bip_paper/combined_correlation_plot_legendbottom.pdf", plot = correlation_plots_wlegendbottom, device = "pdf", height = 35, width = 25)

correlation_plot_combined <- ggplot(
  combined_rgs |> arrange(desc(Pathway)) |> mutate(trait = reorder(trait, 1 - rg)),
  aes(
    x = Variable, y = rg,
    ymin = rg + se * qnorm(0.025), ymax = rg + se * qnorm(0.975),
    shape = Pathway, colour = Pathway, size = N
  )
) +
  geom_linerange(aes(ymin = rg + se * qnorm(0.1), ymax = rg + se * qnorm(0.9)), linewidth = 0.75) +
  geom_pointrange() +
  facet_grid(rows = vars(trait), cols = vars(Analysis)) +
  #scale_x_discrete("Phenotype", limits = c("BD", rev(bip_names))) +
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

ggsave("/local1/home/pazweifel/plots/bip_paper/Correlation_comparison.svg", plot = correlation_plot, width = 10, height = 20, device = "svg") 

target_rg |>
    pivot_wider(names_from = c("Pathway"), values_from = c("rg", "se")) |>
    mutate(rgD = rg_Independent - rg_Common) |>
    select(-se_Common) |>
    arrange(desc(abs(rgD))) |>
    filter(abs(rgD) > se_Independent) |> select(-c(Nca, Nco, Neff, Significant))


rg_ordered_grouped_df <- target_rg %>% filter(Pathway == "Independent") %>% group_by(trait) %>% arrange(-rg, .by_group = TRUE) %>% subset(select = c(trait, BD, rg))

rg_ordered_grouped_df

vec <- filter(rg_ordered_grouped_df, trait == "ADHD") %>% .$BD %>% .[1:(4-1)]

dynamic_model_1 <- "BD=~NA*Clin+Comm+BDI+BDII
{external_trait}~~{trait1}
{external_trait}~~{trait2}
{external_trait}~~{trait3}
"

dynamic_model_2 <- "BD~~1*BD
Clin~~a*Clin
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

dynamic_model_1 <- "BD=~NA*Clin+Comm+BDI+BDII
BD~~1*BD
{trait1}~~{external_trait}
{trait2}~~{external_trait}
{trait3}~~{external_trait}
"

dynamic_model_2 <- "Clin~~a*Clin
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
    #tempsem <- semPlotModel_GSEM(i[["Neu"]])
      #semPaths(tempsem, whatLabels = "std", layout = "circle", edge.color = "black", sizeMan = 12, sizeLat = 12, nCharNodes = 0, residuals = TRUE, style = "mx",
        #edge.label.cex = 1.2,
        #rotation = 1,
        #esize = 2,
        #edge.label.position = 0.7)
    print(i[["Neu"]]$modelfit)
    
    }
print(trait_common.fit[["Neu"]]$modelfit)    
print(trait_independent.fit[["Neu"]]$modelfit)  

a <- fitlist[[1]][["Neu"]]$modelfit
a$version <- 1
a

independent_modelfit <- bind_rows(lapply(traits, function(trait) {
    trait_independent.fit[[trait]]$modelfit %>% mutate(version = "independent")
    }),
    .id = "trait")
            
common_modelfit <- bind_rows(lapply(traits, function(trait) {
    trait_common.fit[[trait]]$modelfit %>% mutate(version = "common")
    }),
    .id = "trait")


modelfit_list <- list(0,0,0)
for (i in 1:3){
    temp.modelfit <- bind_rows(
        lapply(traits, function(trait) {
            temp <- fitlist[[i]][[trait]]$modelfit
            #print(temp)
            temp$version <- as.character(i)
            #print(temp)
            temp}),
        .id = "trait")
    modelfit_list[[i]] <- temp.modelfit
    }
all_modelfits <- bind_rows(modelfit_list, independent_modelfit, common_modelfit)



all_modelfits$version <- factor(all_modelfits$version, levels = c("common", "independent", "1", "2", "3"))

ggplot(all_modelfits, aes(x = version, y = AIC, fill = trait, group = trait, color = trait)) +
geom_point() + theme_bw() + scale_y_log10() + geom_line()

ggplot(all_modelfits, aes(x = version, y = chisq, fill = trait, group = trait, color = trait)) +
geom_point() + theme_bw() + scale_y_log10() + geom_line()

ggplot(all_modelfits, aes(x = version, y = CFI, fill = trait, group = trait, color = trait)) +
geom_point() + theme_bw() + geom_line()

ggplot(all_modelfits, aes(x = version, y = SRMR, fill = trait, group = trait, color = trait)) +
geom_point() + theme_bw() + geom_line()

colnames(trait.qtrait)

their_follow_up_model <- "BD=~NA*Clin+Comm+BDI+BDII
BD~~1*BD
Clin~~a*Clin
a > 0.0001
BD~~{trait}
BDII ~ b*{trait}
"

#extracting the relevant traits
follow_up_names <- trait.qtrait[trait.qtrait$Unconstrained_paths != "None",]$Trait
traits[follow_up_names]

their_follow_up_models.fit <- lapply(traits[follow_up_names], function(trait) {
                                usermodel(covstruc = trait.covstructs[[trait]],
                                          estimation = "DWLS",
                                          model = str_glue(their_follow_up_model),
                                          imp_cov = TRUE
    )     
})



their_follow_up_models.fit[["MDD"]]$modelfit
trait_common.fit[["MDD"]]

for (trait in follow_up_names){
    #currentmodel
    temp <- their_follow_up_models.fit[[trait]]
    #create the semmodel used in the semPaths function
    tempsem <- semPlotModel_GSEM(temp)
    #extract the standardized estimates and the errors from results
    est_temp <- temp$results$STD_Genotype
    se_temp <- temp$results$STD_Genotype_SE

    labels_temp <- paste(round(as.numeric(est_temp), 2),"\n (", round(as.numeric(se_temp), 2),")")
    #create the object in which the plot can be deposited
    svg(filename = paste0("/local1/home/pazweifel/plots/bip_paper/",trait,"_follow_up_model.svg"), width = 11, height = 10)
    semPaths(tempsem, whatLabels = "std", layout = "circle", edge.color = "black", sizeMan = 12, sizeLat = 12, nCharNodes = 0, residuals = TRUE, style = "mx",
        edge.label.cex = 1,
        rotation = 1,
        esize = 2,
        edge.label.position = 0.5,
        edgeLabels = labels_temp[1:11])
    dev.off()
}

#make a named list out of follow_up_names to use for lapply
names(follow_up_names) <- follow_up_names

fit_table_follow <- bind_rows(
                lapply(follow_up_names, function(trait) {
                    temp_follow <- their_follow_up_models.fit[[trait]]$modelfit
                    temp_follow$Model <- "Follow-up"
                    temp_follow
                    }),
                    .id = "Trait")

fit_table_common <- bind_rows(
                        lapply(follow_up_names, function(trait) {
                            temp_follow <- trait_common.fit[[trait]]$modelfit
                            temp_follow$Model <- "Common"
                            temp_follow
                            }),
                    .id = "Trait")

fit_table_independent <- bind_rows(
                            lapply(follow_up_names, function(trait) {
                                temp_follow <- trait_independent.fit[[trait]]$modelfit
                                temp_follow$Model <- "Independent"
                                temp_follow
                                }),
                    .id = "Trait")
                            

follow_up_complete <- bind_rows(fit_table_follow, fit_table_common, fit_table_independent) %>% subset(select = -c(chisq, df, p_chisq)) %>% pivot_longer(cols = -c(Trait, Model), names_to = "Metric", values_to = "Value")

follow_up_complete

modelfit_follow_up_plot <- ggplot(follow_up_complete, aes(x = Trait, y = Value, fill = Model)) +
geom_col(position = position_dodge(), color = "black") + 
theme_bw() + 
facet_wrap(facets = vars(Metric), nrow = 1, ncol = 3, scales = "free") +
scale_fill_manual(values = c("slategray4","orange", "darkorchid")) +
scale_y_continuous(expand = expansion(mult = c(0,0.05))) +
theme(axis.title.y = element_blank())

modelfit_follow_up_plot

ggsave("/local1/home/pazweifel/plots/bip_paper/modelfit_follow_up.svg", modelfit_follow_up_plot, width = 8, height = 6)

model_fit_follow_up_wide <- pivot_wider(follow_up_complete, values_from = "Value", names_from = "Model")

model_fit_follow_up_wide %>% select(Trait, Metric, Common, everything()) %>% mutate(Common = round(Common, 3), `Follow-up` = round(`Follow-up`, 3), Independent = round(Independent, 3))

##### qtrait with adjusted thresholds
trait.qtrait_adjusted <- bind_rows(lapply(traits, function(trait) {
    QTrait(
    LDSCoutput = trait.covstructs[[trait]],
    indicators = sample_names,
    traits = trait,
    mresid = 0.1,
    mresidthreshold = 0.05,
    lsrmr = 0.1,
    lsrmrthreshold = 0.05
    )
    }),
    .id = "Trait")

follow_up_model_core <- "F1=~NA*Clin+Comm+BDI+BDII
F1~~1*F1
Clin~~a*Clin
a > 0.0001
F1~~{trait}
"
internal_1_model <- "\n{internal1} ~ b*{trait}"
internal_2_model <- "\n{internal2} ~ c*{trait}"
internal_3_model <- "\n{internal3} ~ d*{trait}"
internal_4_model <- "\n{internal4} ~ e*{trait}"

internal_models <- c(internal_1_model, internal_2_model, internal_3_model, internal_4_model)

#extracting the relevant traits
follow_up_names_adjusted <- trait.qtrait_adjusted[trait.qtrait_adjusted$Unconstrained_paths != "None",]$Trait
traits[follow_up_names_adjusted]
follow_up_table_adjusted <- trait.qtrait_adjusted[trait.qtrait_adjusted$Unconstrained_paths != "None",]

their_follow_up_models.fit_adjusted <- lapply(traits[follow_up_names_adjusted], function(trait) {
                                outliers <- str_split_1(follow_up_table_adjusted[follow_up_table_adjusted$Trait == trait,][1,]$Unconstrained_paths, pattern = ",")
                                #assign outliers to the variables
                                internal1 <- outliers[1]
                                internal2 <- outliers[2]
                                internal3 <- outliers[3]
                                internal4 <- outliers[4]

                                combined_model <- follow_up_model_core
                                for (i in 1:length(outliers)){
                                    combined_model <- paste0(combined_model, internal_models[i])
                                    
                                }
                                #fit the follow up model with the pasted construct
                                usermodel(covstruc = trait.covstructs[[trait]],
                                          estimation = "DWLS",
                                          model = str_glue(combined_model),
                                          imp_cov = TRUE
    )     
})



for (trait in follow_up_names){
    #currentmodel
    temp <- their_follow_up_models.fit_adjusted[[trait]]
    #create the semmodel used in the semPaths function
    tempsem <- semPlotModel_GSEM(temp)
    #extract the standardized estimates and the errors from results
    est_temp <- temp$results$STD_Genotype
    se_temp <- temp$results$STD_Genotype_SE
   
    labels_temp <- paste(round(as.numeric(est_temp), 2),"\n (", round(as.numeric(se_temp), 2),")")
    #create the object in which the plot can be deposited
    svg(filename = paste0("/local1/home/pazweifel/plots/bip_paper/",trait,"adj_threshold_follow_up_model.svg"), width = 11, height = 10)
    semPaths(tempsem, whatLabels = "std", layout = "circle", edge.color = "black", sizeMan = 12, sizeLat = 12, nCharNodes = 0, residuals = TRUE, style = "mx",
        edge.label.cex = 1,
        rotation = 1,
        esize = 2,
        edge.label.position = 0.5,
        edgeLabels = labels_temp[1:(length(est_temp)-1)])
    dev.off()
}

bias_model_adhd <- "
Liability=~NA*Clin+Comm+ADHD
Bias=~NA*ADHD+Comm
Liability ~~ 0*Bias
Clin ~~ 0*Clin
Comm ~~ 0*Comm
"

#gives a negative cfi, huge standard error for adhd estimators

bias_model_yazheng_adhd <- "
Liability =~ NA*Clin+ADHD
Bias =~ NA*ADHD
Liability ~~ 0*Bias
Clin ~~ 0*Clin
ADHD ~~ 0*ADHD
"

usermodel(trait.covstructs$ADHD, estimation = "DWLS", model = bias_model_adhd, std.lv = TRUE, imp_cov = TRUE, CFIcalc = TRUE)

usermodel(trait.covstructs$ADHD, estimation = "DWLS", model = bias_model_yazheng_adhd, std.lv = TRUE, imp_cov = TRUE, CFIcalc = TRUE)

#calcualte one big covariance matrix involving the clinical bipolar trait and all external ones
covstruct_bias_all <- here::here("matrices/Clinical_all_external.R")

pop_prev <- c(bip_prevalences["Clin"], other_prevalences)

samp_prev <- as.vector(ifelse(is.na(pop_prev), NA, 0.5))

paths <- c(samples[1])

for (name in names(other_prevalences)){
    if (name == "MDD"){
        paths <- c(paths, "/local1/scratch/pazweifel/jolien_paper_sumstats/munged_sumstats/mdd_munged.sumstats.gz")
        }
    else{
        paths <- c(paths, paste0(other_basepath, name, ".sumstats.gz"))
        }
}

    
if (!file.exists(covstruct_bias_all)){
    covstruct_bias_all_ldsc <- ldsc(
        traits = paths,
        trait.names = names(pop_prev),
        sample.prev = samp_prev,
        population.prev = pop_prev,
        ld = ld, 
        wld = weights,
        stand = TRUE
        )
    dput(covstruct_bias_all_ldsc, covstruct_bias_all, control = c("all", "digits17"))
}

bias_all.covstruct <- dget("matrices/Clinical_all_external.R")

#similar to what we had above, one common bias 
model <- "
Liability =~ NA*Clin + ADHD + ALC + AN + ASD + BMI + CAD + Chrono + CRP + Height + Neu + PD + PTSD + SCZ + Smoking + T2D + PartEM + PartMH + PartSF
Bias =~ NA*ADHD + ALC + AN + ASD + BMI + CAD + Chrono + CRP + Height + Neu + PD + PTSD + SCZ + Smoking + T2D + PartEM + PartMH + PartSF

Liability ~~ 0*Bias
Bias ~~ 1*Bias
Liability ~~ 1*Liability

Clin ~~ 0*Clin
"
#idea: liability factor captures pleiotropic factors between phenotypes and clinical, "true" genetic correlation, bias factor captures shared bias variants between external phenotypes, and freely estimated residual variances captures phenotype specific variances

all_bias.fit <- usermodel(covstruc = bias_all.covstruct, estimation = "DWLS", model = model, CFIcalc = TRUE, imp_cov = TRUE)

write.table(x = all_bias.fit$results,
            file = "/local1/scratch/pazweifel/all_bias_bip.fit.tsv",
            sep = "\t",
            row.names = FALSE,
            col.names = TRUE,
            quote = FALSE)
            

# 1. Define initial phenotypes

bias_phenotypes <- c("ADHD", "ALC", "AN", "ASD", "BMI", "CAD", "Chrono", 

                     "CRP", "Height", "Neu", "PD", "PTSD", "SCZ", 

                     "Smoking", "T2D", "PartEM", "PartMH", "PartSF")



# Helper function for the model string

build_model <- function(phenos) {

  paste0("

    Liability =~ NA*Clin

    Liability ~~ 1*Liability

    Bias =~ NA*", paste(phenos, collapse = " + "), "

    Bias ~~ 1*Bias

    Liability ~~ 0*Bias

    Clin ~~ 0*Clin

    Clin ~~ ", paste(phenos, collapse = " + ")

  )

}



# 2. Fit the initial model

fit1 <- usermodel(

  bias_all.covstruct, 

  estimation = "DWLS", 

  model = build_model(bias_phenotypes)

)

fit1

# 3. Extract and check the loadings

# Note: lhs is the Latent Factor (Bias), op is the loading (=~), rhs is the phenotype

loadings_df <- subset(fit1$results, lhs == "Bias" & op == "=~")



# 4. Filter phenotypes based on the 0.3 threshold in STD_Genotype

keep_phenos <- loadings_df$rhs[loadings_df$STD_Genotype >= 0.3]



# 5. Fit the pruned model

fit2 <- usermodel(

  bias_all.covstruct, 

  estimation = "DWLS", 

  model = build_model(keep_phenos)

)



fit2

fitlist_yazheng <- list(fit1, fit2)

saveRDS(object = fitlist_yazheng,
        file = "/local1/scratch/pazweifel/bip_bias_fits.R")
        

aaa <- readRDS(file = "/local1/scratch/pazweifel/bip_bias_fits.R")

#calcualte covariance matrix including the selected phenotypes after yazheng analysis, clinical and community
keep_phenos <- keep_phenos[keep_phenos != "ALC"]

covstruct_bias_selected <- here::here("matrices/Clinical_community_selected_external.R")

pop_prev <- c(bip_prevalences["Clin"], bip_prevalences["Comm"], other_prevalences[keep_phenos])

samp_prev <- as.vector(ifelse(is.na(pop_prev), NA, 0.5))

paths <- c(samples[1], samples[2])

for (name in keep_phenos){
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

covstruct.bias_clincom <- dget("matrices/Clinical_community_selected_external.R")

# 1. Define and clean the phenotype list
# Assuming keep_phenos was defined by your previous 0.3 threshold filtering
keep_phenos <- keep_phenos[keep_phenos != "ALC"]

# ---------------------------------------------------------
# 2. Define the NULL Model String
# Commu is NOT loaded by Bias; rG is estimated between Liability and Liab_commu
# ---------------------------------------------------------
model_null <- paste0("
  Liability =~ NA*Clin
  Liability ~~ 1*Liability
  
  Liab_commu =~ NA*Comm
  Liab_commu ~~ 1*Liab_commu
  
  Bias =~ NA*", paste(keep_phenos, collapse = " + "), "
  Bias ~~ 1*Bias
  
  Liability ~~ 0*Bias
  Liab_commu ~~ 0*Bias
  Clin ~~ 0*Clin
  Comm ~~ 0*Comm
  
  # Estimate rG between the two Latent Factors
  Liability ~~ Liab_commu
  
  # Free residual covariances
  Clin ~~ ", paste(keep_phenos, collapse = " + ")
)

# ---------------------------------------------------------
# 3. Define the FULL Model String
# Same as Null, but adding Bias =~ Commu
# ---------------------------------------------------------
model_full <- paste0("
  Liability =~ NA*Clin
  Liability ~~ 1*Liability
  
  Liab_commu =~ NA*Comm
  Liab_commu ~~ 1*Liab_commu

  # Bias now loads onto the pruned list AND Comm
  Bias =~ NA*", paste(keep_phenos, collapse = " + "), " + Comm
  Bias ~~ 1*Bias
  
  Liability ~~ 0*Bias
  Liab_commu ~~ 0*Bias
  Clin ~~ 0*Clin
  Comm ~~ 0*Comm
  
  # Estimate rG between the two Latent Factors
  Liability ~~ Liab_commu
  
  # Free residual covariances
  Clin ~~ ", paste(keep_phenos, collapse = " + ")
)

# ---------------------------------------------------------
# 4. Fit the Models
# ---------------------------------------------------------
fit_null <- usermodel(covstruct.bias_clincom, estimation = "DWLS", model = model_null)
fit_full <- usermodel(covstruct.bias_clincom, estimation = "DWLS", model = model_full)

# ---------------------------------------------------------
# 5. Compare results
# ---------------------------------------------------------

# Extract rG (Liability ~~ Liab_commu) from both
rg_null <- subset(fit_null$results, lhs == "Liability" & op == "~~" & rhs == "Liab_commu")
rg_full <- subset(fit_full$results, lhs == "Liability" & op == "~~" & rhs == "Liab_commu")

# Extract the loading of Bias on Commu from the Full model
loading_full <- subset(fit_full$results, lhs == "Bias" & op == "=~" & rhs == "Comm")

# Combine for comparison table
comparison_table <- data.frame(
  Model = c("Null", "Full", "Full (Loading)"),
  Relationship = c("rG(Liab, Comm)", "rG(Liab, Comm)", "Bias -> Comm"),
  Unstand_Est = c(rg_null$Unstand_Est, rg_full$Unstand_Est, loading_full$Unstand_Est),
  STD_Genotype = c(rg_null$STD_Genotype, rg_full$STD_Genotype, loading_full$STD_Genotype),
  P_Value = c(rg_null$p_value, rg_full$p_value, loading_full$p_value)
)

print(comparison_table)

# Optional: Compare Model Fit Indices
cat("\nModel Fit Comparison:\n")
cat("Null Model AIC:", fit_null$modelfit$AIC, "\n")
cat("Full Model AIC:", fit_full$modelfit$AIC, "\n")

fit_full

Unstand_Est <-c(rg_null$Unstand_Est, rg_full$Unstand_Est, loading_full$Unstand_Est)
STD_Genotype <- c(rg_null$STD_Genotype, rg_full$STD_Genotype, loading_full$STD_Genotype)

Unstand_Est
STD_Genotype

#filter out all traits where a follow up model was fit and which outlier traits were used
follow_up_traits <- follow_up_table_adjusted %>% filter(Unconstrained_paths != "None") %>% subset(select = c(Trait, Unconstrained_paths))

bdII_traits <- c("ADHD", "ASD", "MDD", "Neu", "Smoking")
Comm_traits <- c("ADHD", "AN", "ASD", "MDD", "Neu", "Smoking")

int_ext_traits <- list("BDII" = bdII_traits, "Comm" = Comm_traits)

gwas_by_sub_model <- "
Liability_clin =~ NA*Clin
Liability_other =~ NA*{internal_trait}
Bias =~ NA*{external_trait} + {internal_trait}

Liability_clin ~~ 0*Bias
Liability_other ~~ 0*Bias

{internal_trait} ~~ 0*{internal_trait}
{external_trait} ~~ 0*{external_trait}
{external_trait} ~~ 0*{internal_trait}

Clin ~~ {external_trait}
Clin ~~ 0*{internal_trait}
Clin ~~ 0*Clin
"

#Community
names(Comm_traits) <- Comm_traits

comm_gwas_by_sub.fits <- lapply(Comm_traits, function(external_trait){
    internal_trait <- "Comm"

    usermodel(trait.covstructs[[external_trait]], estimation = "DWLS", std.lv = TRUE, model = str_glue(gwas_by_sub_model))
    })

#Community
names(bdII_traits) <- bdII_traits

bdII_gwas_by_sub.fits <- lapply(bdII_traits, function(external_trait){
    internal_trait <- "BDII"

    usermodel(trait.covstructs[[external_trait]], estimation = "DWLS", std.lv = TRUE, model = str_glue(gwas_by_sub_model))
    })

comm_gwas_by_sub.fits$Neu

community_table <- bind_rows(
    lapply(Comm_traits, function(trait){
        try(
            a <- filter(comm_gwas_by_sub.fits[[trait]]$results, lhs == "Liability_clin" & rhs == "Liability_other")
            )
        try(tibble(rg = a$STD_Genotype, se = a$STD_Genotype_SE))
        }),
    .id = "External_trait")

community_table$Internal_trait <- "Community"

bdii_table <- bind_rows(
    lapply(bdII_traits, function(trait){
        try(
            a <- filter(bdII_gwas_by_sub.fits[[trait]]$results, lhs == "Liability_clin" & rhs == "Liability_other")
            )
        try(tibble(rg = a$STD_Genotype, se = a$STD_Genotype_SE))
        }),
    .id = "External_trait")

bdii_table$Internal_trait <- "BDII"

combined_table <- bind_rows(community_table, bdii_table)

combined_table



#traits: MDD, PTSD, Neu, AN, ASD, ADHD, Smoking
ambivalent_correlation_traits <- c("MDD", "PTSD", "Neu", "AN", "ASD", "ADHD", "Smoking")


two_factor_selected_traits_model <- paste0(
"Bias =~ NA*",paste(ambivalent_correlation_traits, collapse = " + "),"
Liability =~ NA*Clin + ",paste(ambivalent_correlation_traits, collapse = " + "),"
Bias ~~ 0*Liability
Clin ~~ 0*Clin
MDD ~~ a*MDD
a> 0.001")
    


two_factor_selected_traits_model
covstruct_2_factor_bias <- dget("matrices/Clinical_all_external.R")

two_factor_extcorr.fit <- usermodel(covstruct_2_factor_bias, model = two_factor_selected_traits_model, std.lv = TRUE, estimation = "DWLS")

two_factors_semobject <- semPlotModel_GSEM(two_factor_extcorr.fit)

se <- two_factor_extcorr.fit$results$STD_Genotype_SE
est <- two_factor_extcorr.fit$results$STD_Genotype
est.se <- paste0(round(as.numeric(est), 2),"\n (",round(as.numeric(se), 2),")")

pdf("/local1/home/pazweifel/plots/bip_paper/2_factor_extcorr_model.pdf", width = 30, height = 35)
semPaths(two_factors_semobject, layout = "circle", edge.color = "black", sizeMan = 8, sizeLat = 8, nCharNodes = 0, residuals = TRUE, style = "mx",
        edge.label.cex = 0.5,
        esize = 1,
        edgeLabels = est.se)
dev.off()

semPaths(two_factors_semobject, layout = "circle", edge.color = "black", sizeMan = 8, sizeLat = 8, nCharNodes = 0, residuals = TRUE, style = "mx",
        edge.label.cex = 0.5,
        esize = 1,
        edgeLabels = est.se)

two_factor_extcorr.fit

ambivalent_correlation_traits <- c("MDD", "PTSD", "Neu", "AN", "ASD", "ADHD", "Smoking")

#create a path for a covstruct object where you store the covariance structure 
covstruct_ambivalent_correlations_path <- here::here("matrices",paste("Clinical_ambivalent_external", "R", sep = "."))

#extract population prevalence 
pop_prev <- as.vector(all_prevalences[c("Clin", ambivalent_correlation_traits)])
#fill the samp_prev vector depending on the pop_prev vector, if NA use NA and if not use 0.5
samp_prev <- as.vector(ifelse(is.na(pop_prev), NA, 0.5))
#check if the covstruct file already exists to not do the analysis unnecessarily, first time no file there
paths <- c(samples[1])

for (name in ambivalent_correlation_traits){
    if (name == "MDD"){
        paths <- c(paths, "/local1/scratch/pazweifel/jolien_paper_sumstats/munged_sumstats/mdd_munged.sumstats.gz")
        }
    else{
        paths <- c(paths, paste0(other_basepath, name, ".sumstats.gz"))
        }
}

sample_names <- c("Clin", ambivalent_correlation_traits)

if (!file.exists(covstruct_ambivalent_correlations_path)){
    covstruct_ambi_corr_ldsc <- ldsc(
        traits = paths,
        trait.names = sample_names,
        sample.prev = samp_prev,
        population.prev = pop_prev,
        ld = ld,
        wld = weights,
        stand = TRUE
    )
    #store the file in the previously defined path, now the path is, contorl tell how to deparse
    dput(covstruct_ambi_corr_ldsc, covstruct_ambivalent_correlations_path, control = c("all", "digits17"))
}

covstruct_ambivalent_correlations <- dget("matrices/Clinical_ambivalent_external.R")

sumstats_path <- "/local1/scratch/pazweifel/sumstats_ambits/gsem_ready/"
ref <- "/local1/scratch/pazweifel/jolien_paper_sumstats/reference.1000G.maf.0.005.txt"

traits <- c("Clin", "MDD", "PTSD", "Neu", "AN", "ASD", "ADHD", "Smoking")
se.logit <- c(T, T, T, F, T, T, T, F)
linprob <- c(F, F, F, T, F,F,F,F)
ols <- c(F, F, F, F, F, F, F, T)
            
path_clin <- paste0(sumstats_path,"BD_Clin.txt")
path_rest <- paste0(sumstats_path, traits[2:length(traits)],".txt")
all_paths <- c(path_clin, path_rest)

all_paths

#run the sumstats function
dir.create("bipolar_bias", showWarnings = FALSE)

bip_ambivalent_bias_sumstats_path <- "bipolar_bias/clin_ambivalent_bias_combined_sumstats.R"

if (!file.exists(bip_ambivalent_bias_sumstats_path)){
    bias_sumstats <- sumstats(
             files = all_paths,
             ref = ref,
             trait.names = traits,
             se.logit = se.logit,
             OLS = ols, 
             linprob = linprob,
             parallel = TRUE,
             cores = 10)
    
    saveRDS(bias_sumstats, file = bip_ambivalent_bias_sumstats_path) 
    }
     

sumstats_ambivalent_bias <- readRDS("bipolar_bias/clin_ambivalent_bias_combined_sumstats.R")

nrow(sumstats_ambivalent_bias)

bias_gwas <- readRDS("/local1/home/pazweifel/src/bipolar_bias/GWAS_Bias_factor.R")

head(bias_gwas[[1]])

bias_gwas_filtered <- bias_gwas[[1]] %>% filter(error == 0 & warning == 0 & Pval_Estimate < 0.1)

nrow(bias_gwas[[1]])
nrow(bias_gwas_filtered)

library(qqman)

manhattan(bias_gwas_filtered, chr = "CHR", bp = "BP", p = "Pval_Estimate", snp = "SNP")

qq(bias_gwas_filtered$Pval_Estimate)

qq(filter(bias_gwas[[1]], warning == 0 & error == 0)$Pval_Estimate)

path_to_plots <- "/local1/home/pazweifel/plots/bip_paper/gwas_bias_factor"

pdf(paste0(path_to_plots, "_manhattan.pdf"), height = 10, width = 12)
manhattan(bias_gwas_filtered, chr = "CHR", bp = "BP", p = "Pval_Estimate", snp = "SNP")
dev.off()

pdf(paste0(path_to_plots, "_qq.pdf"), height = 10, width = 12)
qq(filter(bias_gwas[[1]], warning == 0 & error == 0)$Pval_Estimate)
dev.off()

bias_gwas_clean <- filter(bias_gwas[[1]], warning == 0 & error == 0)
head(bias_gwas_clean)

#calculate sample size for factor
##Calculate Implied Sample Size for Factor
#restrict to MAF >= 10%

CorrelatedFactors1<-subset(bias_gwas_clean, bias_gwas_clean$MAF >= .1)

N_hat_F1<-mean(1/((2*CorrelatedFactors1$MAF*(1-CorrelatedFactors1$MAF))*CorrelatedFactors1$SE^2))

N_hat_F1

nrow(bias_gwas_clean)

#prepare the sumstats file for further use in gsem
bias_gwas_gsemready <- bias_gwas_clean %>% rename(Beta = est, P = Pval_Estimate, Z = Z_Estimate) %>% mutate(N = N_hat_F1)

head(bias_gwas_gsemready)

data.table::fwrite(x = bias_gwas_gsemready, file = "/local1/scratch/pazweifel/sumstats_ambits/factor_sumstats/Ambivalent_correlation_bias_factor.txt", sep = "\t")

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

### create the new covstruct including the bias factor sumstats
covstruct_ambivalent_bias_factor <- "matrices/ambivalent_bias_factor.R"

path_factor <- "/local1/scratch/pazweifel/sumstats_ambits/munged/Ambivalent_correlation_bias_factor.sumstats.gz"

combined_paths <- c(samples, path_factor)

pop.prev <- c(all_prevalences[bip_names], NA)

samp.prev <- c(ifelse(is.na(pop.prev), NA, 0.5))

if (!file.exists(covstruct_ambivalent_bias_factor)){
    covstruct_ambifac <- ldsc(
        traits = combined_paths,
        trait.names = c(bip_names, "Bias_GWAS"),
        sample.prev = samp.prev,
        population.prev = pop.prev,
        ld = ld,
        wld = weights,
        stand = TRUE
    )
    #store the file in the previously defined path, now the path is, contorl tell how to deparse
    dput(covstruct_ambifac, covstruct_ambivalent_bias_factor, control = c("all", "digits17"))
}

covstruct_ambivalent_factor <- dget("matrices/ambivalent_bias_factor.R")

gbs_traits <- c("BDII", "Comm")
names(gbs_traits) <- gbs_traits

gwas_by_sub_bias.fit <- lapply(gbs_traits, function(trait){
    internal_trait <- trait

    usermodel(covstruct_ambivalent_factor, estimation = "DWLS", model = str_glue(gwas_by_sub_model), std.lv = TRUE, CFIcalc = TRUE)
})

gwas_by_sub_bias.fit[[1]]

gwas_by_sub_bias.fit[[2]]

bias_corrected_liability_correlations <- bind_rows(
    lapply(gbs_traits, function(trait){
        current_df_line <- gwas_by_sub_bias.fit[[trait]]$results %>% filter(lhs == "Liability_clin", rhs == "Liability_other")
        tibble(External_trait = "Bias_Factor", rg = round(current_df_line$STD_Genotype, 3), se = round(as.double(current_df_line$STD_Genotype_SE), 3), Internal_trait = trait)
        }))

bias_corrected_liability_correlations

#load the non-munged summary stats for which the error are missing and the original summary stats
######original_sumstats <- c("daner_adhd_meta_filtered_NA_iPSYCH23_PGC11_sigPCs_woSEX_2ell6sd_EUR_Neff_70.meta.gz",  "iPSYCH-PGC_ASD_Nov2017.gz",  "PGC3_SCZ_wave3.european.autosome.public.v3.vcf.tsv.gz",  "pgcAN2.2019-07.vcf.tsv.gz",  "pgc-panic2019.vcf.tsv.gz")
#####original_path <- "/local1/scratch/pazweifel/sumstats_ambits/original_sumstats/"
####original_paths <- paste0(original_path, original_sumstats)

###bad_sumstats <- c("ADHD.txt", "ASD.txt", "SCZ.txt", "AN.txt", "PD.txt")
##bad_path <- "/local1/scratch/pazweifel/sumstats_ambits/non_munged/"
#bad_paths <- paste0(bad_path, bad_sumstats)

#names(original_paths) <- c("ADHD", "ASD", "SCZ", "AN", "PD")
###names(bad_paths) <- names(original_paths)traits <- names(bad_paths)
##names(traits) <- traits

##sumstats with overhead
#overhead_sumstats <- c("SCZ", "AN", "PD")

#original_datasets <- lapply(traits, function(trait) {
    ######if (trait %in% overhead_sumstats){
     #####   fread(original_paths[[trait]], skip = "CHROM", data.table = FALSE)
   #### }
    ###else{
       # fread(original_paths[[trait]], data.table = FALSE)
    ##}
   # }
#######)

##bad_datasets <- lapply(traits, function(trait) {
#  #  fread(bad_paths[[trait]], data.table = FALSE)
    #})

#for (trait in bad_datasets){
#    print(nrow(trait))#
#    }

#for (trait in bad_datasets){
#    print("SNP" %in% colnames(trait))
#    }

#for (trait in original_datasets){
#    print("SNP" %in% colnames(trait))
#    }

##check how it is named for the ones not calling it SNP
#head(original_datasets$SCZ)
# its "ID", so rename it
#original_datasets_curated <- lapply(traits, function(trait){
#    if (!("SNP" %in% colnames(original_datasets[[trait]]))){
#        rename(original_datasets[[trait]], SNP = ID)
#    }
#    else{
#        original_datasets[[trait]]
#    }
#    }
#)

#for (trait in original_datasets_curated){
#    print("SNP" %in% colnames(trait))
#    }

#now match on the snp column and combined the bad_one with the original one, so that the new df now also contains the error column
#good_sumstats <- lapply(traits, function(trait) {
#    inner_join(x = bad_datasets[[trait]], y = select(original_datasets_curated[[trait]], SNP, SE), by = join_by(SNP))
#    }
#)

#for (trait in good_sumstats){
#    print(sum(is.na(trait$SE)))
#    }

#for (trait in good_sumstats){
#    print(nrow(trait))
#    }

#for (trait in bad_datasets){
#    print(nrow(trait))
#    }

#head(good_sumstats$AN)

#no NA values in the SE column and all variatns could be assigned a standard error. Now save the corrected summary statistics to a GSEM ready folder.

#new_folder_path <- "/local1/scratch/pazweifel/sumstats_ambits/gsem_ready/"
#
#for (trait in traits){
#    write.table(
#        good_sumstats[[trait]],
#        file = paste0(new_folder_path, trait, ".txt"),
#        sep = "\t",
#        col.names = TRUE,
#        row.names = FALSE,
#        quote = FALSE)
#    }

#cleanup of the clinical dataset..

#clinical_bd_bad <- fread("/local1/scratch/pazweifel/heritability_analysis/ldsc_bipolar_comclin/bip2024_eur_clinical_no23andMe_wbeta_use_neff.txt", data.table = FALSE)

#head(clinical_bd_bad)

#clinical_bd_bad_selected <- clinical_bd_bad %>% select(-c(Nca, Nco, N, BETA)) %>% rename(N = NEFF)

#head(clinical_bd_bad_selected)

#clinical_bd_good <- clinical_bd_bad_selected %>% filter(!str_detect(SNP, "^[0-9]"))

#head(clinical_bd_good)
#nrow(clinical_bd_good)
#nrow(clinical_bd_bad)

#clinical_path <- "/local1/scratch/pazweifel/sumstats_ambits/gsem_ready/BD_Clin.txt"

#if (!file.exists(clinical_path)){
#    write.table(x = clinical_bd_good,
#                file = clinical_path,
#                sep = "\t",
#                col.names = TRUE,
#                row.names = FALSE,
#                quote = FALSE)
#}

###calcualte one big covariance matrix involving the clinical bipolar trait and all external ones
#covstruct_bias_selected <- here::here("matrices/Clinical_selected_external.R")

#selected_traits <- c("ADHD", "ALC", "AN", "ASD", "BMI", "CAD", "Chrono", "CRP", "Height", "Neu", "PartEM", "PartMH", "PartSF", "PD", "PTSD", "SCZ", "Smoking", "T2D")

#pop_prev <- c(bip_prevalences["Clin"], other_prevalences[selected_traits])

#samp_prev <- as.vector(ifelse(is.na(pop_prev), NA, 0.5))

#paths <- c(samples[1])
#other_basepath <- "/local1/scratch/pazweifel/sumstats_ambits/munged/"

#names(pop_prev)

#for (name in selected_traits){
#    paths <- c(paths, paste0(other_basepath, name, ".sumstats.gz"))
#}

    
#if (!file.exists(covstruct_bias_selected)){
#    covstruct_bias_selected_ldsc <- ldsc(
#        traits = paths,
#        trait.names = names(pop_prev),
#        sample.prev = samp_prev,
#        population.prev = pop_prev,
#        ld = ld, 
#        wld = weights,
#        stand = TRUE
#        )
#    dput(covstruct_bias_selected_ldsc, covstruct_bias_selected, control = c("all", "digits17"))
#}

#covstruct.bias_selected <- dget("matrices/Clinical_selected_external.R")

#usermodel(covstruc = bias_all.covstruct, estimation = "DWLS", model = model, CFIcalc = TRUE, imp_cov = TRUE)

#sumstats_path <- "/local1/scratch/pazweifel/sumstats_ambits/gsem_ready/"
#ref <- "/local1/scratch/pazweifel/jolien_paper_sumstats/reference.1000G.maf.0.005.txt"

#traits <- c("Clin", "ADHD", "ALC", "AN", "ASD", "BMI", "CAD", "Chrono", "CRP", "Height", "Neu", "PartEM", "PartMH", "PartSF", "PD", "PTSD", "SCZ", "Smoking", "T2D")
#se.logit <- c(T,T,F,T,T,F,T,F,F,F,F,T,T,T,T,T,T,F,T)
#linprob <- c(F,F,T,F,F,F,F,F,F,F,T,F,F,F,F,F,F,F,F)
#ols <- c(F,F,F,F,F,T,F,T,T,T,F,F,F,F,F,F,F,T,F)
#
#path_clin <- paste0(sumstats_path,"BD_Clin.txt")
#path_rest <- paste0(sumstats_path, traits[2:length(traits)],".txt")
#all_paths <- c(path_clin, path_rest)

#run the sumstats function
##dir.create("bipolar_bias", showWarnings = FALSE)

#bip_bias_sumstats_path <- "bipolar_bias/clin_bias_combined_sumstats.R"

#if (!file.exists(bip_bias_sumstats_path)){
#    bias_sumstats <- sumstats(files = all_paths,
#             ref = ref,
#             trait.names = traits,
#             se.logit = se.logit,
#             OLS = ols, 
#             linprob = linprob,
#             parallel = TRUE,
#             cores = 10)
#    
#    saveRDS(bias_sumstats, file = bip_bias_sumstats_path) 
#    }
     

gwas_by_sub_bias.fit[["Comm"]]

traits <- other_names
names(traits) <- traits
traits

base.results.table <- mutate(tibble(base.fit$results), Model = "Base", External_phenotype = "None", p_value = as.character(p_value)) %>% relocate(Model, External_phenotype)

base.modelfit.table <- mutate(tibble(base.fit$modelfit), Model = "Base", External_phenotype = "None") %>% relocate(Model, External_phenotype)

anchor.results.table <- mutate(tibble(anchor.fit$results), Model = "Anchor", External_phenotype = "None", p_value = as.character(p_value)) %>% relocate(Model, External_phenotype)

anchor.modelfit.table <- mutate(tibble(anchor.fit$modelfit), Model = "Anchor", External_phenotype = "None") %>% relocate(Model, External_phenotype)

trait_common.results.table <- bind_rows(
    lapply(traits, function(trait){
        temp <- trait_common.fit[[trait]]$results
        mutate(tibble(temp), Model = "Common", External_phenotype = trait, p_value = as.character(temp$p_value)) %>% relocate(Model, External_phenotype)
        }))

trait_common.modelfit.table <- bind_rows(
    lapply(traits, function(trait){
        temp <- trait_common.fit[[trait]]$modelfit
        mutate(tibble(temp), Model = "Common", External_phenotype = trait) %>% relocate(Model, External_phenotype)
        }))

nrow(trait_common.results.table)

trait_independent.results.table <- bind_rows(
    lapply(traits, function(trait){
        temp <- trait_independent.fit[[trait]]$results
        mutate(tibble(temp), Model = "Independent", External_phenotype = trait, p_value = as.character(temp$p_value)) %>% relocate(Model, External_phenotype)
        }))

trait_independent.modelfit.table <- bind_rows(
    lapply(traits, function(trait){
        temp <- trait_independent.fit[[trait]]$modelfit
        mutate(tibble(temp), Model = "Independent", External_phenotype = trait) %>% relocate(Model, External_phenotype)
        }))

follow_up_names <- names(their_follow_up_models.fit_adjusted)
names(follow_up_names) <- follow_up_names

trait_follow_up.results.table <- bind_rows(
    lapply(traits, function(trait){
        temp <- their_follow_up_models.fit_adjusted[[trait]]$results
        mutate(tibble(temp), Model = "Follow-up", External_phenotype = trait, p_value = as.character(temp$p_value)) %>% relocate(Model, External_phenotype)
        }))

trait_follow_up.modelfit.table <- bind_rows(
    lapply(traits, function(trait){
        temp <- their_follow_up_models.fit_adjusted[[trait]]$modelfit
        mutate(tibble(temp), Model = "Follow-up", External_phenotype = trait) %>% relocate(Model, External_phenotype)
        }))

comm_mask <- sapply(comm_gwas_by_sub.fits, is.null)
comm_gwas_by_sub.fits.filtered <- comm_gwas_by_sub.fits[!comm_mask]

comm_gwasbysub.results.table <- bind_rows(
    lapply(traits, function(trait){
        temp <- comm_gwas_by_sub.fits.filtered[[trait]]$results
        mutate(tibble(temp), p_value = as.character(temp$p_value), Model = "Comm_GWAS_by_sub", External_phenotype = trait) %>% relocate(Model, External_phenotype)
        
        }))

comm_gwasbysub.modelfit.table <- bind_rows(
    lapply(traits, function(trait){

        temp <- comm_gwas_by_sub.fits[[trait]]$modelfit
        mutate(tibble(temp), Model = "Comm_GWAS_by_sub", External_phenotype = trait) %>% relocate(Model, External_phenotype)
        
        }))

bdII_mask <- sapply(bdII_gwas_by_sub.fits, is.null)
bdII_gwas_by_sub.fits.filtered <- bdII_gwas_by_sub.fits[!bdII_mask]

bdII_gwasbysub.results.table <- bind_rows(
    lapply(traits, function(trait){
        temp <- bdII_gwas_by_sub.fits.filtered[[trait]]$results
        mutate(tibble(temp), p_value = as.character(temp$p_value), Model = "BDII_GWAS_by_sub", External_phenotype = trait) %>% relocate(Model, External_phenotype)
        
        }))

bdII_gwasbysub.modelfit.table <- bind_rows(
    lapply(traits, function(trait){

        temp <- bdII_gwas_by_sub.fits.filtered[[trait]]$modelfit
        mutate(tibble(temp), Model = "BDII_GWAS_by_sub", External_phenotype = trait) %>% relocate(Model, External_phenotype)
        
        }))

common_bias.results.table <- mutate(tibble(two_factor_extcorr.fit$results), Model = "Common_bias", External_phenotype = "None", p_value = as.character(p_value)) %>% relocate(Model, External_phenotype)

common_bias.modelfit.table <- mutate(tibble(two_factor_extcorr.fit$modelfit), Model = "Common_bias", External_phenotype = "None") %>% relocate(Model, External_phenotype)

com_bias.traits <- names(gwas_by_sub_bias.fit)
names(com_bias.traits) <- com_bias.traits

common_bias_gwasbysub.results.table <- bind_rows(
    lapply(com_bias.traits, function(trait){
        temp <- gwas_by_sub_bias.fit[[trait]]$results
        mutate(tibble(temp), p_value = as.character(temp$p_value), Model = paste0(trait,"_GWAS_by_sub"), External_phenotype = "Common_Bias") %>% relocate(Model, External_phenotype)
        
        }))

common_bias_gwasbysub.modelfit.table <- bind_rows(
    lapply(com_bias.traits, function(trait){

        temp <- gwas_by_sub_bias.fit[[trait]]$modelfit
        mutate(tibble(temp), Model = paste0(trait,"_GWAS_by_sub"), External_phenotype = "Common_Bias") %>% relocate(Model, External_phenotype)
        
        }))

combined.table.results.goodies <- bind_rows(base.results.table,
                                    anchor.results.table,
                                    trait_common.results.table,
                                    trait_independent.results.table,
                                    trait_follow_up.results.table,
                                    common_bias.results.table)

combined.table.results.badies <- bind_rows(comm_gwasbysub.results.table,
                                           bdII_gwasbysub.results.table,
                                           common_bias_gwasbysub.results.table)
    

combined.table.results <- bind_rows(combined.table.results.goodies, combined.table.results.badies) %>% mutate(op = paste0("`", op))

head(combined.table.results)

write_csv(combined.table.results, file = "/local1/home/pazweifel/plots/bip_paper/combined.results.table.csv")

combined.table.modelfit.goodies <- bind_rows(base.modelfit.table,
                                    anchor.modelfit.table,
                                    trait_common.modelfit.table,
                                    trait_independent.modelfit.table,
                                    trait_follow_up.modelfit.table,
                                    common_bias.modelfit.table)

combined.table.modelfit.badies <- bind_rows(comm_gwasbysub.modelfit.table,
                                           bdII_gwasbysub.modelfit.table,
                                           common_bias_gwasbysub.modelfit.table)
    

combined.table.modelfit <- bind_rows(combined.table.modelfit.goodies, combined.table.modelfit.badies)

write_csv(x = combined.table.modelfit, file = "/local1/home/pazweifel/plots/bip_paper/combined.modelfit.table.csv")

combined.table.modelfit

follow_up_phenos <- filter(combined.table.modelfit, Model == "Follow-up")$External_phenotype

temp <- combined.table.modelfit %>% filter(External_phenotype %in% follow_up_phenos & (Model == "Follow-up" | Model == "Independent" | Model == "Common")) %>% select(-c(chisq, df, p_chisq, AIC)) %>% arrange(External_phenotype) 

paper_table <- temp %>%
  select(Model, External_phenotype, CFI, SRMR) %>%   # drop AIC
  pivot_longer(cols = c(CFI, SRMR),
               names_to = "Metric",
               values_to = "Value") %>%
  mutate(
    Model = factor(Model, levels = c("Common", "Independent", "Follow-up")),
    Metric = factor(Metric, levels = c("CFI", "SRMR")),
    Value = round(Value, 3)
  ) %>%
  pivot_wider(names_from = Model, values_from = Value) %>%
  arrange(External_phenotype, Metric)

paper_table

write_csv(file = "/local1/home/pazweifel/plots/bip_paper/follow_up_table.csv", x = paper_table, col_names = TRUE)

follow_up_table_adjusted

#load the s-ldsc covstruct
sldsc_covstruct <- readRDS("S_LDSC_community_liability_Neu.RData")

#specify the model syntax for a common factor model
model<-"
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
"

#specify the F1 factor variance as the parameter of interest
params<-c("Liability_comm~~Liability_comm")

#use unit variance identification
std.lv=TRUE

#estimate enrichment using the enrich function
#enrich_gwasbysub<-enrich(s_covstruc=sldsc_covstruct,model=model,params=params,std.lv=std.lv)

#specify the model syntax for a common factor model
model<-"
Liability_clin =~ NA*Clin
"

#specify the F1 factor variance as the parameter of interest
params<-c("Liability_clin~~Liability_clin")

#use unit variance identification
std.lv=TRUE

#estimate enrichment using the enrich function
#enrich_clin<-enrich(s_covstruc=sldsc_covstruct,model=model,params=params,std.lv=std.lv)

gwasbysub_gwas <- readRDS("bipolar_bias/GWAS_Comm_liability.R")

##Calculate Implied Sample Size for Factor
#restrict to MAF >= 10%
gwasbysub_gwas_N <- lapply(gwasbysub_gwas, function(gwas){
        subby<-subset(gwas, gwas$MAF >= .1)
        
        N_hat_F1<-round(mean(1/((2*subby$MAF*(1-subby$MAF))*subby$SE^2)), 0)
        mutate(tibble(gwas), N = N_hat_F1) %>% rename(P = Pval_Estimate, Z = Z_Estimate, Beta = est) %>% select(SNP, CHR, BP, MAF, A1, A2, Beta, SE, Z, P, N)
    })

path <- "../plots/bip_paper/gwasbysub_bias_factor_neuro_manhattan.pdf"

if (!file.exists(path)){
    pdf(file = path)
    manhattan(gwasbysub_gwas[[1]], chr = "CHR", snp = "SNP", p = "Pval_Estimate", bp = "BP")
    dev.off()
    }

path <- "../plots/bip_paper/gwasbysub_community_factor_neuro_manhattan.pdf"
if (!file.exists(path)){
    pdf(file = path)
    manhattan(gwasbysub_gwas[[2]], chr = "CHR", snp = "SNP", p = "Pval_Estimate", bp = "BP")
    dev.off()
    }

path <- "../plots/bip_paper/gwasbysub_clinical_factor_neuro_manhattan.pdf"
if (!file.exists(path)){
    pdf(file = path)
    manhattan(gwasbysub_gwas[[3]], chr = "CHR", snp = "SNP", p = "Pval_Estimate", bp = "BP")
    dev.off()
    }

gwasbysub_filepaths <- c("/local1/scratch/pazweifel/sumstats_ambits/factor_sumstats/gwasbysub_neuro_bias.tsv", "/local1/scratch/pazweifel/sumstats_ambits/factor_sumstats/gwasbysub_neuro_community_liab.tsv", "/local1/scratch/pazweifel/sumstats_ambits/factor_sumstats/gwasbysub_neuro_clinical_liab.tsv")

for (i in c(1:length(gwasbysub_gwas_N))){
    if (!file.exists(gwasbysub_filepaths[i])){
        write_tsv(gwasbysub_gwas_N[[i]], file = gwasbysub_filepaths[i])
        }}

head(gwasbysub_gwas_N[[3]])

cluster_annotations <- read.xlsx("bipolar_bias/cluster_annotations.xlsx") %>% mutate(Name = as.character(Cluster.ID)) %>% select(-c(Cluster.ID))

head(cluster_annotations)

annotation_order_data_labels <- c(
  "Miscellaneous",
  "Microglia",
  "Vascular",
  "Fibroblast",
  "Oligodendrocyte precursor",
  "Committed oligodendrocyte precursor",
  "Oligodendrocyte",
  "Bergmann glia",
  "Astrocyte",
  "Ependymal",
  "Choroid plexus",
  "Deep-layer near-projecting",
  "Deep-layer corticothalamic and 6b",
  "Hippocampal CA1-3",
  "Upper-layer intratelencephalic",
  "Deep-layer intratelencephalic",
  "Amygdala excitatory",
  "Hippocampal CA4",
  "Hippocampal dentate gyrus",
  "Medium spiny neuron",
  "Eccentric medium spiny neuron",
  "Splatter",
  "MGE interneuron",
  "LAMP5-LHX6 and Chandelier",
  "CGE interneuron",
  "Upper rhombic lip",
  "Cerebellar inhibitory",
  "Lower rhombic lip",
  "Mammillary body",
  "Thalamic excitatory",
  "Midbrain-derived inhibitory"
)



#load the data
samples <- c("Clinical", "Liability_clinical", "Liability_Community", "Community", "Bias_Neuro")
names(samples) <- samples

pathlist <- c("bipolar_bias/Clinical_sldsc.h2.cell_type_results.txt", "bipolar_bias/gwasbysub_neuro_clinical_liab.h2.cell_type_results.txt", "bipolar_bias/gwasbysub_neuro_community_liab.h2.cell_type_results.txt", "bipolar_bias/Community_sldsc.h2.cell_type_results.txt", "bipolar_bias/gwasbysub_neuro_bias.h2.cell_type_results.txt")
names(pathlist) <- c("Clinical", "Liability_clinical", "Liability_Community", "Community", "Bias_Neuro")

enrich_table <- bind_rows(lapply(samples, function(sample){
    a <- read_table(pathlist[sample])
    a %>% mutate(Phenotype = sample, Name = str_remove(Name, pattern = "Cluster"), Bonfer = Coefficient_P_value*nrow(a), sign = Bonfer < 0.05, BH = p.adjust(Coefficient_P_value, method = "fdr")) %>% left_join(y = cluster_annotations, by = "Name") %>% mutate(Supercluster= factor(Supercluster, levels = annotation_order_data_labels)) %>% arrange(Supercluster, Name)
    })) %>% mutate(Phenotype = factor(Phenotype, levels = unique(Phenotype)), Name = factor(Name, levels = unique(Name))) 

head(enrich_table)

#add a shape column to make the legend more interpretable
x <- rep(x = c(0:5), length.out = 31)
names(x) <- annotation_order_data_labels

enrich_table <- enrich_table %>% mutate(shape = x[Supercluster])

p_cut <- max(enrich_table$Coefficient_P_value[enrich_table$BH <= 0.05], na.rm = TRUE)

library(Polychrome)
cluster_cols <- setNames(
  Polychrome::palette36.colors(31),
  levels(enrich_table$Supercluster)
)

enrich_table$Phenotype <- recode(enrich_table$Phenotype, Liability_clinical = "Liability[BD]", Liability_Community = "Liability[Comm]", Bias_Neuro = "Bias[Neu]")

sldsc_plot <- ggplot(enrich_table, aes(x = Name, y = -log10(Coefficient_P_value))) +
geom_point(aes(color = Supercluster), size = 4, shape = 17) +
scale_color_manual(values = cluster_cols) +
geom_hline(yintercept = -log10(0.05/length(unique(enrich_table$Name))), color = "red", linetype = "dashed") +
geom_hline(yintercept = -log10(p_cut), color = "blue", linetype = "dashed") +
labs(y = expression(-log[10](P)), color = "Category") + 
facet_wrap(vars(Phenotype), nrow = 5, labeller = label_parsed) +
theme_bw(base_size = 25) +
theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
panel.background = element_blank(), axis.line = element_line(colour = "black"), axis.text.x = element_blank(), axis.ticks.x = element_blank(),
strip.background = element_blank(), strip.placement = "outside", axis.title.x = element_blank(), legend.position = "bottom") 

sldsc_plot

ggsave("/local1/home/pazweifel/plots/bip_paper/sldsc_combined.pdf", plot = sldsc_plot, device = "pdf", width = 30, height = 20)

samples <- c("/local1/scratch/pazweifel/heritability_analysis/ldsc_correlation_bipolar/sumstats/bip2024_eur_clinical_no23andMe_wbeta_use_neff.sumstats.gz",
             "/local1/scratch/pazweifel/heritability_analysis/ldsc_correlation_bipolar/sumstats/bip2024_eur_community_no23andMe_wbeta_use_neff.sumstats.gz",
             "/local1/scratch/pazweifel/heritability_analysis/ldsc_correlation_bipolar/sumstats/pgc-bip2021-BDI_cleaned_neff.sumstats.gz",
             "/local1/scratch/pazweifel/heritability_analysis/ldsc_correlation_bipolar/sumstats/pgc-bip2021-BDII_cleaned_neff.sumstats.gz")

### create the new covstruct including the bias factor sumstats
covstruct_ambivalent_bias_factor_2 <- "matrices/ambivalent_bias_factor_clin_com.R"
orig_names <- c("Clin", "Comm")
path_factor <- "/local1/scratch/pazweifel/sumstats_ambits/munged/Ambivalent_correlation_bias_factor.sumstats.gz"

combined_paths <- c(samples[1:2], path_factor)

pop.prev <- c(all_prevalences[orig_names], NA)

samp.prev <- c(ifelse(is.na(pop.prev), NA, 0.5))

if (!file.exists(covstruct_ambivalent_bias_factor_2)){
    covstruct_ambifac <- ldsc(
        traits = combined_paths,
        trait.names = c(orig_names, "Combined_Bias"),
        sample.prev = samp.prev,
        population.prev = pop.prev,
        ld = ld,
        wld = weights,
        stand = TRUE
    )
    #store the file in the previously defined path, now the path is, contorl tell how to deparse
    dput(covstruct_ambifac, covstruct_ambivalent_bias_factor_2, control = c("all", "digits17"))
}

sumstats_path <- "/local1/scratch/pazweifel/sumstats_ambits/gsem_ready/"
factor_path <- "/local1/scratch/pazweifel/sumstats_ambits/factor_sumstats/Ambivalent_correlation_bias_factor.txt"
ref <- "/local1/scratch/pazweifel/jolien_paper_sumstats/reference.1000G.maf.0.005.txt"

traits <- c("Clin", "Comm", "Combined_Bias")
se.logit <- c(T, T, F)
linprob <- c(F, F, F)
ols <- c(F, F, T)
            
path_clin <- paste0(sumstats_path,"BD_Clin.txt")
path_comm <- paste0(sumstats_path, "Comm.txt")
all_paths <- c(path_clin, path_comm, factor_path)

#run the sumstats function
dir.create("bipolar_bias", showWarnings = FALSE)

bip_ambivalent_bias_sumstats_path <- "bipolar_bias/clin_comm_combinedbias_sumstats.R"

if (!file.exists(bip_ambivalent_bias_sumstats_path)){
    bias_sumstats <- sumstats(
             files = all_paths,
             ref = ref,
             trait.names = traits,
             se.logit = se.logit,
             OLS = ols, 
             linprob = linprob,
             parallel = TRUE,
             cores = 3)
    
    saveRDS(bias_sumstats, file = bip_ambivalent_bias_sumstats_path) 
    }
     

combi_factors <- readRDS("bipolar_bias/GWAS_combined_bias_comm_liability.R")

##Calculate Implied Sample Size for Factor
#restrict to MAF >= 10%
gwasbysub_gwas_N <- lapply(combi_factors, function(gwas){
        subby<-subset(gwas, gwas$MAF >= .1)
        
        N_hat_F1<-round(mean(1/((2*subby$MAF*(1-subby$MAF))*subby$SE^2)), 0)
        mutate(tibble(gwas), N = N_hat_F1) %>% rename(P = Pval_Estimate, Z = Z_Estimate, Beta = est) %>% select(SNP, CHR, BP, MAF, A1, A2, Beta, SE, Z, P, N)
    })

path <- "../plots/bip_paper/gwasbysub_bias_factor_common_bias_manhattan.pdf"

if (!file.exists(path)){
    pdf(file = path)
    manhattan(gwasbysub_gwas[[1]], chr = "CHR", snp = "SNP", p = "Pval_Estimate", bp = "BP")
    dev.off()
    }

path <- "../plots/bip_paper/gwasbysub_community_factor_common_bias_manhattan.pdf"
if (!file.exists(path)){
    pdf(file = path)
    manhattan(gwasbysub_gwas[[2]], chr = "CHR", snp = "SNP", p = "Pval_Estimate", bp = "BP")
    dev.off()
    }

path <- "../plots/bip_paper/gwasbysub_clinical_factor_common_bias_manhattan.pdf"
if (!file.exists(path)){
    pdf(file = path)
    manhattan(gwasbysub_gwas[[3]], chr = "CHR", snp = "SNP", p = "Pval_Estimate", bp = "BP")
    dev.off()
    }

gwasbysub_filepaths <- c("/local1/scratch/pazweifel/sumstats_ambits/factor_sumstats/gwasbysub_common_bias_bias.tsv", "/local1/scratch/pazweifel/sumstats_ambits/factor_sumstats/gwasbysub_common_bias_community_liab.tsv", "/local1/scratch/pazweifel/sumstats_ambits/factor_sumstats/gwasbysub_common_bias_clinical_liab.tsv")

for (i in c(1:length(gwasbysub_gwas_N))){
    if (!file.exists(gwasbysub_filepaths[i])){
        write_tsv(gwasbysub_gwas_N[[i]], file = gwasbysub_filepaths[i])
        }}

#load the data
samples <- c("Clinical", "Liability_clinical", "Liability_Community", "Community", "Bias_Common")
names(samples) <- samples

pathlist <- c("bipolar_bias/Clinical_sldsc.h2.cell_type_results.txt", "bipolar_bias/gwasbysub_common_bias_clinical_liab.h2.cell_type_results.txt", "bipolar_bias/gwasbysub_common_bias_community_liab.h2.cell_type_results.txt", "bipolar_bias/Community_sldsc.h2.cell_type_results.txt", "bipolar_bias/gwasbysub_common_bias_bias.h2.cell_type_results.txt")
names(pathlist) <- c("Clinical", "Liability_clinical", "Liability_Community", "Community", "Bias_Common")

enrich_table_common <- bind_rows(lapply(samples, function(sample){
    a <- read_table(pathlist[sample])
    a %>% mutate(Phenotype = sample, Name = str_remove(Name, pattern = "Cluster"), Bonfer = Coefficient_P_value*nrow(a), sign = Bonfer < 0.05, BH = p.adjust(Coefficient_P_value, method = "fdr")) %>% left_join(y = cluster_annotations, by = "Name") %>% mutate(Supercluster= factor(Supercluster, levels = annotation_order_data_labels)) %>% arrange(Supercluster, Name)
    })) %>% mutate(Phenotype = factor(Phenotype, levels = unique(Phenotype)), Name = factor(Name, levels = unique(Name))) 

p_cut_common <- max(enrich_table_common$Coefficient_P_value[enrich_table_common$BH <= 0.05], na.rm = TRUE)

library(Polychrome)
cluster_cols <- setNames(
  Polychrome::palette36.colors(31),
  levels(enrich_table_common$Supercluster)
)

enrich_table_common$Phenotype <- recode(enrich_table_common$Phenotype, Liability_clinical = "Liability[BD]", Liability_Community = "Liability[Comm]", Bias_Common = "Bias[BF]")

sldsc_plot_combias <- ggplot(enrich_table_common, aes(x = Name, y = -log10(Coefficient_P_value))) +
geom_point(aes(color = Supercluster), size = 4, shape = 17) +
scale_color_manual(values = cluster_cols) +
geom_hline(yintercept = -log10(0.05/length(unique(enrich_table_common$Name))), color = "red", linetype = "dashed") +
geom_hline(yintercept = -log10(p_cut_common), color = "blue", linetype = "dashed") +
labs(y = expression(-log[10](P)), color = "Category") + 
facet_wrap(vars(Phenotype), nrow = 5, labeller = label_parsed) +
theme_bw(base_size = 25) +
theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
panel.background = element_blank(), axis.line = element_line(colour = "black"), axis.text.x = element_blank(), axis.ticks.x = element_blank(),
strip.background = element_blank(), strip.placement = "outside", axis.title.x = element_blank(), legend.position = "bottom") 

sldsc_plot_combias

sldsc_plot_combias_wol <- sldsc_plot_combias + theme(legend.position = "None")
sldsc_plot_wol <- sldsc_plot + theme(legend.position = "None")

library(cowplot)

shared_legend <- get_legend(sldsc_plot_combias)

combined_wol <- plot_grid(plotlist = list(sldsc_plot_wol,sldsc_plot_combias_wol), nrow = 2)

combined_wl <- plot_grid(plotlist = list(combined_wol, shared_legend), nrow = 2, rel_heights = c(1, 0.2))

ggsave("/local1/home/pazweifel/plots/bip_paper/combined_enrichment_plot.pdf", combined_wl, device = pdf, height = 30, width = 25)

ggsave("/local1/home/pazweifel/plots/bip_paper/sldsc_combined_commonbias.pdf", plot = sldsc_plot_combias, device = "pdf", width = 30, height = 20)

write_csv(trait.qtrait_adjusted, file = "/local1/home/pazweifel/plots/bip_paper/qtrait_results_adjusted.csv")

follow_up_table_adjusted

colnames(follow_up_table_adjusted)

follow_up_table_short <- trait.qtrait_adjusted %>% select(c(Trait, rGF1Trait_CPM, SErGF1Trait_CPM, QTrait_CPM, df_CPM, lSRMR_CPM, rGF1Trait_FUM, SErGF1Trait_FUM, QTrait_FUM, df_FUM, lSRMR_FUM, Unconstrained_paths))
follow_up_table_short

follow_up_table_short_paper <- follow_up_table_short %>% mutate(rg_common_combined = paste0(round(rGF1Trait_CPM, 2), " (", round(SErGF1Trait_CPM, 2), ")"), Qtrait_common_combined = paste0(round(QTrait_CPM, 2), " (", df_CPM, ")")) %>%
mutate(rg_fu_combined = paste0(round(rGF1Trait_FUM, 2), " (", round(SErGF1Trait_FUM, 2), ")"), Qtrait_fu_combined = paste0(round(QTrait_FUM, 2), " (", df_FUM, ")"), lSRMR_CPM = round(lSRMR_CPM,2), lSRMR_FUM = round(lSRMR_FUM, 2)) %>%
select(c(Trait, rg_common_combined, Qtrait_common_combined, lSRMR_CPM, rg_fu_combined, Qtrait_fu_combined, lSRMR_FUM, Unconstrained_paths)) %>% arrange(Unconstrained_paths) %>%
mutate(lSRMR_FUM = ifelse(is.na(lSRMR_FUM), "-", lSRMR_FUM), rg_fu_combined = ifelse(rg_fu_combined == "NA (NA)", "-", rg_fu_combined),  Qtrait_fu_combined = ifelse(Qtrait_fu_combined == "NA (NA)", "-", Qtrait_fu_combined))

follow_up_table_short_paper

knitr::kable(follow_up_table_short_paper, format = "pipe")

write.table(x = follow_up_table_short_paper, file = "follow_up_table_bipolar.tsv", sep = "\t", row.names = FALSE)

unique(enrich_table$Phenotype)
unique(enrich_table_common$Phenotype)

write_csv(enrich_table, "/local1/home/pazweifel/plots/bip_paper/enrichments_neuroticism.csv")
write_csv(enrich_table_common, "/local1/home/pazweifel/plots/bip_paper/enrichments_bias_facotr.csv")

#need Clinical, Community, Liability Community Neu, Liability Community common bias, ADHD, MDD, PTSD, Clinical MDD and Help

path_original_base <- "/local1/scratch/pazweifel/sumstats_ambits/munged/"
path_factor_base <- "/local1/scratch/pazweifel/sumstats_ambits/factor_sumstats_munged/"

reg_sumstats <- c("Clin", "Comm", "ADHD", "MDD", "PTSD")
factor_sumstats <- c("gwasbysub_common_bias_community_liab", "gwasbysub_neuro_community_liab")

samples <- c("/local1/scratch/pazweifel/heritability_analysis/ldsc_correlation_bipolar/sumstats/bip2024_eur_clinical_no23andMe_wbeta_use_neff.sumstats.gz",
             "/local1/scratch/pazweifel/heritability_analysis/ldsc_correlation_bipolar/sumstats/bip2024_eur_community_no23andMe_wbeta_use_neff.sumstats.gz",
             "/local1/scratch/pazweifel/heritability_analysis/ldsc_correlation_bipolar/sumstats/pgc-bip2021-BDI_cleaned_neff.sumstats.gz",
             "/local1/scratch/pazweifel/heritability_analysis/ldsc_correlation_bipolar/sumstats/pgc-bip2021-BDII_cleaned_neff.sumstats.gz")

samples_correlations <- c(samples[1:2])
for (phenotype in reg_sumstats[3:length(reg_sumstats)]){
    if (phenotype == "MDD"){
            samples_correlations <- c(samples_correlations, "/local1/scratch/pazweifel/jolien_paper_sumstats/munged_sumstats/mdd_munged.sumstats.gz")
        } else {
            samples_correlations <- c(samples_correlations, paste0(other_basepath,phenotype,".sumstats.gz"))
        }}

samples_factors <- paste0(path_factor_base, factor_sumstats, ".sumstats.gz")

samples_correlations_all <- c(samples_correlations, samples_factors)

covstruct_path_correlations <- "matrices/correlations.R"

trait_names <- c(reg_sumstats, "Comm_cBias_corrected", "Comm_Neu_corrected")
pop.prev <- c(all_prevalences[reg_sumstats], NA, NA)
samp.prev <- ifelse(is.na(pop.prev), NA, 0.5)

if (!file.exists(covstruct_path_correlations)){
    covstruct <- ldsc(
        traits = samples_correlations_all,
        sample.prev = samp.prev,
        population.prev = pop.prev,
        ld = ld,
        wld = weights,
        stand = TRUE,
        trait.names = trait_names)
    dput(covstruct, covstruct_path_correlations, control = c("all", "digits17"))
    }   

correlation_covstruct <- dget("matrices/correlations.R")

sstand <- correlation_covstruct$S_Stand
vstand <- correlation_covstruct$V_Stand

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
expanded_df$values <- round(sstand_vectorized, 4)
expanded_df$se <- round(error_vectorized, 4)

Comm_phenos <- c("Comm", "Comm_cBias_corrected", "Comm_Neu_corrected")

filtered_corr_df <- expanded_df %>% filter(y %in% Comm_phenos & !x %in% Comm_phenos)

comm_corr_plot <- ggplot(filtered_corr_df, aes(x = x, y = values, fill = y, group = y)) +
geom_col(color = "black", position = position_dodge()) +
geom_errorbar(aes(ymin = values - se, ymax = values + se), position = position_dodge(width = 0.9), width = 0.5) +
theme_bw(base_size = 25) +
scale_fill_manual(values = c("grey", "darkorchid", "orange")) +
theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
panel.background = element_blank(), axis.line = element_line(colour = "black"), axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "bottom") +
labs(y = expression(r[g]), x = "Phenotype", fill = "Community dataset") +
scale_y_continuous(expand = expansion(mult = c(0,0.05))) 

comm_corr_plot

ggsave("/local1/home/pazweifel/plots/bip_paper/comm_phenos_correlations.pdf", plot = comm_corr_plot, device = "pdf", width = 15, height = 10)

write_csv(filtered_corr_df, "/local1/home/pazweifel/plots/bip_paper/community_cleaned_correlations.csv")

combined_table
bias_corrected_liability_correlations

base_covstruct <- dget("matrices/BD.R")

sstand <- base_covstruct$S_Stand
vstand <- base_covstruct$V_Stand

#extract the standard error
vstand_diag <- sapply(c(1:10), function(index) vstand[index, index])
                      
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
expanded_df$values <- round(sstand_vectorized, 3)
expanded_df$se <- round(error_vectorized, 3)
observed <- filter(expanded_df, x == "Clin"  &( y == "Comm" | y== "BDII")) %>% rename(Internal_trait = y, rg = values) %>% mutate(External_trait = "Observed") %>% select(-c(x))

combined_table <- combined_table %>% mutate(se = as.double(se))
gwasbysub_liability_correlations <- bind_rows(combined_table, bias_corrected_liability_correlations, observed) %>% mutate(Internal_trait = recode(Internal_trait, "Comm" = "Community")) %>% arrange(desc(rg)) %>%
mutate(External_trait = factor(External_trait, levels = unique(External_trait))) %>% arrange(Internal_trait, desc(rg)) %>% mutate(Internal_trait = factor(Internal_trait, levels = unique(Internal_trait)))


gwasbysub_liability_correlations 

correlations_gwasbysub_plot <- ggplot(gwasbysub_liability_correlations, aes(x = Internal_trait, y = rg, fill = External_trait, group = External_trait)) +
geom_col(color = "black", position = position_dodge()) +
geom_errorbar(aes(ymin = rg - se, ymax = rg + se), position = position_dodge(width = 0.9), width = 0.5) +
theme_bw(base_size = 25) +
theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
panel.background = element_blank(), axis.line = element_line(colour = "black"), axis.text.x = element_text(angle = 45, hjust = 1)) +
labs(y = expression(r[g]), x = "Internal outlier", fill = "External Phenotype")

correlations_gwasbysub_plot

ggsave("/local1/home/pazweifel/plots/bip_paper/gwasbysub_correlations_liability.pdf", correlations_gwasbysub_plot, device = "pdf", width = 10, height = 8)
