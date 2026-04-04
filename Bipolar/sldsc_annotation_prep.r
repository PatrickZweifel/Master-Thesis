#load the required packages

library(lavaanPlot)
library(lavaan)
library(semPlot)

#function to convert your gsem object to an object structure the semPaths function recognizes

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

#naming should be compatible with the main analysis file

for (target in other_names){
    #define the two sem objects
    current_fit_common <- semPlotModel_GSEM(gsem.object = target_common.fits[[target]])
    current_fit_independent <- semPlotModel_GSEM(gsem.object = target_independent.fits[[target]])
    
    #extract the std estimate and the corresponding se to plot them both, for common and independent model
    se <- target_common.fits[[trait]]$results$STD_Genotype_SE
    est <- target_common.fits[[trait]]$results$STD_Genotype
    est.se <- paste0(round(as.numeric(est), 2),"\n (",round(as.numeric(se), 2),")")

    se_in <- target_independent.fits[[target]]$results$STD_Genotype_SE
    est_in <- target_independent.fits[[target]]$results$STD_Genotype
    est.se_in <- paste0(round(as.numeric(est_in), 2),"\n (",round(as.numeric(se_in), 2),")")
    
    #save the common model as an svg (or whatever type you need) object, remove the last label as overlapping and known as factor set to variance of 1
    #use the custom labels including the se
    svg(paste0("/path/to/plot/folder/",trait,"_common_model.svg"), width = 11, height = 12)
    semPaths(current_fit_common, layout = "tree", edge.color = "black", sizeMan = 12, sizeLat = 12, nCharNodes = 0, residuals = TRUE, style = "mx",
        edge.label.cex = 1.2,
        esize = 2,
        edgeLabels = est.se[1:(length(est.se)-1)])
    dev.off()
    
    #save the independent model as an svg object
    svg(paste0("/path/to/plot/folder/",trait,"_independent_model.svg"), width = 12, height = 11)
    semPaths(current_fit_independent, layout = "circle", edge.color = "black", sizeMan = 12, sizeLat = 12, nCharNodes = 0, residuals = TRUE, style = "mx",
        edge.label.cex = 1.2,
        rotation = 1,
        esize = 2,
        edgeLabels = est.se_in)
    dev.off()
}