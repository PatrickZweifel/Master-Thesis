library(fmsb)
library(data.table)

set.seed(42)

#### -------------------------
#### 1. Load data
#### -------------------------
prs_file <- "/local1/scratch/pazweifel/PRS/MDD_Adams_25/score.all_score"

df_prs <- fread(prs_file, sep=" ")
df_y_multi <- fread("/local1/hdata/ukbb/phenos/allvariable.version2.txt")

df_prs[, FID := as.character(FID)]
df_prs[, IID := as.character(IID)]
df_y_multi[, FID := as.character(FID)]
df_y_multi[, IID := as.character(IID)]

setkey(df_prs, FID, IID)
setkey(df_y_multi, FID, IID)

prs_cols <- setdiff(colnames(df_prs), c("FID","IID"))
cov_cols <- c(paste0("PC",1:20), "array")
y_cols  <- setdiff(colnames(df_y_multi), c("FID","IID",cov_cols))

#### -------------------------
#### 2. Merge only necessary rows
#### -------------------------

# Keep only individuals that appear in PRS file
df_y_multi <- df_y_multi[df_prs[, .(FID, IID)], nomatch=0]
# Merge all into single table
df_all <- merge(df_prs, df_y_multi, by=c("FID","IID"))
# colnames(df_all)

#### -------------------------
#### 3. Function to detect binary phenotype
#### -------------------------
is_binary <- function(x) {
  ux <- unique(na.omit(x))
  length(ux) == 2 && all(sort(ux) %in% c(0,1))
}

#### -------------------------
#### 4. Loop over phenotypes × PRS and compute ΔR²
#### -------------------------

results <- list()

for (y in y_cols) {
  for (prs in prs_cols) {
    dt <- df_all[, c("FID","IID", y, prs, cov_cols), with=FALSE]
    dt <- na.omit(dt)

    # Skip degenerate cases
    if (nrow(dt) < 100) next 
    if (all(dt[[y]] == dt[[y]][1])) next  # constant outcome
    
    # Model formulas
  cov_formula <- paste0("`", cov_cols, "`", collapse=" + ")
  form1 <- as.formula(paste0("`", y, "` ~ ", cov_formula))
  form2 <- as.formula(paste0("`", y, "` ~ ", cov_formula, " + `", prs, "`"))

    if (is_binary(dt[[y]])) {
      #### -------------------------
      #### Logistic regression + Nagelkerke R²
      #### -------------------------
      fit1 <- glm(form1, data=dt, family=binomial)
      fit2 <- glm(form2, data=dt, family=binomial)

      r2_1 <- NagelkerkeR2(fit1)$R2
      r2_2 <- NagelkerkeR2(fit2)$R2
      delta_r2 <- r2_2 - r2_1

      results[[length(results)+1]] <- data.table(
        phenotype = y,
        prs = prs,
        type = "binary",
        R2_null = r2_1,
        R2_full = r2_2,
        delta_R2 = delta_r2,
        N = nrow(dt)
      )

    } else {
      #### -------------------------
      #### Linear regression + ordinary R²
      #### -------------------------
      fit1 <- lm(form1, data=dt)
      fit2 <- lm(form2, data=dt)

      r2_1 <- summary(fit1)$r.squared
      r2_2 <- summary(fit2)$r.squared
      delta_r2 <- r2_2 - r2_1

      results[[length(results)+1]] <- data.table(
        phenotype = y,
        prs = prs,
        type = "continuous",
        R2_null = r2_1,
        R2_full = r2_2,
        delta_R2 = delta_r2,
        N = nrow(dt)
      )
    }
  }
}

#### -------------------------
#### 5. Final results table
#### -------------------------
results_dt <- rbindlist(results)
results_dt
write.csv(results_dt, file = "/local1/scratch/pazweifel/PRS/MDD_Adams_25/R_squared_ukbb.csv")
print("Done!")