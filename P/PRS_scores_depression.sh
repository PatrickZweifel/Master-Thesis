#!/bin/bash

sumstats="/local1/scratch/pazweifel/jolien_paper_sumstats/cleaned_sumstats/mdd_eur_neff_2.gz"
genotype=/local1/scratch/caina/ukb/genos/ukb_imp_v3_chr#.unrelatedbritishqced.maf001geno9.biallelic
ld=/local1/hdata/REF/eur_ref_ld_chr/#
destination_folder="/local1/scratch/pazweifel/PRS/MDD_Adams_25/score"

/local1/home/pazweifel/bin/PRSice_linux \
        --no-regress \
        --out "$destination_folder" \
        --base "$sumstats" \
        --target $genotype \
        --snp ID \
        --a1 EA \
        --a2 NEA \
        --bp POS \
        --pvalue PVAL \
        --stat BETA \
        --beta \
        --thread 16 \
        --bar-levels 5e-8,5e-5,1e-5,1e-4,0.001,0.01,0.05,0.1,0.2,0.3,0.4,0.5,1 \
        --fastscore
