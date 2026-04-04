#!/bin/bash

sumstats="/local1/home/pazweifel/src/common_factor_jolien/mental_common_factor_filtered.tsv.gz"
genotype=/local1/scratch/caina/ukb/genos/ukb_imp_v3_chr#.unrelatedbritishqced.maf001geno9.biallelic
ld=/local1/hdata/REF/eur_ref_ld_chr/#
destination_folder="/local1/scratch/pazweifel/PRS/jolien_common_factor/score"

/local1/home/pazweifel/bin/PRSice_linux \
	--no-regress \
	--out "$destination_folder" \
	--base "$sumstats" \
	--target $genotype \
	--snp SNP \
	--a1 A1 \
	--a2 A2 \
	--bp BP \
	--pvalue Pval_Estimate \
	--stat est \
	--beta \
	--base-maf MAF:0.05 \
	--thread 16 \
	--bar-levels 5e-8,5e-5,1e-5,1e-4,0.001,0.01,0.05,0.1,0.2,0.3,0.4,0.5,1 \
	--fastscore
