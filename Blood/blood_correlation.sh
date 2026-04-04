#!/bin/bash
#get all the filenames and put them in a list
files1=($(ls "/local1/scratch/pazweifel/heritability_analysis/ldsc_correlation_blood/sumstats"))
files2=("${files1[@]}")
#define max number of jobs
max_jobs=10
#path to files
sumstatspath="/local1/scratch/pazweifel/heritability_analysis/ldsc_correlation_blood/sumstats"
for file1 in "${files1[@]}"; do
        #remove the current element from list 2 so that you dont make unnecessary comparisons
        files2=( "${files2[@]/$file1}" )
        #remove the suffix so that you only have the filename
        filename1=${file1%".sumstats.gz"}
        for file2 in "${files2[@]}"; do
                #extract the filename as for file1
                filename2=${file2%".sumstats.gz"}
                #do the comparison for the two sumstat files
                /local1/home/pazweifel/bin/ldsc/ldsc.py \
                --rg "${sumstatspath}/$file1,${sumstatspath}/$file2" \
		--ref-ld-chr "/local1/hdata/REF/baseline_v2.2/Variant_effects_baselineLD_v2.2_baselineLD." \
                --w-ld-chr "/local1/scratch/pazweifel/heritability_analysis/1000G_Phase3_weights_hm3_no_MHC/weights.hm3_noMHC." \
                --out "/local1/scratch/pazweifel/heritability_analysis/ldsc_correlation_blood/${filename1}_$filename2"

	done
done

