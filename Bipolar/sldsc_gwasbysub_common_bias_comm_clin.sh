#!/bin/bash

for file in /local1/scratch/pazweifel/sumstats_ambits/factor_sumstats/gwasbysub_common_bias*; do
 #extract filename
        filename=$(basename "$file" .tsv)

        /local1/home/pazweifel/bin/ldsc/munge_sumstats.py \
                --sumstats "$file" \
                --snp SNP \
                --a1 A1 \
                --a2 A2 \
                --p P \
                --signed-sumstats Beta,0 \
                --N-col N \
                --out "/local1/scratch/pazweifel/sumstats_ambits/factor_sumstats_munged/$filename"

        /local1/home/pazweifel/bin/ldsc/ldsc.py \
                --h2-cts "/local1/scratch/pazweifel/sumstats_ambits/factor_sumstats_munged/${filename}.sumstats.gz" \
                --ref-ld-chr /local1/hdata/REF/baseline_v2.2/Variant_effects_baselineLD_v2.2_baselineLD. \
                --w-ld-chr /local1/scratch/pazweifel/heritability_analysis/1000G_Phase3_weights_hm3_no_MHC/weights.hm3_noMHC. \
                --ref-ld-chr-cts /local1/scratch/pazweifel/single_cell_Duncan_og/single_cell_Duncan_og.local.ldcts \
                --out /local1/home/pazweifel/src/bipolar_bias/${filename}.h2 &
done
