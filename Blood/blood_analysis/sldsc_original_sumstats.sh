#!/bin/bash

for file in /local1/scratch/pazweifel/sumstats_blood/munged_sumstats/*.sumstats.gz; do
 #extract filename
        filename=$(basename "$file" .sumstats.gz)

        /local1/home/pazweifel/bin/ldsc/ldsc.py \
                --h2-cts "$file" \
                --ref-ld-chr /local1/hdata/REF/baseline_v2.2/Variant_effects_baselineLD_v2.2_baselineLD. \
                --w-ld-chr /local1/scratch/pazweifel/heritability_analysis/1000G_Phase3_weights_hm3_no_MHC/weights.hm3_noMHC. \
                --ref-ld-chr-cts /local1/scratch/pazweifel/single_cell_1k1k/single_cell_1k1k.local.ldcts \
                --out /local1/home/pazweifel/src/blood_analysis/${filename}.h2 &
done
