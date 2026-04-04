#!/bin/bash


/local1/home/pazweifel/bin/ldsc/ldsc.py \
  --h2-cts /local1/scratch/pazweifel/heritability_analysis/ldsc_correlation_bipolar/sumstats/bip2024_eur_community_no23andMe_wbeta_use_neff.sumstats.gz \
  --ref-ld-chr /local1/hdata/REF/baseline_v2.2/Variant_effects_baselineLD_v2.2_baselineLD. \
  --w-ld-chr /local1/scratch/pazweifel/heritability_analysis/1000G_Phase3_weights_hm3_no_MHC/weights.hm3_noMHC. \
  --ref-ld-chr-cts /local1/scratch/pazweifel/single_cell_Duncan_og/single_cell_Duncan_og.local.ldcts \
  --out /local1/home/pazweifel/src/bipolar_bias/Community_sldsc.h2 \
  --samp-prev 0.5 \
  --pop-prev 0.02
