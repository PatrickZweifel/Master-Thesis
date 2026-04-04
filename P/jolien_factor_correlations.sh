#!/bin/bash

LDSC=/local1/home/pazweifel/bin/ldsc
TARGET=/local1/home/pazweifel/src/common_factor_jolien/mental_common_factor_wN.sumstats.gz
OTHER_DIR=/local1/scratch/pazweifel/sumstats_covariates/combined_munged_sumstats
REF_LD=/local1/hdata/REF/eur_ref_ld_chr/
W_LD=/local1/hdata/REF/eur_w_ld_chr/

OUTDIR=/local1/home/pazweifel/src/common_factor_jolien/rg_results
OUTPREFIX=${OUTDIR}/mental_common_factor_wN_vs_all

mkdir -p "${OUTDIR}"

target_base=$(basename "${TARGET}")
RG_LIST="${TARGET}"

for f in "${OTHER_DIR}"/*.sumstats.gz; do
    [[ "$(basename "$f")" == "$target_base" ]] && continue
    RG_LIST="${RG_LIST},${f}"
done

python "${LDSC}/ldsc.py" \
    --rg "${RG_LIST}" \
    --ref-ld-chr "${REF_LD}" \
    --w-ld-chr "${W_LD}" \
    --out "${OUTPREFIX}"

awk '
BEGIN {capture=0; OFS="\t"; print "trait\trg\tse\tz\tp"}
/Summary of Genetic Correlation Results/ {capture=1; next}
capture && /^p1/ {next}
capture && NF==0 {exit}
capture {
    p2=$2; rg=$3; se=$4; z=$5; p=$6
    gsub(".*/","",p2)
    gsub(".sumstats.gz","",p2)
    print p2, rg, se, z, p
}
' "${OUTPREFIX}.log" > "${OUTDIR}/rg_summary.tsv"
