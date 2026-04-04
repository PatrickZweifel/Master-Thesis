#!/bin/bash
#get all the filenames and put them in a list
files1=($(ls "/local1/scratch/pazweifel/jolien_paper_sumstats/cleaned_sumstats"))

#initialize and fill the prevalence associative arrays
declare -A SNP
declare -A A1
declare -A A2
declare -A sumstats
declare -A center
declare -A info
declare -A P


SNP[asd]="SNP"
SNP[adhd]="SNP"
SNP[bip]="SNP"
SNP[scz]="ID"
SNP[mdd]="ID"

A1[asd]="A1"
A1[adhd]="A1"
A1[bip]="A1"
A1[scz]="A1"
A1[mdd]="EA"

A2[asd]="A2"
A2[adhd]="A2"
A2[bip]="A2"
A2[scz]="A2"
A2[mdd]="NEA"

sumstats[asd]="OR"
sumstats[adhd]="OR"
sumstats[bip]="OR"
sumstats[scz]="BETA"
sumstats[mdd]="BETA"

center[asd]=1
center[adhd]=1
center[bip]=1
center[scz]=0
center[mdd]=0

info[asd]="INFO"
info[adhd]="INFO"
info[bip]="INFO"
info[scz]="IMPINFO"
info[mdd]="IMPINFO"

P[asd]="P"
P[adhd]="P"
P[bip]="P"
P[scz]="PVAL"
P[mdd]="PVAL"
#path to files
sumstatspath="/local1/scratch/pazweifel/jolien_paper_sumstats/cleaned_sumstats/"
for file1 in "${files1[@]}"; do
        #remove the suffix so that you only have the filename
        filename1=${file1%"_eur_neff.gz"}
        #do the comparison for the two sumstat files
        /local1/home/pazweifel/bin/ldsc/munge_sumstats.py \
                --sumstats "${sumstatspath}${file1}" \
                --a1 ${A1[${filename1}]} \
                --a2 ${A2[${filename1}]} \
		--snp ${SNP[${filename1}]} \
                --p ${P[${filename1}]} \
                --signed-sumstats ${sumstats[${filename1}]},${center[${filename1}]} \
		--info ${info[${filename1}]} \
		--N-col N_EFF \
		--ignore N \
                --out "/local1/scratch/pazweifel/jolien_paper_sumstats/munged_sumstats/${filename1}_munged" &
done
