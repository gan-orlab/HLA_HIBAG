#!/bin/bash

#Extract SNPs hard calls in reference panel from VCF; add sex, pheno; remove duplicates

FILE=$1
REF=$2
ANCESTRY=$3


plink --bfile bfile/$FILE --snps-only just-acgt --make-bed --out bfile/common_$FILE

plink --bfile bfile/common_$FILE --extract RData/csv/$REF/HLA_snp.csv --update-sex txt_data/sex/$FILE.sex --pheno txt_data/pheno/$FILE.pheno --make-bed --out pre_FRED

Rscript FRED.R RData/csv/$REF/unique_HLA.csv pre_FRED.bim

mv a pre_FRED.bim

plink --bfile pre_FRED --extract RData/csv/$REF/HLA_snp.csv --make-bed --out bfile/final_$FILE'_'$ANCESTRY

rm pre_FRED*

    
