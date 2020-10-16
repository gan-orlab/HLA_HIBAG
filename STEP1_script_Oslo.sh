#!/bin/bash

#Extract SNPs hard calls in reference panel from VCF; add sex, pheno; remove duplicates

FILE=$1
REF=$2
ANCESTRY=$3


plink --const-fid --vcf  vcfFiles/$FILE.vcf.gz --keep MERGED_unrelated.fam --make-bed --out bfile/$FILE

awk '{$2=$1":"$4; print $0}' bfile/$FILE.bim > a
mv a bfile/$FILE.bim

awk '{print $2}' bfile/$FILE.bim | sort | uniq -d  > txt_data/dup/$FILE.dup

plink --bfile bfile/$FILE --exclude txt_data/dup/$FILE.dup --make-bed --out bfile/unique_$FILE

plink --bfile bfile/unique_$FILE --qual-scores vcfFiles/$FILE.info 7 1 1 --qual-threshold 0.8 --make-bed --out bfile/$FILE'_hardcall'

cut -f 3-6 bfile/$FILE'_hardcall'.bim > bfile/temp.bim
cut -f 1-2 bfile/$FILE'_hardcall'.bim | cut -d ":" -f 1-2 > bfile/temp2.bim
paste -d "\t" bfile/temp2.bim bfile/temp.bim > bfile/$FILE'_hardcall'.bim

rm bfile/temp*

plink --bfile bfile/$FILE'_hardcall' --extract RData/csv/$REF/HLA_snp.csv --update-sex txt_data/sex/$FILE.sex --pheno txt_data/pheno/$FILE.pheno --make-bed --out pre_FRED

Rscript FRED.R RData/csv/$REF/unique_HLA.csv pre_FRED.bim


mv a pre_FRED.bim

plink --bfile pre_FRED --extract RData/csv/$REF/HLA_snp.csv --make-bed --out bfile/final_$FILE'_'$ANCESTRY

rm pre_FRED*

