#!/bin/bash

#Extract SNPs hard calls in reference panel from VCF; add sex, pheno; remove duplicates

FILE=$1
REF=$2
ANCESTRY=$3


awk '{print 0" "$1"_"$2" "$6}' ~/runs/eyu8/QC_NIH/QC_$FILE/*/data_here/FILTERED.fam > txt_data/pheno/$FILE.pheno
awk '{print 0" "$1"_"$2" "$5}' ~/runs/eyu8/QC_NIH/QC_$FILE/*/data_here/FILTERED.fam > txt_data/sex/$FILE.sex

plink --const-fid --vcf  vcfFiles/$FILE.vcf.gz  --keep MERGED_unrelated.fam --make-bed --out bfile/$FILE

awk '{print $2}' bfile/$FILE.bim | sort | uniq -d  > txt_data/dup/$FILE.dup

plink --bfile bfile/$FILE --exclude txt_data/dup/$FILE.dup --make-bed --out bfile/unique_$FILE

plink --bfile bfile/unique_$FILE --qual-scores vcfFiles/$FILE.info 7 1 1 --qual-threshold 0.8 --make-bed --out bfile/$FILE'_hardcall'

cut -f 3-6 bfile/$FILE'_hardcall'.bim > bfile/temp.bim
cut -f 1-2 bfile/$FILE'_hardcall'.bim | cut -d ":" -f 1-2 > bfile/temp2.bim
paste -d "\t" bfile/temp2.bim bfile/temp.bim > bfile/$FILE'_hardcall'.bim

rm bfile/temp*

plink --bfile bfile/$FILE'_hardcall' --extract RData/csv/$REF/HLA_snp.csv --update-sex txt_data/sex/$FILE.sex --pheno txt_data/pheno/$FILE.pheno --make-bed --out pre_exc

Rscript exc.R RData/csv/$REF/unique_HLA.csv pre_exc.bim

mv a pre_exc.bim

plink --bfile pre_exc --extract RData/csv/$REF/HLA_snp.csv --make-bed --out bfile/final_$FILE'_'$ANCESTRY

rm pre_exc*

    
