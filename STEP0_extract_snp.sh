#!/bin/bash

RDATA=$1
REF=$2

mkdir tmp
Rscript script/extract_snp.R $RDATA

cat tmp/*position* > tmp/HLA_snp_position.csv
cat tmp/*allele* > tmp/HLA_snp_allele.csv
paste tmp/HLA_snp_position.csv tmp/HLA_snp_allele.csv | sed 's/"//g'|tr '/' '\t' > tmp/dup.txt

sort -u tmp/dup.txt > tmp/unique_HLA.csv

awk '{print "6:"$1}' tmp/unique_HLA.csv > tmp/HLA_snp.csv

mv tmp/ RData/csv/$REF/
