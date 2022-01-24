#!/bin/bash

module load gcc/7.3.0 r-bundle-bioconductor/3.9


FILE=bfile/final_$1_$2

SUFFIX=$1_$2

REF=$3
core=$4
allele=$2

mkdir $SUFFIX
mkdir log/$SUFFIX

Rscript ~/runs/eyu8/data/HLA_typing/HIBAG/script/HIBAG_script_Psychchip.R $allele $FILE $SUFFIX $REF $core

rm -r csv/$SUFFIX

mv $SUFFIX csv/

sed -i -- 's/"//g' csv/$SUFFIX/*
