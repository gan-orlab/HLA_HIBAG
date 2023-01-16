#!/bin/bash

module load singularity

FILE=bfile/final_$1_$2

SUFFIX=$1_$2

REF=$3
core=$4
allele=$2

mkdir $SUFFIX
mkdir log/$SUFFIX

singularity exec -B /lustre03,/project,/scratch,/cvmfs ~/runs/eyu8/soft/HIBAG.sif Rscript script/HIBAG_script_Psychchip.R $allele $FILE $SUFFIX $REF $core

rm -r csv/$SUFFIX

mv $SUFFIX csv/

sed -i -- 's/"//g' csv/$SUFFIX/*
