#!/bin/bash

module load singularity

FILE=bfile/final_$1_$2

SUFFIX=$1_$2

REF=$3
core=$4


mkdir -p $SUFFIX
mkdir -p log/$SUFFIX

singularity exec -B /lustre03,/project,/scratch,/cvmfs ~/runs/eyu8/soft/HIBAG.sif Rscript script/HIBAG_script.R A $FILE $SUFFIX $REF $core  
singularity exec -B /lustre03,/project,/scratch,/cvmfs ~/runs/eyu8/soft/HIBAG.sif Rscript script/HIBAG_script.R B $FILE $SUFFIX $REF $core  
singularity exec -B /lustre03,/project,/scratch,/cvmfs ~/runs/eyu8/soft/HIBAG.sif Rscript script/HIBAG_script.R C $FILE $SUFFIX $REF $core  
singularity exec -B /lustre03,/project,/scratch,/cvmfs ~/runs/eyu8/soft/HIBAG.sif Rscript script/HIBAG_script.R DPB1 $FILE $SUFFIX $REF $core  
singularity exec -B /lustre03,/project,/scratch,/cvmfs ~/runs/eyu8/soft/HIBAG.sif Rscript script/HIBAG_script.R DQA1 $FILE $SUFFIX $REF $core  
singularity exec -B /lustre03,/project,/scratch,/cvmfs ~/runs/eyu8/soft/HIBAG.sif Rscript script/HIBAG_script.R DQB1 $FILE $SUFFIX $REF $core  
singularity exec -B /lustre03,/project,/scratch,/cvmfs ~/runs/eyu8/soft/HIBAG.sif Rscript script/HIBAG_script.R DRB1 $FILE $SUFFIX $REF $core  

rm -r csv/$SUFFIX

mv $SUFFIX csv/

sed -i -- 's/"//g' csv/$SUFFIX/*
