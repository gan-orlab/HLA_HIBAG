#!/bin/bash

allel=$1

module load metal

echo "MARKER Pos" > meta.script
echo "ALLELE Alt Ref" >> meta.script
echo "PVALUE aa_p" >> meta.script
echo "EFFECT aa_b" >> meta.script
echo "WEIGHT HLA_ntotal" >> meta.script
echo "SCHEME STDERR" >> meta.script
echo "STDERR aa_StdErr" >> meta.script
echo "FREQLABEL HLA_freq" >> meta.script
echo "AVERAGEFREQ ON" >> meta.script
echo "MINMAXFREQ ON" >> meta.script
echo "PROCESS      result/HLA-${allel}_IPDGC_Euro_aa.txt" >> meta.script
echo "PROCESS      result/HLA-${allel}_NIND_Euro_aa.txt" >> meta.script
echo "PROCESS      result/HLA-${allel}_APDGC_Euro_aa.txt" >> meta.script
echo "PROCESS      result/HLA-${allel}_Oslo_Euro_aa.txt" >> meta.script
echo "PROCESS      result/HLA-${allel}_MCGILL_Euro_aa.txt" >> meta.script
echo "PROCESS      result/HLA-${allel}_PPMI_Euro_aa.txt" >> meta.script
echo "PROCESS      result/HLA-${allel}_NGRC_Euro_aa.txt" >> meta.script
echo "PROCESS      ukbb/HLA-${allel}_ukbb_PD_aa.txt" >> meta.script
echo "OUTFILE	METAL-HLA-${allel} _aa.tbl" >> meta.script
echo "ANALYZE HETEROGENEITY" >> meta.script
echo "QUIT" >> meta.script


metal meta.script

mv METAL-HLA-${allel}1_aa.tbl METAL-HLA-${allel}_aa_noproxy.tbl

mv METAL-HLA-${allel}1_aa.tbl.info METAL-HLA-${allel}_aa_noproxy.tbl.info
