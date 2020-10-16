#!/bin/bash

allel=$1

module load metal

echo "MARKER Haplotype" > meta.script
echo "PVALUE p" >> meta.script
echo "EFFECT b" >> meta.script
echo "WEIGHT N" >> meta.script
echo "SCHEME STDERR" >> meta.script
echo "STDERR StdErr" >> meta.script
echo "FREQLABEL Haplo_freq" >> meta.script
echo "AVERAGEFREQ ON" >> meta.script
echo "MINMAXFREQ ON" >> meta.script
echo "PROCESS      result/haplotype_${allel}_NIND_Euro.txt" >> meta.script
echo "PROCESS      result/haplotype_${allel}_IPDGC_Euro.txt" >> meta.script
echo "PROCESS      result/haplotype_${allel}_APDGC_Euro.txt" >> meta.script
echo "PROCESS      result/haplotype_${allel}_Oslo_Euro.txt" >> meta.script
echo "PROCESS      result/haplotype_${allel}_MCGILL_Euro.txt" >> meta.script
echo "PROCESS      result/haplotype_${allel}_PPMI_Euro.txt" >> meta.script
echo "PROCESS      result/haplotype_${allel}_NGRC_Euro.txt" >> meta.script
echo "PROCESS      ukbb/haplotype_${allel}_ukbb_PD.txt" >> meta.script
echo "EFFECT b_adjusted" >> meta.script
echo "STDERR StdErr_adjusted" >> meta.script
echo "PROCESS      ukbb/haplotype_${allel}_ukbb_Proxy.txt" >> meta.script
echo "OUTFILE	METAL-haplotype_${allel} .tbl" >> meta.script
echo "ANALYZE HETEROGENEITY" >> meta.script
echo "QUIT" >> meta.script


metal meta.script

cut -f1,4- METAL-haplotype_${allel}1.tbl | sed 's/MarkerName/Allele/'  > METAL-haplotype_${allel}.tbl

rm METAL-haplotype_${allel}1.tbl
mv METAL-haplotype_${allel}1.tbl.info METAL-haplotype_${allel}.tbl.info
