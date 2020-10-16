#!/bin/bash

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
echo "PROCESS      result/haplotype_DRB4_DRB1_DQA1_DQB1_NIND_Euro.txt" >> meta.script
echo "PROCESS      result/haplotype_DRB4_DRB1_DQA1_DQB1_IPDGC_Euro.txt" >> meta.script
echo "PROCESS      result/haplotype_DRB4_DRB1_DQA1_DQB1_APDGC_Euro.txt" >> meta.script
echo "PROCESS      result/haplotype_DRB4_DRB1_DQA1_DQB1_Oslo_Euro.txt" >> meta.script
echo "PROCESS      result/haplotype_DRB4_DRB1_DQA1_DQB1_MCGILL_Euro.txt" >> meta.script
echo "PROCESS      result/haplotype_DRB4_DRB1_DQA1_DQB1_PPMI_Euro.txt" >> meta.script
echo "PROCESS      result/haplotype_DRB4_DRB1_DQA1_DQB1_NGRC_Euro.txt" >> meta.script
echo "PROCESS      ukbb/haplotype_DRB4_DRB1_DQA1_DQB1_ukbb_PD.txt" >> meta.script
echo "EFFECT b_adjusted" >> meta.script
echo "STDERR StdErr_adjusted" >> meta.script
echo "PROCESS      ukbb/haplotype_DRB4_DRB1_DQA1_DQB1_ukbb_Proxy.txt" >> meta.script
echo "OUTFILE	METAL-haplotype_DRB4_DRB1_DQA1_DQB1 .tbl" >> meta.script
echo "ANALYZE HETEROGENEITY" >> meta.script
echo "QUIT" >> meta.script


metal meta.script

cut -f1,4- METAL-haplotype_DRB4_DRB1_DQA1_DQB11.tbl | sed 's/MarkerName/Allele/'  > METAL-haplotype_DRB4_DRB1_DQA1_DQB1.tbl

rm METAL-haplotype_DRB4_DRB1_DQA1_DQB11.tbl
mv METAL-haplotype_DRB4_DRB1_DQA1_DQB11.tbl.info METAL-haplotype_DRB4_DRB1_DQA1_DQB1.tbl.info
