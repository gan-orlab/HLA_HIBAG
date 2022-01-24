#!/bin/bash

module load metal

echo "MARKER Haplotype" > meta.script
echo "PVALUE haplo_p" >> meta.script
echo "EFFECT haplo_b" >> meta.script
echo "WEIGHT Haplo_ntotal" >> meta.script
echo "SCHEME STDERR" >> meta.script
echo "STDERR haplo_StdErr" >> meta.script
echo "FREQLABEL Haplo_freq" >> meta.script
echo "AVERAGEFREQ ON" >> meta.script
echo "MINMAXFREQ ON" >> meta.script
echo "PROCESS      haplotype_H13_H33_D57_NIND_Euro.txt" >> meta.script
echo "PROCESS      haplotype_H13_H33_D57_IPDGC_Euro.txt" >> meta.script
echo "PROCESS      haplotype_H13_H33_D57_APDGC_Euro.txt" >> meta.script
echo "PROCESS      haplotype_H13_H33_D57_Oslo_Euro.txt" >> meta.script
echo "PROCESS      haplotype_H13_H33_D57_MCGILL_Euro.txt" >> meta.script
echo "PROCESS      haplotype_H13_H33_D57_PPMI_Euro.txt" >> meta.script
echo "PROCESS      haplotype_H13_H33_D57_NGRC_Euro.txt" >> meta.script
echo "PROCESS      ukbb/haplotype_H13_H33_D57_ukbb_PD.txt" >> meta.script
echo "EFFECT b_adjusted" >> meta.script
echo "STDERR StdErr_adjusted" >> meta.script
echo "PROCESS      ukbb/haplotype_H13_H33_D57_ukbb_Proxy.txt" >> meta.script
echo "OUTFILE   METAL-haplotype_H13_H33_D57_dominant .tbl" >> meta.script
echo "ANALYZE HETEROGENEITY" >> meta.script
echo "QUIT" >> meta.script


metal meta.script

cut -f1,4- METAL-haplotype_H13_H33_D57_dominant1.tbl | sed 's/MarkerName/Allele/'  > METAL-haplotype_H13_H33_D57_dominant.tbl

rm METAL-haplotype_H13_H33_D57_dominant1.tbl
mv METAL-haplotype_H13_H33_D57_dominant1.tbl.info METAL-haplotype_H13_H33_D57_dominant.tbl.info
