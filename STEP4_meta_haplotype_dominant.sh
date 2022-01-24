#!/bin/bash

allel=$1

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
echo "PROCESS      haplotype_${allel}_NIND_Euro_dominant.txt" >> meta.script
echo "PROCESS      haplotype_${allel}_IPDGC_Euro_dominant.txt" >> meta.script
echo "PROCESS      haplotype_${allel}_APDGC_Euro_dominant.txt" >> meta.script
echo "PROCESS      haplotype_${allel}_Oslo_Euro_dominant.txt" >> meta.script
echo "PROCESS      haplotype_${allel}_MCGILL_Euro_dominant.txt" >> meta.script
echo "PROCESS      haplotype_${allel}_PPMI_Euro_dominant.txt" >> meta.script
echo "PROCESS      haplotype_${allel}_NGRC_Euro_dominant.txt" >> meta.script
echo "PROCESS      ukbb/haplotype_${allel}_ukbb_PD_dominant.txt" >> meta.script
echo "EFFECT b_adjusted" >> meta.script
echo "STDERR StdErr_adjusted" >> meta.script
echo "PROCESS      ukbb/haplotype_${allel}_ukbb_Proxy_dominant.txt" >> meta.script
echo "OUTFILE   METAL-haplotype_${allel}_dominant .tbl" >> meta.script
echo "ANALYZE HETEROGENEITY" >> meta.script
echo "QUIT" >> meta.script


metal meta.script

cut -f1,4- METAL-haplotype_${allel}_dominant1.tbl | sed 's/MarkerName/Allele/'  > METAL-haplotype_${allel}_dominant.tbl

rm METAL-haplotype_${allel}_dominant1.tbl
mv METAL-haplotype_${allel}_dominant1.tbl.info METAL-haplotype_${allel}_dominant.tbl.info
