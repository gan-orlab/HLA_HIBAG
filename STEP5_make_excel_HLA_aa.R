#module load gcc/7.3 r-bundle-bioconductor/3.9

library(data.table)
library(readr)
library(openxlsx)

args <- commandArgs(trailingOnly = TRUE)

hla <- as.data.frame(fread(args[1]))
pos9 <- as.data.frame(fread(args[2]))


hla <- hla[,c(1:4,6:11,15)]
pos9 <- pos9[,c(1,8:11,15)]



result <- merge(hla, pos9, by = "MarkerName", all.x = T )
names(result) <- c("Amino Acid", "Ref", "Alt", "Carrier Freq", "MinFreq", "MaxFreq", "Effect", "StdErr", "P-value", "Direction", "HetPVal",
	"Effect", "StdErr", "P-value", "Direction", "HetPVal")
result[,3] <- "Other"

write.xlsx(result,paste0(args[3]))
