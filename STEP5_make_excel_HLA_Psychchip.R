#module load gcc/7.3 r-bundle-bioconductor/3.9

library(data.table)
library(readr)
library(openxlsx)

args <- commandArgs(trailingOnly = TRUE)

hla <- as.data.frame(fread(args[1]))
rs <- as.data.frame(fread(args[2]))
dqa1_0301 <- as.data.frame(fread(args[3]))
dqa1_0303 <- as.data.frame(fread(args[4]))
prs <- as.data.frame(fread(args[5]))

hla <- hla[,c(1,2,4:9,13)]
rs <- rs[,c(1,6:9,13)]
names(rs) <- c("Allele","V1","V2","V3","V4","V5")
dqa1_0301 <- dqa1_0301[,c(1,6:9,13)]
names(dqa1_0301) <- c("Allele","V10","V20","V30","V40","V50")
dqa1_0303 <- dqa1_0303[,c(1,6:9,13)]
names(dqa1_0303) <- c("Allele","V100","V200","V300","V400","V500")
prs <- prs[,c(1,6:9,13)]
names(prs) <- c("Allele","V1000","V2000","V3000","V4000","V5000")

result <- merge(merge(merge(merge(hla, rs, by = "Allele", all.x = T ), dqa1_0301, by = "Allele", all.x = T ),
	dqa1_0303, by = "Allele", all.x = T ), prs, by = "Allele", all.x = T )
names(result) <- c("Allele", "Carrier Freq", "MinFreq", "MaxFreq", "Effect", "StdErr", "P-value", "Direction", "HetPVal",
 "Effect", "StdErr", "P-value", "Direction", "HetPVal",
 "Effect", "StdErr", "P-value", "Direction", "HetPVal",
 "Effect", "StdErr", "P-value", "Direction", "HetPVal",
 "Effect", "StdErr", "P-value", "Direction", "HetPVal")


write.xlsx(result,paste0(args[6]))
