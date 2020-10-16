#Remove duplicates and tracks snps from the bim file compoared to the reference panel

args <- commandArgs(trailingOnly = TRUE)


hlaSNP <- read.csv(args[1], sep = "\t", stringsAsFactors = FALSE, header= FALSE )

bimSNP <- read.csv(args[2], sep = "\t", stringsAsFactors = FALSE, header= FALSE )

#For each snp in refrence panel
for(i in 1:nrow(hlaSNP)){
    #Check if SNP is part of the bim file using chr position information
    if(hlaSNP[i,]$V1 %in% bimSNP$V4){
        #Get all alleles at position
        snp <- bimSNP[bimSNP$V4 == hlaSNP[i,]$V1,]
        #For each allele
        for(j in 1:nrow(snp)){
            pos <- snp[j,]
            allele1 <- pos$V5
            allele2 <- pos$V6
#In bim file, change indels into substitution
            if(substr(allele1,2,nchar(allele1)) == substr(allele2,2,nchar(allele2))){
                allele1 <- substr(allele1,1,1)
                allele2 <- substr(allele2,1,1)
                bimSNP[bimSNP$V4 == hlaSNP[i,]$V1,]$V5[j] <- allele1
                bimSNP[bimSNP$V4 == hlaSNP[i,]$V1,]$V6[j] <- allele2
            }

            if(allele1 == hlaSNP[i,]$V2 & allele2 == hlaSNP[i,]$V3){

            } else if (allele2 == hlaSNP[i,]$V2 & allele1 == hlaSNP[i,]$V3){
#Try flipping strand
            } else{
                if(pos$V5 == "A"){
                    pos$V5 = "T"
                } else if(pos$V5 == "T"){
                    pos$V5 = "A"
                } else if(pos$V5 == "G"){
                    pos$V5 = "C"
                } else if(pos$V5 == "C"){
                    pos$V5 = "G"
                }
                if(pos$V6 == "A"){
                    pos$V6 = "T"
                } else if(pos$V6 == "T"){
                    pos$V6 = "A"
                } else if(pos$V6 == "G"){
                    pos$V6 = "C"
                } else if(pos$V6 == "C"){
                    pos$V6 = "G"
                }
#If flipping doesn't work, label as "FRED" to exclude later                
                if(pos$V5 == hlaSNP[i,]$V2 & pos$V6 == hlaSNP[i,]$V3){
                    bimSNP[bimSNP$V4 == hlaSNP[i,]$V1,] <- pos
                } else if (pos$V6 == hlaSNP[i,]$V2 & pos$V5 == hlaSNP[i,]$V3){
                    bimSNP[bimSNP$V4 == hlaSNP[i,]$V1,] <- pos
                } else{
                    bimSNP[bimSNP$V4 == hlaSNP[i,]$V1,]$V2[j] <- "FRED"
                }               
            }
        }
    }
}

write.table(bimSNP, file='a', quote=FALSE, sep='\t', col.names = FALSE, row.names = FALSE)
