#Extract snp and allele from gene

library(HIBAG)
args <- commandArgs(trailingOnly = TRUE)

model.list <- get(load(args[1]))

gene <- names(model.list)
for (n in gene){

    write.csv(model.list[[n]]$snp.position,paste("tmp/HLA-",n,"_snp_position.csv",sep=""),col.names=FALSE,row.names=FALSE)
    write.csv(model.list[[n]]$snp.allele,paste("tmp/HLA-",n,"_snp_allele.csv",sep=""),col.names=FALSE,row.names=FALSE)
}
