#Extract snp and allele from gene
d <- getwd()


.libPaths(c("~/runs/eyu8/library/HIBAG",
            "/cvmfs/soft.computecanada.ca/easybuild/software/2017/avx512/Compiler/gcc7.3/r-bundle-bioconductor/3.9",
            "/cvmfs/soft.computecanada.ca/easybuild/software/2017/avx512/Compiler/gcc7.3/r/3.6.0/lib64/R/library"))

library("HIBAG")
setwd(d)
args <- commandArgs(trailingOnly = TRUE)


model.list <- get(load(args[1]))

gene <- names(model.list)
for (n in gene){

    write.csv(model.list[[n]]$snp.position,paste("tmp/HLA-",n,"_snp_position.csv",sep=""),col.names=FALSE,row.names=FALSE)
    write.csv(model.list[[n]]$snp.allele,paste("tmp/HLA-",n,"_snp_allele.csv",sep=""),col.names=FALSE,row.names=FALSE)
}
