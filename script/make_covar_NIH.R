library(data.table)

args <- commandArgs(trailingOnly = TRUE)

FILE <- args[1]

pc <- as.data.frame(fread(paste0("~/runs/go_lab/dbGap/",FILE,"/cleaned/",FILE,"_pre-impute.eigenvec")))[,1:11]

names(pc) <- c("ID","pc1","pc2","pc3","pc4","pc5","pc6","pc7","pc8","pc9","pc10")

age <- as.data.frame(fread(paste0("~/runs/eyu8/data/HLA_typing/HIBAG/txt_data/",FILE,"/",FILE,".age")))

names(age) <- c("ID","age")

pheno <- as.data.frame(fread(paste0("~/runs/eyu8/data/HLA_typing/HIBAG/bfile/final_",FILE,"_Euro.fam")))[,c(2,5,6)]

names(pheno) <- c("ID","sex","phenotype")

covar <- merge(pheno,merge(age,pc))
names(covar)[1] <- "sample.id"
write.csv(covar,paste0("~/runs/eyu8/data/HLA_typing/HIBAG/txt_data/",FILE,"/",FILE,"_covar.txt"), row.names = F, quote = F)


make_HLA <- function(allele){
    HLA <- as.data.frame(fread(paste0("~/runs/eyu8/data/HLA_typing/HIBAG/csv/",FILE,"_Euro/HLA-",allele,"_",FILE,"_Euro.csv")))
    names(HLA)[1] <- c("ID")
    HLA_covar <- merge(covar,HLA)
    write.csv(HLA_covar,paste0("~/runs/eyu8/data/HLA_typing/HIBAG/txt_data/",FILE,"/HLA-",allele,"_",FILE,"_Euro.logit"), row.names = F, quote = F)
}

#for(allele in c("A","B","C","DPB1","DQA1","DQB1","DRB1")){
#    make_HLA(allele)
#}
