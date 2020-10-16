library(data.table)

pc <- as.data.frame(fread(paste0("~/runs/go_lab/PPMI/cleaned/PPMI_pre-impute.eigenvec")))[,1:11]

names(pc) <- c("ID","pc1","pc2","pc3","pc4","pc5","pc6","pc7","pc8","pc9","pc10")

age <- as.data.frame(fread("~/runs/eyu8/data/HLA_typing/HIBAG/txt_data/PPMI/PPMI.age"))

names(age) <- c("ID","age")

pheno <- as.data.frame(fread("~/runs/eyu8/data/HLA_typing/HIBAG/bfile/final_PPMI_Euro.fam"))[,c(2,5,6)]

names(pheno) <- c("ID","sex","phenotype")

covar <- merge(pheno,merge(age,pc))
names(covar)[1] <- "sample.id"

write.csv(covar,"~/runs/eyu8/data/HLA_typing/HIBAG/txt_data/PPMI/PPMI_covar.txt", row.names = F, quote = F)


make_HLA <- function(allele){
    HLA <- as.data.frame(fread(paste0("~/runs/eyu8/data/HLA_typing/HIBAG/csv/PPMI_Euro/HLA-",allele,"_PPMI_Euro.csv")))
    names(HLA)[1] <- c("ID")
    HLA_covar <- merge(covar,HLA)
    write.csv(HLA_covar,paste0("~/runs/eyu8/data/HLA_typing/HIBAG/txt_data/PPMI/HLA-",allele,"_PPMI_Euro.logit"), row.names = F, quote = F)
}

#for(allele in c("A","B","C","DPB1","DQA1","DQB1","DRB1")){
#    make_HLA(allele)
#}
