library(data.table)

pc <- as.data.frame(fread("~/runs/go_lab/omnix/MCGILL-RBD/cleaned/MCGILL-RBD_pre-impute.eigenvec"))[,1:11]

names(pc) <- c("ID","pc1","pc2","pc3","pc4","pc5","pc6","pc7","pc8","pc9","pc10")

age <- as.data.frame(fread("~/runs/eyu8/data/HLA_typing/HIBAG/txt_data/MCGILL-RBD/MCGILL-RBD.age"))

names(age) <- c("ID","age")

pheno <- as.data.frame(fread("~/runs/eyu8/data/HLA_typing/HIBAG/bfile/final_MCGILL-RBD_Euro.fam"))[,c(2,5,6)]

names(pheno) <- c("ID","sex","phenotype")

covar <- merge(pheno,merge(age,pc))
names(covar)[1] <- "sample.id"

write.csv(covar,"~/runs/eyu8/data/HLA_typing/HIBAG/txt_data/MCGILL-RBD/MCGILL-RBD_covar.txt", row.names = F, quote = F)

