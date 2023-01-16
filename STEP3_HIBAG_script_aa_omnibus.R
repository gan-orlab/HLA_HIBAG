library(data.table)

library(HIBAG)
library(readr)
library(lmtest)

#arg 1 is FILE arg2 is REGION
args <- commandArgs(trailingOnly = TRUE)

FILE <- args[1]
REGION <- args[2]


if(FILE != "ukbb"){

    covar <- as.data.frame(fread(paste0("/lustre03/project/6004655/COMMUN/runs/eyu8/data/HLA_typing/HIBAG/txt_data/", FILE, "/", FILE, "_covar.txt")))

    A <- read.csv(file=paste0("csv/", FILE, "_", REGION, "/HLA-A_", FILE, "_", REGION, ".csv"), sep=",",stringsAsFactors=FALSE)
    B <- read.csv(file=paste0("csv/", FILE, "_", REGION, "/HLA-B_", FILE, "_", REGION, ".csv"), sep=",",stringsAsFactors=FALSE)
    C <- read.csv(file=paste0("csv/", FILE, "_", REGION, "/HLA-C_", FILE, "_", REGION, ".csv"), sep=",",stringsAsFactors=FALSE)
    DPB1 <- read.csv(file=paste0("csv/", FILE, "_", REGION, "/HLA-DPB1_", FILE, "_", REGION, ".csv"), sep=",",stringsAsFactors=FALSE)
    DQA1 <- read.csv(file=paste0("csv/", FILE, "_", REGION, "/HLA-DQA1_", FILE, "_", REGION, ".csv"), sep=",",stringsAsFactors=FALSE)
    DQB1 <- read.csv(file=paste0("csv/", FILE, "_", REGION, "/HLA-DQB1_", FILE, "_", REGION, ".csv"), sep=",",stringsAsFactors=FALSE)
    DRB1 <- read.csv(file=paste0("csv/", FILE, "_", REGION, "/HLA-DRB1_", FILE, "_", REGION, ".csv"), sep=",",stringsAsFactors=FALSE)

} else if(FILE == "ukbb" && REGION == "PD"){

    PD <- as.data.frame(fread("/lustre03/project/6004655/COMMUN/runs/eyu8/data/HLA_typing/HIBAG/ukbb/ukbb_PD_covar.txt"))
    control_PD <- as.data.frame(fread("/lustre03/project/6004655/COMMUN/runs/eyu8/data/HLA_typing/HIBAG/ukbb/ukbb_control_PD_covar.txt"))

    PD$phenotype <- 2
    control_PD$phenotype <- 1

    covar <- rbind(PD, control_PD)

    A <- as.data.frame(fread("ukbb_imp_HLA_Euro/HLA-A_ukbb_imp_HLA_Euro.csv"))
    B <- as.data.frame(fread("ukbb_imp_HLA_Euro/HLA-B_ukbb_imp_HLA_Euro.csv"))
    C <- as.data.frame(fread("ukbb_imp_HLA_Euro/HLA-C_ukbb_imp_HLA_Euro.csv"))
    DPB1 <- as.data.frame(fread("ukbb_imp_HLA_Euro/HLA-DPB1_ukbb_imp_HLA_Euro.csv"))
    DQA1 <- as.data.frame(fread("ukbb_imp_HLA_Euro/HLA-DQA1_ukbb_imp_HLA_Euro.csv"))
    DQB1 <- as.data.frame(fread("ukbb_imp_HLA_Euro/HLA-DQB1_ukbb_imp_HLA_Euro.csv"))
    DRB1 <- as.data.frame(fread("ukbb_imp_HLA_Euro/HLA-DRB1_ukbb_imp_HLA_Euro.csv"))


} else if(FILE == "ukbb" && REGION == "Proxy"){

    Proxy <- as.data.frame(fread("/lustre03/project/6004655/COMMUN/runs/eyu8/data/HLA_typing/HIBAG/ukbb/ukbb_proxy_covar.txt"))
    control_Proxy <- as.data.frame(fread("/lustre03/project/6004655/COMMUN/runs/eyu8/data/HLA_typing/HIBAG/ukbb/ukbb_control_proxy_covar.txt"))

    Proxy$phenotype <- 2
    control_Proxy$phenotype <- 1

    covar <- rbind(Proxy, control_Proxy)

    A <- as.data.frame(fread("ukbb_imp_HLA_Euro/HLA-A_ukbb_imp_HLA_Euro.csv"))
    B <- as.data.frame(fread("ukbb_imp_HLA_Euro/HLA-B_ukbb_imp_HLA_Euro.csv"))
    C <- as.data.frame(fread("ukbb_imp_HLA_Euro/HLA-C_ukbb_imp_HLA_Euro.csv"))
    DPB1 <- as.data.frame(fread("ukbb_imp_HLA_Euro/HLA-DPB1_ukbb_imp_HLA_Euro.csv"))
    DQA1 <- as.data.frame(fread("ukbb_imp_HLA_Euro/HLA-DQA1_ukbb_imp_HLA_Euro.csv"))
    DQB1 <- as.data.frame(fread("ukbb_imp_HLA_Euro/HLA-DQB1_ukbb_imp_HLA_Euro.csv"))
    DRB1 <- as.data.frame(fread("ukbb_imp_HLA_Euro/HLA-DRB1_ukbb_imp_HLA_Euro.csv"))

}

HLA <- list(A=merge(covar,A),B=merge(covar,B),C=merge(covar,C),
            DPB1=merge(covar,DPB1),DQA1=merge(covar,DQA1),DQB1=merge(covar,DQB1),DRB1=merge(covar,DRB1))

poorSamples <- lapply(HLA,function(x) x[x$prob<0.5,]$"sample.id")
poorSamples <- unlist(poorSamples, use.names=FALSE)
poorSamples <- poorSamples[duplicated(poorSamples)]
filtered_HLA <- lapply(HLA,function(x) x[!(x$"sample.id" %in% poorSamples),])

for(i in 1:7){
    message(paste0("Analysing HLA-",names(HLA)[i]))
    hla_group <- filtered_HLA[[i]]
    hla_group$phenotype[hla_group$phenotype == -9] <- NA
    hla_group$phenotype <- as.factor(hla_group$phenotype - 1)
    hla <- hlaAllele(hla_group$"sample.id", H1=hla_group$allele1, H2=hla_group$allele2, locus=names(HLA[i]), assembly="hg19", prob=hla_group$prob)
    name <- names(HLA[i])
    aa <- hlaConvSequence(hla, code="P.code.merge")
    hla_aa <- aa$value
    hla_aa_table <- summary(aa)

    #Remove low prob
    filtered_HLA_gene <- hla_aa[hla_aa$prob > 0.5,]

    #Increment position
    increment <- hla_aa_table[1,"Pos"] - 1
    hla_aa_table[,"Pos"] <- hla_aa_table[,"Pos"] - increment

    hla_pos <- hla_aa_table[,"Pos"]
    HLA_allele1 <- filtered_HLA_gene[, c("sample.id","allele1")]
    HLA_allele2 <- filtered_HLA_gene[, c("sample.id","allele2")]
    names(HLA_allele1) <- c("ID","allele")
    names(HLA_allele2) <- c("ID","allele")
    HLA_allele <- rbind(HLA_allele1, HLA_allele2)

    hla_aa_result <- lapply(hla_pos, function(pos){
            HLA_aa_code <- HLA_allele
            HLA_aa_code$allele <- substr(HLA_aa_code$allele, pos, pos)

            HLA_count <- as.data.frame.matrix(table(HLA_aa_code))
            names(HLA_count)[names(HLA_count) == "-"] <- "Ref"
            names(HLA_count)[names(HLA_count) == "*"] <- "Amb"
            names(HLA_count)[names(HLA_count) == "."] <- "CNV"
            HLA_allele_name <- names(HLA_count)
            HLA_count$"sample.id" <- rownames(HLA_count)
            if(any(HLA_count$All == 0)){
                HLA_count[HLA_count$All == 0,]$All <- 1
            }
            if(any(HLA_count$All == 2)){
                HLA_count[HLA_count$All == 2,]$All <- 0
            }
            HLA_count[HLA_count == 2] <- 1
            HLA_meta_data <- hla_group[, !(names(hla_group) %in% c("allele1", "allele2", "prob", "matching"))]
            HLA_carrier <- merge(HLA_meta_data, HLA_count)


            HLA_count <- HLA_aa_code
            if(length(na.omit(HLA_count[HLA_count == "-"])) == 0){
                return(NULL)
            }
            HLA_count[HLA_count == "-"] <- "Ref"
            HLA_count[HLA_count == "*"] <- "Amb"
            HLA_count[HLA_count == "."] <- "CNV"

            names(HLA_count)[1] <- "sample.id"
            HLA_meta_data <- hla_group[, !(names(hla_group) %in% c("allele1", "allele2", "prob", "matching"))]
            HLA_glm <- merge(HLA_meta_data, HLA_count)
            HLA_glm$sex <- as.factor(HLA_glm$sex)
            HLA_glm$allele <- as.factor(HLA_glm$allele)
            HLA_glm$allele <- relevel(HLA_glm$allele, ref="Ref")

            if(length(levels(HLA_glm$allele)) == 1){
                        return(NULL)
            }
            f <- formula(phenotype ~ allele + age + sex + pc1 + pc2 + pc3 + pc4 + pc5 + pc6 + pc7 + pc8 + pc9 + pc10)
            if(FILE == "ukbb"){
                f <- formula(phenotype ~ allele + Townsend + AgeAtRecruit + sex + pc1 + pc2 + pc3 + pc4 + pc5 + pc6 + pc7 + pc8 + pc9 + pc10)
            }
            fnull <- formula(phenotype ~ age + sex + pc1 + pc2 + pc3 + pc4 + pc5 + pc6 + pc7 + pc8 + pc9 + pc10)
            if(FILE == "ukbb"){
                fnull <- formula(phenotype ~ Townsend + AgeAtRecruit + sex + pc1 + pc2 + pc3 + pc4 + pc5 + pc6 + pc7 + pc8 + pc9 + pc10)
            }
            HLA_result <- glm(f,family = binomial, data = HLA_glm)
            null <- glm(fnull,family = binomial, data=HLA_result$model)
            anova_result <- lrtest( null, HLA_result)
            Omnibus_P <- anova_result$Pr[2]

            HLA_table <- lapply(HLA_allele_name, function(allele){

                if(length(levels(as.factor(HLA_carrier[,allele]))) == 1){
                        return(NULL)
                }

                HLA_carrier_freq <-  table(HLA_carrier[, c(allele, "phenotype")])
                HLA_freq <- signif((HLA_carrier_freq[2,1] + HLA_carrier_freq[2,2])/sum(HLA_carrier_freq),3)
                HLA_carrier_freq_case <- signif(HLA_carrier_freq[2,2]/(HLA_carrier_freq[1,2] + HLA_carrier_freq[2,2]),3)
                HLA_carrier_freq_control <- signif(HLA_carrier_freq[2,1]/(HLA_carrier_freq[1,1] + HLA_carrier_freq[2,1]),3)
                HLA_ncase <- HLA_carrier_freq[1,2] + HLA_carrier_freq[2,2]
                HLA_ncontrol <- HLA_carrier_freq[1,1] + HLA_carrier_freq[2,1]
                HLA_ntotal <- sum(HLA_carrier_freq)
                index <- grep(paste0("allele",allele,"$"),rownames(coef(summary(HLA_result))))
                if(length(index) == 0){
                    return(NULL)
                }

                
                CI <- signif(as.data.frame(confint.default(HLA_result))[index,],3)
                var <- signif(t(as.data.frame(coef(summary(HLA_result))[index,])),3)
                var[,3] <- paste(CI[1,1], CI[1,2], sep = "-")
                colnames(var) <- c("b", "StdErr", "95%_CI", "p")
                var <-  as.data.frame(var)
                return(cbind(Pos = paste0("Pos",pos + increment,allele), Ref = substr(aa$reference,pos,pos), Alt = allele, HLA_ntotal, HLA_ncase, HLA_ncontrol, HLA_freq, HLA_carrier_freq_case, HLA_carrier_freq_control, var, Omnibus_P))})

            result <- Reduce(rbind, HLA_table)
            return(result)
            })

    hla_result <- as.data.frame(Reduce(rbind, hla_aa_result))
    write_delim(hla_result,paste0("HLA-",names(HLA)[i],"_",FILE,"_",REGION,"_aa_omnibus.csv"), delim = ",", na = "NA", append = FALSE, col_names = TRUE, escape = "double")
}

