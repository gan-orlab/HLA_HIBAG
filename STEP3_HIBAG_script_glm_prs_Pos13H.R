#module load gcc/7.3 r-bundle-bioconductor/3.9

.libPaths(c("~/runs/eyu8/library/HIBAG",
"/cvmfs/soft.computecanada.ca/easybuild/software/2017/avx512/Compiler/gcc7.3/r-bundle-bioconductor/3.9",
"/cvmfs/soft.computecanada.ca/easybuild/software/2017/avx512/Compiler/gcc7.3/r/3.6.0/lib64/R/library"))

library(HIBAG)
library(readr)
library(data.table)

#arg 1 is FILE arg2 is REGION
args <- commandArgs(trailingOnly = TRUE)



FILE <- args[1]
REGION <- args[2]
if(FILE != "ukbb"){

    prs <- as.data.frame(fread(paste0("~/runs/eyu8/data/PRS/HLA/", FILE,"/", FILE, "_PRS_zscore.txt")))
    prs$sample.id <- prs$IID

    covar <- as.data.frame(fread(paste0("~/runs/eyu8/data/HLA_typing/HIBAG/txt_data/", FILE, "/", FILE, "_covar.txt")))

    A <- read.csv(file=paste0("csv/", FILE, "_", REGION, "/HLA-A_", FILE, "_", REGION, ".csv"), sep=",",stringsAsFactors=FALSE)
    B <- read.csv(file=paste0("csv/", FILE, "_", REGION, "/HLA-B_", FILE, "_", REGION, ".csv"), sep=",",stringsAsFactors=FALSE)
    C <- read.csv(file=paste0("csv/", FILE, "_", REGION, "/HLA-C_", FILE, "_", REGION, ".csv"), sep=",",stringsAsFactors=FALSE)
    DPB1 <- read.csv(file=paste0("csv/", FILE, "_", REGION, "/HLA-DPB1_", FILE, "_", REGION, ".csv"), sep=",",stringsAsFactors=FALSE)
    DQA1 <- read.csv(file=paste0("csv/", FILE, "_", REGION, "/HLA-DQA1_", FILE, "_", REGION, ".csv"), sep=",",stringsAsFactors=FALSE)
    DQB1 <- read.csv(file=paste0("csv/", FILE, "_", REGION, "/HLA-DQB1_", FILE, "_", REGION, ".csv"), sep=",",stringsAsFactors=FALSE)
    DRB1 <- read.csv(file=paste0("csv/", FILE, "_", REGION, "/HLA-DRB1_", FILE, "_", REGION, ".csv"), sep=",",stringsAsFactors=FALSE)
    DRB3 <- read.csv(file=paste0("csv/", FILE, "_DRB3/HLA-DRB3_", FILE, "_DRB3.csv"), sep=",",stringsAsFactors=FALSE)
    DRB4 <- read.csv(file=paste0("csv/", FILE, "_DRB4/HLA-DRB4_", FILE, "_DRB4.csv"), sep=",",stringsAsFactors=FALSE)
    DRB5 <- read.csv(file=paste0("csv/", FILE, "_DRB5/HLA-DRB5_", FILE, "_DRB5.csv"), sep=",",stringsAsFactors=FALSE)

} else if(FILE == "ukbb" && REGION == "PD"){

    prs <- as.data.frame(fread(paste0("~/runs/eyu8/data/PRS/HLA/ukbb/ukbb_PRS_zscore.txt")))
    names(prs)[1] <- "sample.id"

    PD <- as.data.frame(fread("~/runs/eyu8/data/ukbb_dagher/ukbb_PD_covar.txt"))
    control_PD <- as.data.frame(fread("~/runs/eyu8/data/HLA_typing/HIBAG/ukbb/ukbb_control_PD_covar.txt"))

    PD$phenotype <- 2
    control_PD$phenotype <- 1

    names(PD)[1] <- "sample.id"
    names(control_PD)[1] <- "sample.id"

    covar <- rbind(PD, control_PD)

    A <- as.data.frame(fread("ukbb_imp_HLA_Euro/HLA-A_ukbb_imp_HLA_Euro.csv"))
    B <- as.data.frame(fread("ukbb_imp_HLA_Euro/HLA-B_ukbb_imp_HLA_Euro.csv"))
    C <- as.data.frame(fread("ukbb_imp_HLA_Euro/HLA-C_ukbb_imp_HLA_Euro.csv"))
    DPB1 <- as.data.frame(fread("ukbb_imp_HLA_Euro/HLA-DPB1_ukbb_imp_HLA_Euro.csv"))
    DQA1 <- as.data.frame(fread("ukbb_imp_HLA_Euro/HLA-DQA1_ukbb_imp_HLA_Euro.csv"))
    DQB1 <- as.data.frame(fread("ukbb_imp_HLA_Euro/HLA-DQB1_ukbb_imp_HLA_Euro.csv"))
    DRB1 <- as.data.frame(fread("ukbb_imp_HLA_Euro/HLA-DRB1_ukbb_imp_HLA_Euro.csv"))
    DRB3 <- as.data.frame(fread("ukbb_imp_HLA_DRB3/HLA-DRB3_ukbb_imp_HLA_DRB3.csv"))
    DRB4 <- as.data.frame(fread("ukbb_imp_HLA_DRB4/HLA-DRB4_ukbb_imp_HLA_DRB4.csv"))
    DRB5 <- as.data.frame(fread("ukbb_imp_HLA_DRB5/HLA-DRB5_ukbb_imp_HLA_DRB5.csv"))


} else if(FILE == "ukbb" && REGION == "Proxy"){

    prs <- as.data.frame(fread(paste0("~/runs/eyu8/data/PRS/HLA/ukbb/ukbb_PRS_zscore.txt")))
    names(prs)[1] <- "sample.id"

    Proxy <- as.data.frame(fread("~/runs/eyu8/data/ukbb_dagher/ukbb_proxy_covar.txt"))
    control_Proxy <- as.data.frame(fread("~/runs/eyu8/data/HLA_typing/ukbb/ukbb_control_proxy_covar.txt"))

    Proxy$phenotype <- 2
    control_Proxy$phenotype <- 1

    names(Proxy)[1] <- "sample.id"
    names(control_Proxy)[1] <- "sample.id"

    covar <- rbind(Proxy, control_Proxy)

    A <- as.data.frame(fread("ukbb_imp_HLA_Euro/HLA-A_ukbb_imp_HLA_Euro.csv"))
    B <- as.data.frame(fread("ukbb_imp_HLA_Euro/HLA-B_ukbb_imp_HLA_Euro.csv"))
    C <- as.data.frame(fread("ukbb_imp_HLA_Euro/HLA-C_ukbb_imp_HLA_Euro.csv"))
    DPB1 <- as.data.frame(fread("ukbb_imp_HLA_Euro/HLA-DPB1_ukbb_imp_HLA_Euro.csv"))
    DQA1 <- as.data.frame(fread("ukbb_imp_HLA_Euro/HLA-DQA1_ukbb_imp_HLA_Euro.csv"))
    DQB1 <- as.data.frame(fread("ukbb_imp_HLA_Euro/HLA-DQB1_ukbb_imp_HLA_Euro.csv"))
    DRB1 <- as.data.frame(fread("ukbb_imp_HLA_Euro/HLA-DRB1_ukbb_imp_HLA_Euro.csv"))
    DRB3 <- as.data.frame(fread("ukbb_imp_HLA_DRB3/HLA-DRB3_ukbb_imp_HLA_DRB3.csv"))
    DRB4 <- as.data.frame(fread("ukbb_imp_HLA_DRB4/HLA-DRB4_ukbb_imp_HLA_DRB4.csv"))
    DRB5 <- as.data.frame(fread("ukbb_imp_HLA_DRB5/HLA-DRB5_ukbb_imp_HLA_DRB5.csv"))
}

HLA <- list(A=merge(prs,merge(covar,A)),B=merge(prs,merge(covar,B)),C=merge(prs,merge(covar,C)),
    DPB1=merge(prs,merge(covar,DPB1)),DQA1=merge(prs,merge(covar,DQA1)),DQB1=merge(prs,merge(covar,DQB1)),DRB1=merge(prs,merge(covar,DRB1)),
    DRB3=merge(prs,merge(covar,DRB3)),DRB4=merge(prs,merge(covar,DRB4)),DRB5=merge(prs,merge(covar,DRB5)))

poorSamples <- lapply(HLA,function(x) x[x$prob<0.5,]$sample.id)
poorSamples <- unlist(poorSamples, use.names=FALSE)
poorSamples <- poorSamples[duplicated(poorSamples)]
filtered_HLA <- lapply(HLA,function(x) x[!(x$sample.id %in% poorSamples),])

hla_group <- filtered_HLA[[7]]
hla_group$phenotype[hla_group$phenotype == -9] <- NA
hla_group$phenotype <- as.factor(hla_group$phenotype - 1)
hla <- hlaAllele(hla_group$sample.id, H1=hla_group$allele1, H2=hla_group$allele2, locus=names(HLA[7]), assembly="hg19", prob=hla_group$prob)
name <- names(HLA[7])
aa <- hlaConvSequence(hla, code="P.code.merge")
hla_aa <- aa$value
hla_aa_table <- summary(aa)
filtered_HLA_gene <- hla_aa[hla_aa$prob > 0.5,]
increment <- hla_aa_table[1,"Pos"] - 1
hla_aa_table[,"Pos"] <- hla_aa_table[,"Pos"] - increment
hla_pos <- hla_aa_table[,"Pos"]
HLA_allele1 <- filtered_HLA_gene[, c("sample.id","allele1")]
HLA_allele2 <- filtered_HLA_gene[, c("sample.id","allele2")]
names(HLA_allele1) <- c("sample.id","allele")
names(HLA_allele2) <- c("sample.id","allele")
HLA_allele <- rbind(HLA_allele1, HLA_allele2)
HLA_aa_code <- HLA_allele

hla_pos <- function(pos){
    HLA_aa_code$allele <- substr(HLA_aa_code$allele, 9, 9)
    HLA_count <- as.data.frame.matrix(table(HLA_aa_code))
    names(HLA_count)[names(HLA_count) == "-"] <- "Ref"
    names(HLA_count)[names(HLA_count) == "*"] <- "Amb"
    names(HLA_count)[names(HLA_count) == "."] <- "CNV"
    HLA_allele_name <- names(HLA_count)
    HLA_count$sample.id <- rownames(HLA_count)
    if(any(HLA_count$All == 0)){
        HLA_count[HLA_count$All == 0,]$All <- 1
    }
    if(any(HLA_count$All == 2)){
        HLA_count[HLA_count$All == 2,]$All <- 0
    }
    HLA_count[HLA_count == 2] <- 1
    HLA_meta_data <- hla_group[, !(names(hla_group) %in% c("allele1", "allele2", "prob", "matching"))]
    HLA_glm <- merge(HLA_meta_data, HLA_count)
    return(HLA_glm)
}

hla_drb1_pos13H <- hla_pos(13)[, c("sample.id", "H")]
names(hla_drb1_pos13H)[2] <- "DRB1_Pos13H"

for(i in 1:10){

   filtered_HLA_gene <- filtered_HLA[[i]][filtered_HLA[[i]]$prob > 0.5,]
    HLA_allele1 <- filtered_HLA_gene[, c("sample.id","allele1")]
    HLA_allele2 <- filtered_HLA_gene[, c("sample.id","allele2")]
    names(HLA_allele1) <- c("sample.id","allele")
    names(HLA_allele2) <- c("sample.id","allele")
    HLA_allele <- rbind(HLA_allele1, HLA_allele2)
    HLA_count <- as.data.frame.matrix(table(HLA_allele))
    HLA_allele_name <- names(HLA_count)
    HLA_count$sample.id <- rownames(HLA_count)
    HLA_count[HLA_count == 2] <- 1
    HLA_meta_data <- filtered_HLA_gene[, !(names(filtered_HLA_gene) %in% c("allele1", "allele2", "prob", "matching"))]
    HLA_glm <- merge(merge(HLA_meta_data, HLA_count), hla_drb1_pos13H)
    HLA_table <- lapply(HLA_allele_name, function(allele){
        HLA_glm$phenotype[HLA_glm$phenotype == -9] <- NA
        HLA_glm$phenotype <- as.factor(HLA_glm$phenotype - 1)
        HLA_glm$sex <- as.factor(HLA_glm$sex)
        names(HLA_glm)[names(HLA_glm) == allele] <- "h"
        if(all(HLA_glm[,"h"] == HLA_glm[,"DRB1_Pos13H"])){
                return(NULL)
        }
        HLA_glm$h <- as.numeric(HLA_glm$h)
        if(length(levels(factor(HLA_glm$h))) == 1){
        	return(NULL)
        }

        if(FILE != "ukbb"){

            f <- formula(paste("phenotype ~ h + Z * DRB1_Pos13H + age + sex + pc1 + pc2 + pc3 + pc4 + pc5 + pc6 + pc7 + pc8 + pc9 + pc10"))

        } else if(FILE == "ukbb"){

            f <- formula(paste("phenotype ~ h + Z * DRB1_Pos13H + Townsend + AgeAtRecruit + sex + pc1 + pc2 + pc3 + pc4 + pc5 + pc6 + pc7 + pc8 + pc9 + pc10"))

        }

        HLA_result <- glm(f,family = binomial, data = HLA_glm)
        HLA_carrier_freq <-  table(HLA_glm[, c("h", "phenotype")])
        HLA_freq <- signif((HLA_carrier_freq[2,1] + HLA_carrier_freq[2,2])/sum(HLA_carrier_freq), 3)
        HLA_carrier_freq_case <- signif(HLA_carrier_freq[2,2]/(HLA_carrier_freq[1,2] + HLA_carrier_freq[2,2]), 3)
        HLA_carrier_freq_control <- signif(HLA_carrier_freq[2,1]/(HLA_carrier_freq[1,1] + HLA_carrier_freq[2,1]), 3)
        HLA_ncase <- HLA_carrier_freq[1,2] + HLA_carrier_freq[2,2]
        HLA_ncontrol <- HLA_carrier_freq[1,1] + HLA_carrier_freq[2,1]
        HLA_ntotal <- sum(HLA_carrier_freq)
        CI <- signif(as.data.frame(confint.default(HLA_result)), 3)
        var <- signif(as.data.frame(coef(summary(HLA_result))), 3)
        var_list <- lapply(2:nrow(var), function(i){
            summary_stats <- cbind(var[i, 1, drop = FALSE], var[i, 2, drop = FALSE], var[i, 4, drop = FALSE])
            names(summary_stats) <- c("b", "StdErr", "p")
            names(summary_stats) <- paste(row.names(summary_stats), names(summary_stats), sep = "_")
            return(summary_stats)
            })
        if(FILE == "ukbb" && REGION == "Proxy"){

            var_list <- append(var_list[1], var_list)
            var_list[[2]] <- var_list[[2]] * c(2,2,1)
            names(var_list[[2]]) <- c("h_b_adjusted", "h_StdErr_adjusted", "h_p_adjusted")
        }

        return(cbind(allele, HLA_ntotal, HLA_ncase, HLA_ncontrol, HLA_freq, HLA_carrier_freq_case, HLA_carrier_freq_control, Reduce(cbind, var_list)))})
    result <- Reduce(rbind, HLA_table)
    write_delim(result,paste0("HLA-",names(HLA)[i],"_",FILE,"_",REGION,"_prs_Pos13H.txt"), delim = " ", na = "NA", append = FALSE, col_names = TRUE, quote_escape = "double")
}
