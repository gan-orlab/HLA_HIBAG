#module load gcc/7.3 r-bundle-bioconductor/3.9

.libPaths(c("~/runs/eyu8/library/haplo.stats",
"/cvmfs/soft.computecanada.ca/easybuild/software/2017/avx512/Compiler/gcc7.3/r-bundle-bioconductor/3.9",
"/cvmfs/soft.computecanada.ca/easybuild/software/2017/avx512/Compiler/gcc7.3/r/3.6.0/lib64/R/library"))

library(haplo.stats)
library(readr)
library(data.table)

seed <- c(17, 53, 1, 40, 37, 0, 62, 56, 5, 52, 12, 1)
set.seed(seed)

args <- commandArgs(trailingOnly = TRUE)

#args <- c("ukbb","Proxy")

FILE <- args[1]
REGION <- args[2]
if(FILE != "ukbb"){

    covar <- as.data.frame(fread(paste0("~/runs/eyu8/data/HLA_typing/HIBAG/txt_data/", FILE, "/", FILE, "_covar.txt")))

    DQA1 <- read.csv(file=paste0("csv/", FILE, "_", REGION, "/HLA-DQA1_", FILE, "_", REGION, ".csv"), sep=",",stringsAsFactors=FALSE)
    DQB1 <- read.csv(file=paste0("csv/", FILE, "_", REGION, "/HLA-DQB1_", FILE, "_", REGION, ".csv"), sep=",",stringsAsFactors=FALSE)
    DRB1 <- read.csv(file=paste0("csv/", FILE, "_", REGION, "/HLA-DRB1_", FILE, "_", REGION, ".csv"), sep=",",stringsAsFactors=FALSE)
    DRB4 <- read.csv(file=paste0("csv/", FILE, "_DRB4/HLA-DRB4_", FILE, "_DRB4.csv"), sep=",",stringsAsFactors=FALSE)

} else if(FILE == "ukbb" && REGION == "PD"){

    PD <- as.data.frame(fread("~/runs/eyu8/data/ukbb_dagher/ukbb_PD_covar.txt"))
    control_PD <- as.data.frame(fread("~/runs/eyu8/data/HLA_typing/HIBAG/ukbb/ukbb_control_PD_covar.txt"))

    PD$phenotype <- 2
    control_PD$phenotype <- 1

    covar <- rbind(PD, control_PD)

    covar <- covar[!is.na(covar$Townsend),]

    DQA1 <- as.data.frame(fread("ukbb_imp_HLA_Euro/HLA-DQA1_ukbb_imp_HLA_Euro.csv"))
    DQB1 <- as.data.frame(fread("ukbb_imp_HLA_Euro/HLA-DQB1_ukbb_imp_HLA_Euro.csv"))
    DRB1 <- as.data.frame(fread("ukbb_imp_HLA_Euro/HLA-DRB1_ukbb_imp_HLA_Euro.csv"))
    DRB4 <- as.data.frame(fread("ukbb_imp_HLA_DRB4/HLA-DRB4_ukbb_imp_HLA_DRB4.csv"))


} else if(FILE == "ukbb" && REGION == "Proxy"){

    Proxy <- as.data.frame(fread("~/runs/eyu8/data/ukbb_dagher/ukbb_proxy_covar.txt"))
    control_Proxy <- as.data.frame(fread("~/runs/eyu8/data/HLA_typing/HIBAG/ukbb/ukbb_control_proxy_covar.txt"))

    Proxy$phenotype <- 2
    control_Proxy$phenotype <- 1

    covar <- rbind(Proxy, control_Proxy)

    covar <- covar[!is.na(covar$Townsend),]
  
    DQA1 <- as.data.frame(fread("ukbb_imp_HLA_Euro/HLA-DQA1_ukbb_imp_HLA_Euro.csv"))
    DQB1 <- as.data.frame(fread("ukbb_imp_HLA_Euro/HLA-DQB1_ukbb_imp_HLA_Euro.csv"))
    DRB1 <- as.data.frame(fread("ukbb_imp_HLA_Euro/HLA-DRB1_ukbb_imp_HLA_Euro.csv"))
    DRB4 <- as.data.frame(fread("ukbb_imp_HLA_DRB4/HLA-DRB4_ukbb_imp_HLA_DRB4.csv"))
}

HLA <- list(DRB4=merge(covar,DRB4),DRB1=merge(covar,DRB1),
    DQA1=merge(covar,DQA1),DQB1=merge(covar,DQB1))

poorSamples <- lapply(HLA,function(x) x[x$prob<0.5,]$sample.id)
poorSamples <- unlist(poorSamples, use.names=FALSE)
poorSamples <- poorSamples[duplicated(poorSamples)]

filtered_HLA <- lapply(HLA,function(x) x[!(x$sample.id %in% poorSamples), c("sample.id","allele1","allele2")])
##################

HLA.list <- lapply(names(filtered_HLA),function(gene){
        HLA_gene <- filtered_HLA[[gene]]
        names(HLA_gene) <- c("sample.id",paste(gene,c("a1","a2"),sep = "."))
        return(HLA_gene)
        })

HLA_merged <- merge(covar, Reduce(merge, HLA.list))

HLA_covar <- HLA_merged[, !(names(HLA_merged) %in% grep(".a",names(HLA_merged)))]


haplotype_glm <- function(pattern){

    geno <- HLA_merged[,grep(pattern,names(HLA_merged))]
    label <- names(HLA)[grep(pattern,names(HLA))]

    geno.glm <- setupGeno(geno, miss.val=c(0,NA), locus.label=label)

    glm.data <- cbind(HLA_covar,data.frame(geno.glm))

    glm.data$phenotype[glm.data$phenotype == -9] <- NA
    glm.data$phenotype <- glm.data$phenotype - 1
    glm.data$sex <- as.factor(glm.data$sex)

    if(FILE != "ukbb"){

        f <- formula(paste("phenotype ~ geno.glm + age + sex + pc1 + pc2 + pc3 + pc4 + pc5 + pc6 + pc7 + pc8 + pc9 + pc10"))

    } else if(FILE == "ukbb"){

        f <- formula(paste("phenotype ~ geno.glm + Townsend + AgeAtRecruit + sex + pc1 + pc2 + pc3 + pc4 + pc5 + pc6 + pc7 + pc8 + pc9 + pc10"))

    }

    fit.bin <- haplo.glm(f, family = binomial,
        data=glm.data, na.action = "na.geno.keep",
        locus.label=label,
        control = haplo.glm.control(haplo.freq.min=.01, 
            em.c = haplo.em.control(min.posterior = 0.1)))

    haplo <- summary(fit.bin)$haplotype
    haplo_name.list <- lapply(1:(length(names(haplo))-1), function(n){
            paste(names(haplo)[n],haplo[,n],sep = "*")
    })
    haplo$name <- Reduce(function(x, y) paste(x, y, sep = "~"), haplo_name.list)

    CI <- signif(as.data.frame(confint.default(fit.bin)), 3)
    var <- signif(as.data.frame(coef(summary(fit.bin))), 3)
    var <- var[!(row.names(var) %in% c("geno.glm.rare")),]
    var <- var[grep("geno",row.names(var)),]
    var_list <- lapply(1:nrow(var), function(i){
        summary_stats <- cbind(haplo[row.names(var)[i],]$name, haplo[row.names(var)[i],]$hap.freq, var[i, 1, drop = FALSE], var[i, 2, drop = FALSE], paste(CI[i,1], CI[i,2], sep = "-") , var[i, 4, drop = FALSE])
        names(summary_stats) <- c("Haplotype", "Haplo_freq", "b", "StdErr", "95%_CI", "p")
        if(FILE == "ukbb" && REGION == "Proxy"){

            summary_stats$b_adjusted <- summary_stats$b * 2
            summary_stats$StdErr_adjusted <- summary_stats$StdErr * 2
            summary_stats$p_adjusted <- summary_stats$p

        }
        return(summary_stats)
    })

    result <- Reduce(rbind, var_list)
    result$cohort <- FILE
    result$N <- nrow(HLA_merged)
    result$N_case <- nrow(subset(glm.data, phenotype == 1))
    result$N_control <- nrow(subset(glm.data, phenotype == 0))
        filename <- paste(label, collapse = "_")
        write_delim(result,paste0("haplotype_",filename,"_",FILE,"_",REGION,".txt"), delim = " ", na = "NA", append = FALSE, col_names = TRUE, quote_escape = "double")

}

haplotype_glm("DQ|DR")
