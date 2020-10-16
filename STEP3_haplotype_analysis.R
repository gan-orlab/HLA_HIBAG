#module load gcc/7.3 r-bundle-bioconductor/3.9

.libPaths(c("~/runs/eyu8/library/haplo.stats",
"/cvmfs/soft.computecanada.ca/easybuild/software/2017/avx512/Compiler/gcc7.3/r-bundle-bioconductor/3.9",
"/cvmfs/soft.computecanada.ca/easybuild/software/2017/avx512/Compiler/gcc7.3/r/3.6.0/lib64/R/library"))

library(haplo.stats)
library(readr)

seed <- c(17, 53, 1, 40, 37, 0, 62, 56, 5, 52, 12, 1)
set.seed(seed)

args <- commandArgs(trailingOnly = TRUE)

#args <- c("MCGILL","Euro")

FILE <- args[1]
REGION <- args[2]

A <- read.csv(file=paste0("txt_data/",FILE,"/HLA-A_",FILE,"_",REGION,".logit"), sep=",",stringsAsFactors=FALSE)
B <- read.csv(file=paste0("txt_data/",FILE,"/HLA-B_",FILE,"_",REGION,".logit"), sep=",",stringsAsFactors=FALSE)
C <- read.csv(file=paste0("txt_data/",FILE,"/HLA-C_",FILE,"_",REGION,".logit"), sep=",",stringsAsFactors=FALSE)
DPB1 <- read.csv(file=paste0("txt_data/",FILE,"/HLA-DPB1_",FILE,"_",REGION,".logit"), sep=",",stringsAsFactors=FALSE)
DQA1 <- read.csv(file=paste0("txt_data/",FILE,"/HLA-DQA1_",FILE,"_",REGION,".logit"), sep=",",stringsAsFactors=FALSE)
DQB1 <- read.csv(file=paste0("txt_data/",FILE,"/HLA-DQB1_",FILE,"_",REGION,".logit"), sep=",",stringsAsFactors=FALSE)
DRB1 <- read.csv(file=paste0("txt_data/",FILE,"/HLA-DRB1_",FILE,"_",REGION,".logit"), sep=",",stringsAsFactors=FALSE)

covar <- A[, !(names(A) %in% c("allele1","allele2","prob","matching"))]
covar$phenotype <- covar$phenotype - 1

HLA <- list(A=A,B=B,C=C,DPB1=DPB1,DQA1=DQA1,DQB1=DQB1,DRB1=DRB1)


poorSamples <- lapply(HLA,function(x) x[x$prob<0.5,]$ID)
poorSamples <- unlist(poorSamples, use.names=FALSE)

filtered_HLA <- lapply(HLA,function(x) x[!(x$ID %in% poorSamples), c("ID","allele1","allele2")])



HLA.list <- lapply(names(filtered_HLA),function(gene){
	HLA_gene <- filtered_HLA[[gene]]
	names(HLA_gene) <- c("ID",paste(gene,c("a1","a2"),sep = "."))
	return(HLA_gene)
	})
HLA_merged <- merge(covar, Reduce(merge, HLA.list))
HLA_covar <- HLA_merged[, !(names(HLA_merged) %in% grep(".a",names(HLA_merged)))]


haplotype_glm <- function(pattern){

    geno <- HLA_merged[,grep(pattern,names(HLA_merged))]
    label <- names(HLA)[grep(pattern,names(HLA))]

	geno.glm <- setupGeno(geno, miss.val=c(0,NA), locus.label=label)

	glm.data <- cbind(HLA_covar,data.frame(geno.glm))

	fit.bin <- haplo.glm(phenotype ~ geno.glm + sex + age + pc1 + pc2 + pc3 + pc4 + pc5 + pc6 + pc7 + pc8 + pc9 + pc10, family = binomial,
		data=glm.data, na.action = "na.geno.keep",
		locus.label=label,
		control = haplo.glm.control(em.c = haplo.em.control(min.posterior = 0.2)))


	haplo <- summary(fit.bin)$haplotype
	haplo_name.list <- lapply(1:(length(names(haplo))-1), function(n){
		paste(names(haplo)[n],haplo[,n],sep = "_")
	})
	haplo$name <- Reduce(function(x, y) paste(x, y, sep = "_"), haplo_name.list)

	CI <- signif(as.data.frame(confint.default(fit.bin)), 3)
    var <- signif(as.data.frame(coef(summary(fit.bin))), 3)
    var <- var[!(row.names(var) %in% c("geno.glm.rare")),]
    var <- var[grep("geno",row.names(var)),]
    var_list <- lapply(1:nrow(var), function(i){
        summary_stats <- cbind(haplo[row.names(var)[i],]$name, haplo[row.names(var)[i],]$hap.freq, var[i, 1, drop = FALSE], var[i, 2, drop = FALSE], paste(CI[i,1], CI[i,2], sep = "-") , var[i, 4, drop = FALSE])
        names(summary_stats) <- c("Haplotype", "Haplo_freq", "b", "StdErr", "95%_CI", "p")
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

haplotype_glm("DQ")
haplotype_glm("DQA|DR")
haplotype_glm("DQB|DR")
haplotype_glm("DQ|DR")
