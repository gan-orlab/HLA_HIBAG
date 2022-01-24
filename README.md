# Fine mapping of the *HLA* locus in Parkinson’s disease in Europeans

## WORKFLOW 
#### [0. Extract SNPs from HIBAG reference panel](#0)
#### [1. Pre-processing of imputed data](#1)
#### [2. HLA imputation](#2)
#### [3.1 Statistical analysis - HLA allele](#3.1)
#### [3.2 Statistical analysis - HLA haplotype](#3.2)
#### [3.3 Statistical analysis - HLA amino acid](#3.3)
#### [4. Meta-analysis](#4)

## Requirements

- R packages
	- HIBAG
	- data.table
	- parallel
	- readr
	- haplo.stats
 

<a id="0"></a>
## 0. Extract SNPs from HIBAG reference panel
The SNPs and alleles used in a HIBAG reference panel are first extracted to ensure that the imputed SNPs matches with those on the reference panel. Here, duplicated SNPs are filtered into uniques SNPs.
```bash
RDATA=$1
REF=$2
module load gcc/7.3.0 r-bundle-bioconductor/3.9

mkdir tmp
Rscript script/extract_snp.R $RDATA

cat tmp/*position* > tmp/HLA_snp_position.csv
cat tmp/*allele* > tmp/HLA_snp_allele.csv
paste tmp/HLA_snp_position.csv tmp/HLA_snp_allele.csv | sed 's/"//g'|tr '/' '\t' > tmp/dup.txt

sort -u tmp/dup.txt > tmp/unique_HLA.csv

awk '{print "6:"$1}' tmp/unique_HLA.csv > tmp/HLA_snp.csv

mv tmp/ RData/csv/$REF/
```
A small R script is used to read the SNPs from the reference panels. 
```R
#Extract snp and allele from gene
library(HIBAG)
args <- commandArgs(trailingOnly = TRUE)

model.list <- get(load(args[1]))

gene <- names(model.list)
for (n in gene){

    write.csv(model.list[[n]]$snp.position,paste("tmp/HLA-",n,"_snp_position.csv",sep=""),col.names=FALSE,row.names=FALSE)
    write.csv(model.list[[n]]$snp.allele,paste("tmp/HLA-",n,"_snp_allele.csv",sep=""),col.names=FALSE,row.names=FALSE)
```

<a id="1"></a>
## 1. Pre-processing of imputed data
Individual-level data was first imputed using 1000G phase 3 version 5. This step is used to select  hard calls (*r*>0.8) SNPs from the reference panel and add sex and phenotype labels removed during imputation. The R script called exc.R will match and flip alleles accordingly. 
```bash
#Extract SNPs hard calls in reference panel from VCF; add sex, pheno; remove duplicates

FILE=$1
REF=$2
ANCESTRY=$3


awk '{print 0" "$1"_"$2" "$6}' ~/runs/eyu8/QC_NIH/QC_$FILE/*/data_here/FILTERED.fam > txt_data/pheno/$FILE.pheno
awk '{print 0" "$1"_"$2" "$5}' ~/runs/eyu8/QC_NIH/QC_$FILE/*/data_here/FILTERED.fam > txt_data/sex/$FILE.sex

plink --const-fid --vcf  vcfFiles/$FILE.vcf.gz  --keep MERGED_unrelated.fam --make-bed --out bfile/$FILE

awk '{print $2}' bfile/$FILE.bim | sort | uniq -d  > txt_data/dup/$FILE.dup

plink --bfile bfile/$FILE --exclude txt_data/dup/$FILE.dup --make-bed --out bfile/unique_$FILE

#Extract hard calls
plink --bfile bfile/unique_$FILE --qual-scores vcfFiles/$FILE.info 7 1 1 --qual-threshold 0.8 --make-bed --out bfile/$FILE'_hardcall'

cut -f 3-6 bfile/$FILE'_hardcall'.bim > bfile/temp.bim
cut -f 1-2 bfile/$FILE'_hardcall'.bim | cut -d ":" -f 1-2 > bfile/temp2.bim
paste -d "\t" bfile/temp2.bim bfile/temp.bim > bfile/$FILE'_hardcall'.bim

rm bfile/temp*

plink --bfile bfile/$FILE'_hardcall' --extract RData/csv/$REF/HLA_snp.csv --update-sex txt_data/sex/$FILE.sex --pheno txt_data/pheno/$FILE.pheno --make-bed --out pre_exc

Rscript exc.R RData/csv/$REF/unique_HLA.csv pre_exc.bim

mv a pre_exc.bim

plink --bfile pre_exc --extract RData/csv/$REF/HLA_snp.csv --make-bed --out bfile/final_$FILE'_'$ANCESTRY

rm pre_exc*

```
<a id="2"></a>
## 2. HLA imputation
HLA allele imputation is ran using the R package HIBAG with imputed plink genotype bfile as input along with a HIBAG reference panel. 
```bash
Rscript script/HIBAG_script.R $gene $bfile  $SUFFIX  $REF  $core
```
Load necessary packages, create logs, setup multiple CPUs
```R
library(HIBAG)
library(parallel)

args <- commandArgs(trailingOnly = TRUE)

gene <- args[1]
FILE <- args[2]
SUFFIX <- args[3]
REF <- args[4]
core <- args[5]

#Output logs to this file
sink(paste("log/",SUFFIX,"/HLA-",allele,"_",SUFFIX,".log",sep=""), split = TRUE)


#parallel for snow
cl <- makeCluster(as.numeric(core))
```
Load reference panel, load data using hlaBED2Geno function
```R
model.list <- get(load(REF))

#bim files
yourgeno <- hlaBED2Geno(bed.fn=paste(FILE,".bed",sep=""), fam.fn=paste(FILE,".fam",sep=""), bim.fn=paste(FILE,".bim",sep=""))
```
Load HLA gene of interest (e.g. A) and impute alleles with the predict function based on positions of SNPs.
```R
message("Imputation of HLA-",SUFFIX)
model <- hlaModelFromObj(model.list[[gene]])

pred.guess <- predict(model, yourgeno, type="response+prob",match.type="Position",cl=cl)

write.csv(pred.guess$value,paste(SUFFIX,"/HLA-",allele,"_",SUFFIX,".csv",sep=""),row.names=FALSE)

```
<a id="3.1"></a>
## 3.1 Statistical analysis - HLA allele
This will be an example of the testing the assication of the dominant effect of HLA allele on Parkinson's disease with logistic regression . Other covariates can be added as necessary.

Load the HLA alleles and other metadata (e.g. age, sex, phenotype, principal components)
```R
library(readr)
library(data.table)

#arg 1 is FILE arg2 is REGION
args <- commandArgs(trailingOnly = TRUE)

FILE <- args[1]
REGION <- args[2]

#load data depending if its UKB or not
if(FILE != "ukbb"){

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

    PD <- as.data.frame(fread("ukbb_PD_covar.txt"))
    control_PD <- as.data.frame(fread("ukbb_control_PD_covar.txt"))

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
    DRB3 <- as.data.frame(fread("ukbb_imp_HLA_DRB3/HLA-DRB3_ukbb_imp_HLA_DRB3.csv"))
    DRB4 <- as.data.frame(fread("ukbb_imp_HLA_DRB4/HLA-DRB4_ukbb_imp_HLA_DRB4.csv"))
    DRB5 <- as.data.frame(fread("ukbb_imp_HLA_DRB5/HLA-DRB5_ukbb_imp_HLA_DRB5.csv"))


} else if(FILE == "ukbb" && REGION == "Proxy"){

    Proxy <- as.data.frame(fread("ukbb_proxy_covar.txt"))
    control_Proxy <- as.data.frame(fread("ukbb_control_proxy_covar.txt"))

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
    DRB3 <- as.data.frame(fread("ukbb_imp_HLA_DRB3/HLA-DRB3_ukbb_imp_HLA_DRB3.csv"))
    DRB4 <- as.data.frame(fread("ukbb_imp_HLA_DRB4/HLA-DRB4_ukbb_imp_HLA_DRB4.csv"))
    DRB5 <- as.data.frame(fread("ukbb_imp_HLA_DRB5/HLA-DRB5_ukbb_imp_HLA_DRB5.csv"))
}

#Merge allele with covariates
HLA <- list(A=merge(covar,A),B=merge(covar,B),C=merge(covar,C),
    DPB1=merge(covar,DPB1),DQA1=merge(covar,DQA1),DQB1=merge(covar,DQB1),DRB1=merge(covar,DRB1),
    DRB3=merge(covar,DRB3),DRB4=merge(covar,DRB4),DRB5=merge(covar,DRB5))
```
Example data for variable HLA[[A]]:
|`sample.id`|sex|phenotype|age| pc1 | pc2 | pc3 | pc4 | pc5 | allele1 | allele2| prob| matching|
|:-------------:|:---:|----------|:---------:|:---------:|:---------:|:---------:|:---------:|:---------:|:---------:|:---------:|:---------:|:---------:|:---------:|
|123456|1|2|65|0.04|0.01| 0.01| 0.01| 0.006|01:01|01:02| 0.99| 1.65e-5|
|123457|2|1|62|-0.01|0.02| 0.03| 0.01| 0.006|03:01|01:02| 0.97| 6.65e-4|

Filter out samples with 2 alleles with posterior probability less than 0.5
```R
#Exclude samples with 2 alleles with posterior probability less than 0.5
poorSamples <- lapply(HLA,function(x) x[x$prob<0.5,]$sample.id)
poorSamples <- unlist(poorSamples, use.names=FALSE)
poorSamples <- poorSamples[duplicated(poorSamples)]
filtered_HLA <- lapply(HLA,function(x) x[!(x$sample.id %in% poorSamples),])
```
Perform analysis for each gene (e.g. HLA[[A]]) independently. Filter alleles with prob < 0.5.
```R
for(i in 1:10){

    message(paste0("Analysing HLA-",names(HLA)[i]))
    #Exclude allele where posterior probability less than 0.5
    filtered_HLA_gene <- filtered_HLA[[i]][filtered_HLA[[i]]$prob > 0.5,]
    #Extract alleles
    HLA_allele1 <- filtered_HLA_gene[, c("sample.id","allele1")]
    HLA_allele2 <- filtered_HLA_gene[, c("sample.id","allele2")]
    names(HLA_allele1) <- c("sample.id","allele")
    names(HLA_allele2) <- c("sample.id","allele")
    #Marge allele into 1 column
    HLA_allele <- rbind(HLA_allele1, HLA_allele2)
    HLA_count <- as.data.frame.matrix(table(HLA_allele))
    HLA_allele_name <- names(HLA_count)
    HLA_count$sample.id <- rownames(HLA_count)
```

Perform logistic regression 
```R
    #Change homozygote to heterozygote for a dominant model
    HLA_count[HLA_count == 2] <- 1
    HLA_meta_data <- filtered_HLA_gene[, !(names(filtered_HLA_gene) %in% c("allele1", "allele2", "prob", "matching"))]
    #Merge back alleles
    HLA_glm <- merge(HLA_meta_data, HLA_count)
    #Set sex, phenotype as factor; exclude alleles not present in controls
    HLA_table <- lapply(HLA_allele_name, function(allele){
        HLA_glm$phenotype[HLA_glm$phenotype == -9] <- NA
        HLA_glm$phenotype <- as.factor(HLA_glm$phenotype - 1)
        HLA_glm$sex <- as.factor(HLA_glm$sex)
        names(HLA_glm)[names(HLA_glm) == allele] <- "h"
        HLA_glm$h <- as.numeric(HLA_glm$h)
        if(length(levels(as.factor(HLA_glm$h))) == 1 ){
            return(NULL)
        }
        if(length(levels(as.factor(HLA_glm$phenotype))) == 1){
            return(NULL)
        }
        if(FILE != "ukbb"){

            f <- formula(paste("phenotype ~ h + age + sex + pc1 + pc2 + pc3 + pc4 + pc5 + pc6 + pc7 + pc8 + pc9 + pc10"))

        } else if(FILE == "ukbb"){

            f <- formula(paste("phenotype ~ h + Townsend + AgeAtRecruit + sex + pc1 + pc2 + pc3 + pc4 + pc5 + pc6 + pc7 + pc8 + pc9 + pc10"))

        }

        HLA_result <- glm(f,family = binomial, data = HLA_glm)
        #Get frequency and count
        HLA_carrier_freq <-  table(HLA_glm[, c("h", "phenotype")])
        HLA_freq <- signif((HLA_carrier_freq[2,1] + HLA_carrier_freq[2,2])/sum(HLA_carrier_freq),3)
        HLA_carrier_freq_case <- signif(HLA_carrier_freq[2,2]/(HLA_carrier_freq[1,2] + HLA_carrier_freq[2,2]),3)
        HLA_carrier_freq_control <- signif(HLA_carrier_freq[2,1]/(HLA_carrier_freq[1,1] + HLA_carrier_freq[2,1]),3)
        HLA_ncase <- HLA_carrier_freq[1,2] + HLA_carrier_freq[2,2]
        HLA_ncontrol <- HLA_carrier_freq[1,1] + HLA_carrier_freq[2,1]
        HLA_ntotal <- sum(HLA_carrier_freq)
        #Get statistics (beta, se, pvalue)
        var <- signif(as.data.frame(coef(summary(HLA_result))),3)
        if(!grepl("h",row.names(var)[2])){
            return(NULL)
        }
        var_list <- lapply(2:nrow(var), function(i){
            summary_stats <- cbind(var[i, 1, drop = FALSE], var[i, 2, drop = FALSE], var[i, 4, drop = FALSE])
            names(summary_stats) <- c("b", "StdErr", "p")
            names(summary_stats) <- paste(row.names(summary_stats), names(summary_stats), sep = "_")
            return(summary_stats)
        })
```
Adjust beta, standard error for proxy-cases
```R
        if(FILE == "ukbb" && REGION == "Proxy"){

            var_list <- append(var_list[1], var_list)
            var_list[[2]] <- var_list[[2]] * c(2,2,1)
            names(var_list[[2]]) <- c("h_b_adjusted", "h_StdErr_adjusted", "h_p_adjusted")
        }
        return(cbind(allele, HLA_ntotal, HLA_ncase, HLA_ncontrol, HLA_freq, HLA_carrier_freq_case, HLA_carrier_freq_control, Reduce(cbind, var_list)))})
    result <- Reduce(rbind, HLA_table)
    write_delim(result,paste0("HLA-",names(HLA)[i],"_",FILE,"_",REGION,".txt"), delim = " ", na = "NA", append = FALSE, col_names = TRUE, quote_escape = "double")
}
```
<a id="3.2"></a>
## 3.2 Statistical analysis - HLA haplotype
This will be an example of the testing the assication of the dominant effect of HLA haplotypes on Parkinson's disease with logistic regression . Other covariates can be added as necessary.

Load the HLA alleles and other metadata (e.g. age, sex, phenotype, principal components)
```R
library(haplo.stats)
library(readr)

seed <- c(17, 53, 1, 40, 37, 0, 62, 56, 5, 52, 12, 1)
set.seed(seed)

args <- commandArgs(trailingOnly = TRUE)

FILE <- args[1]
REGION <- args[2]

...

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
```
Example data for variable HLA_covar:
|ID|sex|phenotype|age| pc1 | pc2 | pc3 | pc4 | pc5 | A.a1 | A.a2| B.a1| B.a2|
|:-------------:|:---:|----------|:---------:|:---------:|:---------:|:---------:|:---------:|:---------:|:---------:|:---------:|:---------:|:---------:|:---------:|
|123456|1|2|65|0.04|0.01| 0.01| 0.01| 0.006|01:01|01:02| 40:02| 51:01|
|123457|2|1|62|-0.01|0.02| 0.03| 0.01| 0.006|03:01|01:02| 18:01| 56:01|

Perform haplotype estimation and logistic regression 

```R
haplotype_glm <- function(pattern){

    geno <- HLA_merged[,grep(pattern,names(HLA_merged))]
    label <- names(HLA)[grep(pattern,names(HLA))]

	geno.glm <- setupGeno(geno, miss.val=c(0,NA), locus.label=label)

	glm.data <- cbind(HLA_covar,data.frame(geno.glm))
	
	#Regression adjusted for covariates
	fit.bin <- haplo.glm(phenotype ~ geno.glm + sex + age + pc1 + pc2 + pc3 + pc4 + pc5 + pc6 + pc7 + pc8 + pc9 + pc10, family = binomial,
		data=glm.data, na.action = "na.geno.keep",
		locus.label=label,
		control = haplo.glm.control(em.c = haplo.em.control(min.posterior = 0.2)))

	#Rearrange data
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
```
<a id="3.3"></a>
## 3.3 Statistical analysis - HLA amino acid
This will be an example of the testing the assication of the dominant effect of HLA amino acid on Parkinson's disease with logistic regression . Other covariates can be added as necessary. Note: HIBAG only provides amino acid for *HLA-A*, *HLA-B*, *HLA-C*, *HLA-DPB1*, *HLA-DQA1*, *HLA-DQB1*, *HLA-DRB1*.

Load the HLA alleles and other metadata (e.g. age, sex, phenotype, principal components)
```R

library(data.table)
library(HIBAG)
library(readr)

#arg 1 is FILE arg2 is REGION
args <- commandArgs(trailingOnly = TRUE)

FILE <- args[1]
REGION <- args[2]

...

HLA <- list(A=merge(covar,A),B=merge(covar,B),C=merge(covar,C),
            DPB1=merge(covar,DPB1),DQA1=merge(covar,DQA1),DQB1=merge(covar,DQB1),DRB1=merge(covar,DRB1))
```
Filter out samples with 2 alleles with posterior probability less than 0.5
```R
poorSamples <- lapply(HLA,function(x) x[x$prob<0.5,]$"sample.id")
poorSamples <- unlist(poorSamples, use.names=FALSE)
poorSamples <- poorSamples[duplicated(poorSamples)]
filtered_HLA <- lapply(HLA,function(x) x[!(x$"sample.id" %in% poorSamples),])
```
Perform analysis for each gene (e.g. HLA[[A]]) independently. Filter amino acid with prob < 0.5.
```R
for(i in 1:7){
    message(paste0("Analysing HLA-",names(HLA)[i]))
	hla_group <- filtered_HLA[[i]]
    hla_group$phenotype[hla_group$phenotype == -9] <- NA
    hla_group$phenotype <- as.factor(hla_group$phenotype - 1)
    #Reformat allele
    hla <- hlaAllele(hla_group$"sample.id", H1=hla_group$allele1, H2=hla_group$allele2, locus=names(HLA[i]), assembly="hg19", prob=hla_group$prob)
    name <- names(HLA[i])
    #Predict amino acid
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
```
Perform logistic regression for each amino acid. Relabel reference amino acid, ambiguous and indels. 
```R
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
	    HLA_glm <- merge(HLA_meta_data, HLA_count)
	    HLA_table <- lapply(HLA_allele_name, function(allele){
	        HLA_glm$sex <- as.factor(HLA_glm$sex)
	        HLA_glm[,allele] <- as.numeric(HLA_glm[,allele])

	        if(length(levels(as.factor(HLA_glm[,allele]))) == 1){
	        	return(NULL)
	        }
	        f <- formula(paste("phenotype ~ ",allele," + age + sex + pc1 + pc2 + pc3 + pc4 + pc5 + pc6 + pc7 + pc8 + pc9 + pc10"))
	        HLA_result <- glm(f,family = binomial, data = HLA_glm)
	        HLA_carrier_freq <-  table(HLA_glm[, c(allele, "phenotype")])
	        HLA_freq <- signif((HLA_carrier_freq[2,1] + HLA_carrier_freq[2,2])/sum(HLA_carrier_freq),3)
	        HLA_carrier_freq_case <- signif(HLA_carrier_freq[2,2]/(HLA_carrier_freq[1,2] + HLA_carrier_freq[2,2]),3)
	        HLA_carrier_freq_control <- signif(HLA_carrier_freq[2,1]/(HLA_carrier_freq[1,1] + HLA_carrier_freq[2,1]),3)
	        HLA_ncase <- HLA_carrier_freq[1,2] + HLA_carrier_freq[2,2]
	        HLA_ncontrol <- HLA_carrier_freq[1,1] + HLA_carrier_freq[2,1]
	        HLA_ntotal <- sum(HLA_carrier_freq)
	        CI <- signif(as.data.frame(confint.default(HLA_result)),3)
	        var <- signif(as.data.frame(coef(summary(HLA_result))),3)
            if(!grepl(allele,row.names(var)[2])){
                return(NULL)
            }
	        var_list <- lapply(2:nrow(var), function(i){
	            summary_stats <- cbind(var[i, 1, drop = FALSE], var[i, 2, drop = FALSE], paste(CI[i,1], CI[i,2], sep = "-") , var[i, 4, drop = FALSE])
	            names(summary_stats) <- c("b", "StdErr", "95%_CI", "p")
	            if(i == 2){
	            	row.names(summary_stats) <- "aa"
	        	}
	            names(summary_stats) <- paste(row.names(summary_stats), names(summary_stats), sep = "_")
	            return(summary_stats)
	            })
	        return(cbind(Pos = paste0("Pos",pos + increment,allele), Ref = "z", Alt = allele, HLA_ntotal, HLA_ncase, HLA_ncontrol, HLA_freq, HLA_carrier_freq_case, HLA_carrier_freq_control, Reduce(cbind, var_list)))})
	    result <- Reduce(rbind, HLA_table)
	    return(result)
	    })
    
    hla_result <- Reduce(rbind, hla_aa_result)
    write_delim(hla_result,paste0("HLA-",names(HLA)[i],"_",FILE,"_",REGION,"_aa.txt"), delim = " ", na = "NA", append = FALSE, col_names = TRUE, quote_escape = "double")
}

```
<a id="4"></a>
## 4. Meta-analysis
Fixed-effect meta-analysis was performed using METAL with inverse variance weighting

```bash
allel=$1

module load metal

echo "MARKER allele" > meta.script
echo "PVALUE h_p" >> meta.script
echo "EFFECT h_b" >> meta.script
echo "WEIGHT HLA_ntotal" >> meta.script
echo "SCHEME STDERR" >> meta.script
echo "STDERR h_StdErr" >> meta.script
echo "FREQLABEL HLA_freq" >> meta.script
echo "AVERAGEFREQ ON" >> meta.script
echo "MINMAXFREQ ON" >> meta.script
echo "PROCESS      result/HLA-${allel}_NIND_Euro.txt" >> meta.script
echo "PROCESS      result/HLA-${allel}_IPDGC_Euro.txt" >> meta.script
echo "PROCESS      result/HLA-${allel}_APDGC_Euro.txt" >> meta.script
echo "PROCESS      result/HLA-${allel}_Oslo_Euro.txt" >> meta.script
echo "PROCESS      result/HLA-${allel}_MCGILL_Euro.txt" >> meta.script
echo "PROCESS      result/HLA-${allel}_PPMI_Euro.txt" >> meta.script
echo "PROCESS      result/HLA-${allel}_NGRC_Euro.txt" >> meta.script
echo "PROCESS      ukbb/HLA-${allel}_ukbb_PD.txt" >> meta.script
echo "EFFECT h_b_adjusted" >> meta.script
echo "STDERR h_StdErr_adjusted" >> meta.script
echo "PROCESS      ukbb/HLA-${allel}_ukbb_Proxy.txt" >> meta.script
echo "OUTFILE	METAL-HLA-${allel} .tbl" >> meta.script
echo "ANALYZE HETEROGENEITY" >> meta.script
echo "QUIT" >> meta.script


metal meta.script

```
