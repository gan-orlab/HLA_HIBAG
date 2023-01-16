library(HIBAG)
library(parallel)

#args[1] is gene args[2] is reference args[3] is suffix
args <- commandArgs(trailingOnly = TRUE)

gene <- args[1]
FILE <- args[2]
SUFFIX <- args[3]
REF <- args[4]
core <- args[5]

sink(paste("log/",SUFFIX,"/HLA-",gene,"_",SUFFIX,".log",sep=""), split = TRUE)


#parallel for snow
cl <- makeCluster(as.numeric(core))

model.list <- get(load(REF))

#bim files
yourgeno <- hlaBED2Geno(bed.fn=paste(FILE,".bed",sep=""), fam.fn=paste(FILE,".fam",sep=""), bim.fn=paste(FILE,".bim",sep=""))

message("Imputation of HLA ",SUFFIX)
model <- hlaModelFromObj(model.list)

pred.guess <- predict(model, yourgeno,
type="response+prob",match.type="Position",cl=cl)

summary(pred.guess)

write.csv(pred.guess$value,paste(SUFFIX,"/HLA-",gene,"_",SUFFIX,".csv",sep=""),row.names=FALSE)


