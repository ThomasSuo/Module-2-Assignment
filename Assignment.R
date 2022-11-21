getwd()
setwd("/Users/zihansuo/InformaticsAssignment")

# read in RNAseq data
rna.data <- read.table("TCGA.COAD.HiSeqV2",sep="\t",head=TRUE,row.names=1)
dim (rna.data)
rna.data [1:3,1:3]
rna.data <- as.matrix (rna.data)

# read in survival data
surv.data <- read.table("COAD_survival.txt", sep="\t", header=T,row.names=1) 
rownames(surv.data)<-gsub(rownames(surv.data), pattern="-", replace=".")

# read in clinical data
clin.data <- read.table("COAD_clinicalMatrix", sep="\t", header=T, row.names=1)
rownames(clin.data)<-gsub(rownames(clin.data), pattern="-", replace=".")

library(survival)
os.time <- surv.data[colnames(rna.data),"OS.time"]
os.event <- as.numeric(surv.data[colnames(rna.data),"OS"])
cesc.os <- Surv(os.time,os.event)
coxph (cesc.os~rna.data [1,])

#Create a univariate model table
results.univariate <- array (NA, c(nrow(rna.data),4))
results.univariate [1:3,]
colnames (results.univariate) <- c("HR", "LCI", "UCI", "PVAL")
rownames (results.univariate) <- rownames (rna.data)
results.univariate [1:3,]
results.univariate <- as.data.frame (results.univariate)

for (i in 1:nrow(rna.data)){
  coxphmodel <- coxph(cesc.os ~ as.numeric(rna.data[i,]))
  results.univariate$HR[i] <- summary(coxphmodel)$coef[1,2]
  results.univariate$LCI[i] <- summary(coxphmodel)$conf.int[1,3]
  results.univariate$UCI[i] <- summary(coxphmodel)$conf.int[1,4]
  results.univariate$PVAL[i] <- summary(coxphmodel)$coef[1,5]
}
results.univariate<-as.data.frame(results.univariate)
results.univariate$FDR<-p.adjust(results.univariate$PVAL,method="fdr")
results.univariate<-results.univariate[order(results.univariate$FDR, decreasing=F),]
results.univariate [1:10,]

#Select potential genes which FDR <= 0.05 and make a table
PotentialGenes <- results.univariate [which (results.univariate$FDR <= 0.05),]
nrow (PotentialGenes)

#Determine the threshold of median expression value in rna.data
summary (rna.data ["CETN1",])
getElement(summary(rna.data ["CETN1",]),"Median")
expval <- array (NA, c(nrow(rna.data), 2))
colnames (expval) <- c("Mean", "Median")
rownames (expval) <- rownames (rna.data)
expval [1:3,]
expval <- as.data.frame (expval)
for (i in 1:nrow(rna.data)){
  expval$Median[i] <- getElement(summary(rna.data [i,]),"Median")
  expval$Mean[i] <- getElement(summary (rna.data[i,]), "Mean")
}
expval [1:3,]

#Do some statistics in expval median and mean
summary (expval$Median)
summary (expval$Mean)
boxplot (expval$Mean, expval$Median, names = c("Mean", "Median"), col = c("yellow","blue"), ylab = "Expression value (log2(x+1))")

#Make a table to collect HR, Expval$Median, and Expval$Mean, order them in the descending order
ExpvalMedian <- array (NA, c(nrow(PotentialGenes),3))
rownames (ExpvalMedian) <- rownames (PotentialGenes)
colnames (ExpvalMedian) <- c("HR","Median","Mean")
ExpvalMedian <- as.data.frame (ExpvalMedian)
class(PotentialGenes)
PotentialGenes$HR[1]
getElement(summary(rna.data[rownames (ExpvalMedian)[1],]), "Median")
for (i in 1:nrow (PotentialGenes)){
  ExpvalMedian$HR[i] <- PotentialGenes$HR[i]
  ExpvalMedian$Median[i] <- getElement(summary(rna.data[rownames(ExpvalMedian)[i],]), "Median")
  ExpvalMedian$Mean[i] <- getElement(summary (rna.data[rownames(ExpvalMedian)[i],]), "Mean")
}
ExpvalMedian <- ExpvalMedian[order(ExpvalMedian$Median, decreasing = T),]
ExpvalMedian[1:10,]

#We found that "ZNF248" has no research in prognostic or survival before, so that should be a new biomarker candidate
summary (rna.data ["ZNF248",])

#Fit a multivariate model. Factors: Age, Gender, Obesity, colon polyps
clin.data <- clin.data [colnames (rna.data),]
clin.data [1:10,]
colnames (clin.data)
clin.data$colon_polyps_present
clin.data$year_of_initial_pathologic_diagnosis
clin.data$age_at_initial_pathologic_diagnosis
