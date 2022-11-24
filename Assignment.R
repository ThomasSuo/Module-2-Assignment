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

#Fit a multivariate model. Factors: Age, , anatomic, colon polyps
clin.data <- clin.data [colnames (rna.data),]
clin.data [1:10,]
colnames (clin.data)
clin.data$colon_polyps_present
clin.data$year_of_initial_pathologic_diagnosis
clin.data$age_at_initial_pathologic_diagnosis
#Factor1: Age
age <- as.numeric (clin.data$age_at_initial_pathologic_diagnosis)
#Factor2: Polyps
polyps <- clin.data$colon_polyps_present
polyps [which (clin.data$colon_polyps_present=="YES")] <- 1
polyps [which (clin.data$colon_polyps_present=="NO")] <- 0
polyps
#Factor3: Lymph node counts 
lymph.node <- as.numeric(clin.data$lymph_node_examined_count)
lymph.node
clin.data$height
#Factor4: Obesity
height <- clin.data$height /100
clin.data$weight
weight <- clin.data$weight
bmi <- weight / height**2
bmi <- as.numeric(bmi)
##high bmi <- 1; normal bmi <- 0
bmi [which (bmi < 25)] <- 0
bmi [which (bmi >= 25)] <- 1
bmi

#whether the factors associate with survival data??
#1. Age
summary (coxph(cesc.os~age))$coef
summary (coxph (cesc.os~polyps))$coef
summary (coxph (cesc.os~lymph.node))$coef
summary (coxph (cesc.os~bmi))$coef

#Make a multivariate table
results.multivariate <- array (NA, c(nrow (rna.data),4))
colnames (results.multivariate) <- c ("HR", "LCI", "UCI", "PVAL")
rownames (results.multivariate) <- rownames (rna.data)
results.multivariate <- as.data.frame (results.multivariate)
for(i in 1:nrow(rna.data)){
  coxphmodel <- coxph(cesc.os ~ rna.data[i,]+age+polyps+lymph.node)
  results.multivariate$HR[i] <- summary(coxphmodel)$coef[1,2]
  results.multivariate$LCI[i] <- summary(coxphmodel)$conf.int[1,3]
  results.multivariate$UCI[i] <- summary(coxphmodel)$conf.int[1,4]
  results.multivariate$PVAL[i] <- summary(coxphmodel)$coef[1,5]
}
results.multivariate<-as.data.frame(results.multivariate)
results.multivariate$FDR<-p.adjust(results.multivariate$PVAL,method="fdr")
results.multivariate<-results.multivariate[order(results.multivariate$FDR, decreasing=F),]
results.multivariate[1:10,]

#Select Genes which FDR <= 0,05
PotentialGenesMulti <- results.multivariate [which (results.multivariate$FDR <= 0.05),]
nrow (PotentialGenesMulti)

#Make a table for multivariate to collect HR, Expval$Median, and Expval$Mean, order them in the descending order
ExpvalMedianMulti <- array (NA, c(nrow(PotentialGenesMulti),3))
rownames (ExpvalMedianMulti) <- rownames (PotentialGenesMulti)
colnames (ExpvalMedianMulti) <- c("HR","Median","Mean")
ExpvalMedianMulti <- as.data.frame (ExpvalMedianMulti)
ExpvalMedianMulti[1:10,]
ExpvalMedianMulti <- as.data.frame (ExpvalMedianMulti)
class(PotentialGenesMulti)
PotentialGenesMulti$HR[1]
getElement(summary(rna.data[rownames (ExpvalMedianMulti)[1],]), "Median")
for (i in 1:nrow (PotentialGenesMulti)){
  ExpvalMedianMulti$HR[i] <- PotentialGenesMulti$HR[i]
  ExpvalMedianMulti$Median[i] <- getElement(summary(rna.data[rownames(ExpvalMedianMulti)[i],]), "Median")
  ExpvalMedianMulti$Mean[i] <- getElement(summary (rna.data[rownames(ExpvalMedianMulti)[i],]), "Mean")
}
ExpvalMedianMulti <- ExpvalMedianMulti[order(ExpvalMedianMulti$Median, decreasing = T),]
head (ExpvalMedianMulti)
nrow (ExpvalMedian)
nrow (ExpvalMedianMulti)

#Find same genes in these two tables which FDR <= 0.05
overlap <- intersect(rownames(ExpvalMedian),rownames(ExpvalMedianMulti))
gene.summary <- array (NA, c(length(overlap), 5))
rownames (gene.summary) <- overlap
colnames (gene.summary) <- c("ExpVal.Median", "HR.Uni", "pval.Uni", "HR.Multi", "pval.Multi")
gene.summary <- as.data.frame(gene.summary)

for (i in 1:length(overlap)){
  gene.summary$ExpVal.Median[i] <- expval[rownames (gene.summary)[i],]$Median
  gene.summary$HR.Uni[i] <- ExpvalMedian[rownames (gene.summary)[i],]$HR
  gene.summary$pval.Uni[i] <- results.univariate [rownames(gene.summary)[i],]$FDR
  gene.summary$HR.Multi[i] <- ExpvalMedianMulti[rownames (gene.summary)[i],]$HR
  gene.summary$pval.Multi[i] <- results.multivariate [rownames(gene.summary)[i],]$FDR
}
gene.summary [1:10,]

#Finally, according to multivariate results. ANKDD1A was selected for the biomarker study.
gene.summary["ANKDD1A",]
results.univariate["ANKDD1A",]
results.multivariate["ANKDD1A",]
summary (rna.data["ANKDD1A",])

#Stratify patients ANKDD1A expression and draw KM plot
png ("COAD_OS_byANKDD1A.png", width = 6, height = 6, units = 'in', res = 300)
ANKDD1A.high <- as.numeric (rna.data ["ANKDD1A",] > median (rna.data["ANKDD1A",]))
plot (survfit(cesc.os~ANKDD1A.high), col = c("black","red"), lwd =2, mark.time = TRUE, xlab = "OS time (days)", ylab = "proportion surviving", mean = "HIST2H2BE in Colon Cancer")
legend("topright",legend=c("ANKDD1A-high","ANKDD1A-low"),col=c("red","black"),lwd=2)
text(3000,0.1,"HR=1.82 (95%CI 1.37-2.42)")
dev.off()

#DNA Methylation Detection
##Loading methylation annotation download from server
meth.annot <- readRDS("annot450k.rds")
dim(meth.annot)
##Loading DNA methylation data download from CTGA
meth.data <- read.table ("TCGA.COAD.HumanMethylation450", sep = "\t", header = T, row.names = 1)
dim(meth.data)
dim(rna.data)
rna.data2<-rna.data[,which(is.element(colnames(rna.data),colnames(meth.data)))]
meth.data2<-meth.data[,which(is.element(colnames(meth.data),colnames(rna.data2)))]
surv.data2<-surv.data[which(is.element(rownames(surv.data),colnames(rna.data2))),]

#Arrange them into same order
rna.data2<-as.matrix(rna.data2[,order(colnames(rna.data2))])
meth.data2<-as.matrix(meth.data2[,order(colnames(meth.data2))])
surv.data2<-as.data.frame(surv.data2)
nrow(surv.data2)
ncol(rna.data2)
rna.data2[1:3,1:4]
meth.data2[1:3,1:4]
surv.data2[1:3,1:4]

meth.ankdd1a <- rownames(meth.annot[which (meth.annot$UCSC_RefGene_Name == "ANKDD1A"),])
meth.ankdd1a

meth.annot [meth.ankdd1a,]
meth.data.ankdd1a <- meth.data2[meth.ankdd1a,]
meth.data.ankdd1a
#Get rid of NA value
na.count <- apply (meth.data.ankdd1a, 1, function (x) sum (as.numeric(is.na(x))))
na.count
##Exclude the NA values of cg which NA > 50%
exclude<-as.numeric(na.count>0.5*ncol(meth.data.ankdd1a))
exclude
meth.data.ankdd1a <- meth.data.ankdd1a[which (exclude==0),]
meth.data.ankdd1a [,1:4]

#Make a table of cg site methylation correlates to gene expression value and methylation value in high/low expression value 
results<-array(NA,c(nrow(meth.data.ankdd1a),4))
rownames(results)<-rownames(meth.data.ankdd1a)
colnames(results)<-c("Cor.ANKDD1A","Cor.test.ANKDD1A","Mean.high.ANKDD1A","Mean.low.ANKDD1A")
results
ANKDD1A.high2 <- as.numeric (as.numeric(rna.data2["ANKDD1A",])>median (as.numeric (rna.data2["ANKDD1A",])))
class(results)
for (i in 1:nrow (meth.data.ankdd1a))
{ 
  results [i,1] <- cor.test (as.numeric(rna.data2["ANKDD1A",]), as.numeric (meth.data.ankdd1a [i,]), use = "c")$est
  results [i,2] <- cor.test (as.numeric(rna.data2["ANKDD1A",]), as.numeric (meth.data.ankdd1a [i,]), use = "c")$p.value
}
results [,3] <- apply (meth.data.ankdd1a [,which (ANKDD1A.high2 ==1)], 1, mean, na.rm = T)
results [,4] <- apply (meth.data.ankdd1a [,which (ANKDD1A.high2 == 0)], 1, mean, na.rm = T)

#Make a plot methylation vs ANKDD1A gene expression (use the methylation data on cg09291474) 
png("COAD_ANKDD1A_Exp&Methylation.png")
plot (as.numeric(meth.data2["cg09291474",]),as.numeric(rna.data2["ANKDD1A",]), xlab = "cg09291474", ylab = "ANKDD1A RNAseq")
abline(lm(rna.data2["ANKDD1A",]~meth.data2["cg09291474",]))
dev.off()
