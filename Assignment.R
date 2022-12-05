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
coad.os <- Surv(os.time,os.event)


#Create a univariate model table
results.univariate <- array (NA, c(nrow(rna.data),4))
results.univariate [1:3,]
colnames (results.univariate) <- c("HR", "LCI", "UCI", "PVAL")
rownames (results.univariate) <- rownames (rna.data)
results.univariate [1:3,]
results.univariate <- as.data.frame (results.univariate)

for (i in 1:nrow(rna.data)){
  coxphmodel <- coxph(coad.os ~ as.numeric(rna.data[i,]))
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

#Determine the mean and median expression value in rna.data
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

#Make a table to collect HR, Expval$Median, and Expval$Mean, order them in the descending order according to univariate regression
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
ExpvalMedian

#Fit a multivariate model.
clin.data <- clin.data [colnames (rna.data),]
clin.data [1:10,]
colnames (clin.data)
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
#Factor5: Gender
gender <- clin.data$gender
gender [which (gender == "MALE")] <- 1
gender [which (gender == "FEMALE")] <- 0
gender
#Factor6: Pathologic stage
clin.data$pathologic_stage
x3 <- grep ("III", clin.data$pathologic_stage)
x4 <- grep ("IV", clin.data$pathologic_stage)
stage.high <- rep (0, nrow (clin.data))
stage.high [c(x3,x4)] <- 1
stage.high
#whether the factors associate with survival data??
summary (coxph(coad.os~stage.high))$coef
summary (coxph(coad.os~age))$coef
summary (coxph (coad.os~lymph.node))$coef
summary (coxph (coad.os~bmi))$coef
summary (coxph (coad.os~gender))$coef
summary (coxph (coad.os~polyps))

#Make a multivariate table
results.multivariate <- array (NA, c(nrow (rna.data),4))
colnames (results.multivariate) <- c ("HR", "LCI", "UCI", "PVAL")
rownames (results.multivariate) <- rownames (rna.data)
results.multivariate <- as.data.frame (results.multivariate)
for(i in 1:nrow(rna.data)){
  coxphmodel <- coxph(coad.os ~ rna.data[i,]+age+stage.high+lymph.node+gender)
  results.multivariate$HR[i] <- summary(coxphmodel)$coef[1,2]
  results.multivariate$LCI[i] <- summary(coxphmodel)$conf.int[1,3]
  results.multivariate$UCI[i] <- summary(coxphmodel)$conf.int[1,4]
  results.multivariate$PVAL[i] <- summary(coxphmodel)$coef[1,5]
}
results.multivariate<-as.data.frame(results.multivariate)
results.multivariate$FDR<-p.adjust(results.multivariate$PVAL,method="fdr")
results.multivariate<-results.multivariate[order(results.multivariate$FDR, decreasing=F),]
results.multivariate["ANKDD1A",]
head (results.multivariate)
#Select Genes which FDR <= 0.05
PotentialGenesMulti <- results.multivariate [which (results.multivariate$FDR <= 0.05),]
nrow (PotentialGenesMulti)
PotentialGenesMulti
#Make a table for multivariate to collect HR, Expval$Median, and Expval$Mean, order them in the descending order
#find shared genes in both regression analysis
overlap.genes <- PotentialGenesMulti [which (is.element (rownames (PotentialGenesMulti), rownames (PotentialGenes))),]
overlap.genes

ExpvalMedianMulti <- array (NA, c(nrow(overlap.genes),3))
rownames (ExpvalMedianMulti) <- rownames (overlap.genes)
colnames (ExpvalMedianMulti) <- c("HR","Median","Mean")
ExpvalMedianMulti <- as.data.frame (ExpvalMedianMulti)
for (i in 1:nrow (overlap.genes)){
  ExpvalMedianMulti$HR[i] <- overlap.genes$HR[i]
  ExpvalMedianMulti$Median[i] <- getElement(summary(rna.data[rownames(overlap.genes)[i],]), "Median")
  ExpvalMedianMulti$Mean[i] <- getElement(summary (rna.data[rownames(overlap.genes)[i],]), "Mean")
}

#Select high expression genes and ordering by HR
ExpvalMedianMulti <- ExpvalMedianMulti[order(ExpvalMedianMulti$Median, decreasing = T),]
ExpvalMedianMulti.high <- ExpvalMedianMulti[1:10,]
ExpvalMedianMulti.high <- ExpvalMedianMulti.high[order(ExpvalMedianMulti.high$HR, decreasing = T),]
ExpvalMedianMulti.high
write.csv (ExpvalMedianMulti.high, file = "ExpvalMedian.csv", quote=FALSE)

#Finally, according to multivariate results. ANKDD1A was selected for the biomarker study.
results.univariate["ANKDD1A",]
results.multivariate["ANKDD1A",]
summary (rna.data["ANKDD1A",])

#Stratify patients ANKDD1A expression and draw KM plot
png ("COAD_OS_byANKDD1A.png", width = 6, height = 6, units = 'in', res = 300)
ANKDD1A.high <- as.numeric (rna.data ["ANKDD1A",] > median (rna.data["ANKDD1A",]))
plot (survfit(coad.os~ANKDD1A.high), col = c("black","red"), lwd =2, mark.time = TRUE, xlab = "OS time (days)", ylab = "proportion surviving", main = "ANKDD1A in Colon Cancer")
legend("topright",legend=c("ANKDD1A-high","ANKDD1A-low"),col=c("red","black"),lwd=2)
text(3000,0.1,"HR=1.93 (95%CI 1.40-2.64)")
dev.off()

#DNA Methylation Detection
meth.annot <- readRDS("annot450k.rds")
dim(meth.annot)
##Loading DNA methylation data download from TCGA
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
ncol(meth.data.ankdd1a)
#Make a table of cg site methylation correlates to gene expression value and methylation value in high/low expression value 
results<-array(NA,c(nrow(meth.data.ankdd1a),4))
rownames(results)<-rownames(meth.data.ankdd1a)
colnames(results)<-c("Cor.ANKDD1A","Cor.test.ANKDD1A","Mean.high.ANKDD1A","Mean.low.ANKDD1A")
ANKDD1A.high2 <- as.numeric (as.numeric(rna.data2["ANKDD1A",])>median (as.numeric (rna.data2["ANKDD1A",])))
ANKDD1A.high2
for (i in 1:nrow (meth.data.ankdd1a))
{ 
  results [i,1] <- cor.test (as.numeric(rna.data2["ANKDD1A",]), as.numeric (meth.data.ankdd1a [i,]), use = "c")$est
  results [i,2] <- cor.test (as.numeric(rna.data2["ANKDD1A",]), as.numeric (meth.data.ankdd1a [i,]), use = "c")$p.value
}
results [,3] <- apply (meth.data.ankdd1a [,which (ANKDD1A.high2 ==1)], 1, mean, na.rm = T)
results [,4] <- apply (meth.data.ankdd1a [,which (ANKDD1A.high2 == 0)], 1, mean, na.rm = T)
results
write.csv (results, file = "Methylation.csv", quote=FALSE)

#Find out normal and cancer cells
colnames (meth.data.ankdd1a)
colnames145 <- substr (colnames (meth.data.ankdd1a),14,15)
colnames145
colnames145 <- as.numeric(colnames145)
patientCancer <- ifelse (colnames145 <10, "tumour", "normal")
patientCancer <- factor (patientCancer)
group_list = patientCancer
group_list
##Finally in the clin.data. We got 19 normal cells and 280 cancer cells
table (group_list)
#Gene Expression difference 
normal <- which (group_list == "normal")
cancer <- which (group_list == "tumour")

#Combination of ANKDD1A expression data and methylation data
expMeth <- array (NA, dim = c(ncol(rna.data2), 3))
rownames(expMeth) <- colnames (rna.data2)
colnames (expMeth) <- c("ExpVal", "Betavalue", "Type")
expMeth
for (i in 1:ncol(rna.data2)){
  expMeth [i,1] <- as.numeric(rna.data2 ["ANKDD1A",i])
  expMeth [i,2] <- as.numeric (meth.data.ankdd1a ["cg09291474", i])
}
expMeth [normal, 3] <- "Normal"
expMeth [cancer, 3] <- "Cancer"
expMeth <- as.data.frame (expMeth)
expMeth$ExpVal <- as.numeric (expMeth$ExpVal)
expMeth$Betavalue <- as.numeric (expMeth$Betavalue)
expMeth$Type <- factor (expMeth$Type)

#Survival status add into expMeth
death <- rownames(surv.data2[which (surv.data2$OS == 1),])
alive <- rownames(surv.data2[which (surv.data2$OS == 0),])
expMeth [death,4] <- "death"
expMeth [alive,4] <- "alive"
colnames (expMeth)[4] <- "Vital"
#Delete normal sample but death
expMeth <- expMeth [c(-5,-7,-251),]
expMeth <- expMeth [c(-7,-13,-19,-109,-111,-113,-115),]
expMeth$Vital <- factor (expMeth$Vital)
expMeth$Type <- factor (expMeth$Type)

#Cancer patients ANKDD1A Expression Alive vs Death
library (ggplot2)
expMeth.cancer <- expMeth [which (expMeth$Type == "Cancer"),]
ExpVal.death <- expMeth.cancer [which (expMeth.cancer$Vital == "death"), 1]
ExpVal.alive <- expMeth.cancer [which (expMeth.cancer$Vital == "alive"), 1]
Meth.death <- expMeth.cancer [which (expMeth.cancer$Vital == "death"), 2]
Meth.alive <- expMeth.cancer [which (expMeth.cancer$Vital == "alive"), 2]
Survival <- expMeth.cancer$Vital
Survival
ggplot(data=expMeth.cancer, aes(x=Betavalue, y=ExpVal, color = Survival)) + xlab ("cg09291474 DNA Methylation") + ylab ("ANKDD1A Expression Value") + geom_point()+ geom_smooth(method = "lm")
shapiro.test(ExpVal.death)
shapiro.test(Meth.death)
cor.test (ExpVal.death, Meth.death)
shapiro.test(ExpVal.alive)
shapiro.test(Meth.alive)
cor.test (ExpVal.alive, Meth.alive)
t.test (ExpVal.alive,ExpVal.death)

##Use CNV Gistic2 
Gistic2.data <- read.table ("TCGA.COAD.sampleMap_Gistic2_CopyNumber_Gistic2_all_thresholded.by_genes", sep = "\t", head = TRUE, row.names = 1)
head (Gistic2.data)
dim (Gistic2.data)
Gistic2.data <- Gistic2.data [,which(is.element(colnames (Gistic2.data),colnames(rna.data)))]
Gistic2.data ["ANKDD1A",]
Gistic2.data <- as.matrix(Gistic2.data)
amplification.counts<- apply(Gistic2.data,MARGIN=1,function(x)sum(x==1 | x == 2))
deletion.counts <- apply (Gistic2.data,MARGIN=1, function(x)sum(x==-1 | x== -2))
normal.counts <- apply (Gistic2.data,MARGIN=1,function(x)sum(x==0))
CNA.frame <- data.frame(ID=rownames(Gistic2.data),Normal = normal.counts,AmpCount=amplification.counts, DelCount=deletion.counts)
CNA.frame ["ANKDD1A",]
write.csv (CNA.frame ["ANKDD1A",], file = "cnv.csv", quote=FALSE)

#Make an array to collect expression value and cna
ExpCNV <- array (NA, c(ncol(Gistic2.data), 2))
head (ExpCNV)                 
rownames (ExpCNV) <- colnames (Gistic2.data)
colnames (ExpCNV) <- c("CNA.Type", "ExpVal")
nrow (ExpCNV)
ncol (Gistic2.data)
head (rna.data2)
rownames (ExpCNV)[1]
for (i in 1:ncol (Gistic2.data)){
  ExpCNV[i,1] <- Gistic2.data["ANKDD1A",rownames (ExpCNV)[i]]
  ExpCNV[i,2] <- rna.data ["ANKDD1A", rownames (ExpCNV)[i]]
}
ExpCNV <- as.data.frame(ExpCNV)
ExpCNV$CNA.Type <- as.numeric (ExpCNV$CNA.Type)
ExpCNV$CNA.Type [which (ExpCNV$CNA.Type == 1 )] <- "Amplification"
ExpCNV$CNA.Type [which (ExpCNV$CNA.Type == 2 )] <- "Amplification"
ExpCNV$CNA.Type [which (ExpCNV$CNA.Type == 0)] <- "Normal"
ExpCNV$CNA.Type [which (ExpCNV$CNA.Type == -1)] <- "Deletion"
ExpCNV$CNA.Type [which (ExpCNV$CNA.Type == -2)] <- "Deletion"
ExpCNV$CNA.Type <- factor (ExpCNV$CNA.Type, levels = c("Deletion", "Normal", "Amplification"))
ExpCNV$CNA.Type 

png("ANKDD1A Expression vs CNA.png")
ggplot(data=ExpCNV, aes(x=CNA.Type, y=ExpVal, fill = CNA.Type)) + geom_boxplot()+ labs (x = "Type of CNA", y = "ANKDD1A Expression Value")
dev.off()
shapiro.test(ExpCNV[which (ExpCNV$CNA.Type == "Normal"),]$ExpVal)
shapiro.test(ExpCNV[which (ExpCNV$CNA.Type == "Amplification"),]$ExpVal)
shapiro.test(ExpCNV[which (ExpCNV$CNA.Type == "Deletion"),]$ExpVal)
t.test (ExpCNV[which (ExpCNV$CNA.Type == "Deletion"),]$ExpVal, ExpCNV[which (ExpCNV$CNA.Type == "Normal"),]$ExpVal)
t.test (ExpCNV[which (ExpCNV$CNA.Type == "Deletion"),]$ExpVal, ExpCNV[which (ExpCNV$CNA.Type == "Amplification"),]$ExpVal)
t.test (ExpCNV[which (ExpCNV$CNA.Type == "Amplification"),]$ExpVal, ExpCNV[which (ExpCNV$CNA.Type == "Normal"),]$ExpVal)

#Somatic mutation
soma.data <- read.table ("mc3_COAD_mc3.txt", sep = "\t", head = TRUE)
soma.data$sample<-gsub(soma.data$sample, pattern="-", replace=".")
dim (soma.data)
soma.data2 <- soma.data [which (is.element(soma.data$sample, colnames(rna.data))),]
soma.data [1:3,1:5]
sNames.soma <- unique (as.character(soma.data$sample))
sNames.soma
soma.data.ANKDD1A <- soma.data[which (soma.data$gene == "ANKDD1A"), ]
soma.data.ANKDD1A
#Four results
rna.data["ANKDD1A","TCGA.CA.6716.01"]
rna.data["ANKDD1A","TCGA.CA.6717.01"]
rna.data["ANKDD1A","TCGA.CM.6675.01"]
rna.data["ANKDD1A","TCGA.F4.6856.01"]

#rppa
rppa.data <- read.table ("TCGA.COAD.sampleMap_RPPA_RBN", sep = "\t", head = TRUE, row.name = 1)
dim (rppa)
rppa.data <- rppa.data [,which (is.element(colnames (rppa), colnames (rna.data)))]
rppa.data["ANKDD1A",]
##No ANKDD1A results
