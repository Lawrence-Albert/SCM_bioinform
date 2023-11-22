## ---------------------------
## Author: Yukun Li
##
## Date Created: 2023-07-10
##
## Copyright (c) Yukun Li, 2023
## Email: lorenli@mail.ccmu.edu.cn

## Purpose of script:
##
## Selection of Candidate Hub Genes via Machine Learning; Modeling of a SCM Diagnostic Column Line Graph



BiocManager::install("sigFeature")



library(tidyverse)
library(glmnet)
source('msvmRFE.R')   
library(VennDiagram)
library(sigFeature)
library(e1071)
library(caret)
library(randomForest)
library(limma)


inputFile="GSE79962.txt"    
C="C"                        


rt=read.table(inputFile, header=T, sep="\t", check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
data=avereps(data)
data=t(data)
data=data[,read.table("disease.txt", header=F, sep="\t", check.names=F)[,1]]


sample=read.table("sample.txt",sep="\t",header=F,check.names=F,row.names = 1)
data=data[rownames(sample),]
afcon=sum(sample[,1]==C)
group=c(rep("0",afcon),rep("1",nrow(data)-afcon))
group=as.matrix(as.numeric(group))
rownames(group)=rownames(data)
colnames(group)="Type"
input <- as.data.frame(cbind(group,data))
input$Type=as.factor(input$Type)

svmRFE(input, k = 10, halve.above = 100) 
nfold = 10
nrows = nrow(input)
folds = rep(1:nfold, len=nrows)[sample(nrows)]
folds = lapply(1:nfold, function(x) which(folds == x))
results = lapply(folds, svmRFE.wrap, input, k=10, halve.above=100) 
top.features = WriteFeatures(results, input, save=F) 
head(top.features)

write.csv(top.features,"feature_svm.csv")




featsweep = lapply(1:100, FeatSweep.wrap, results, input) 


no.info = min(prop.table(table(input[,1])))
errors = sapply(featsweep, function(x) ifelse(is.null(x), NA, x$error))


pdf("svm-error.pdf",width = 5,height = 5)
PlotErrors(errors, no.info=no.info) 
dev.off()

pdf("svm-accuracy.pdf",width = 5,height = 5)
Plotaccuracy(1-errors,no.info=no.info) 
dev.off()


which.min(errors) 


library(survival)
library(glmnet)
library(ggplot2)
library(ggsci)
library(patchwork)
library(limma)

inputFile="GSE79962.txt"       
C="C"                        


rt=read.table(inputFile, header=T, sep="\t", check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
data=avereps(data)
data=t(data)
data=data[,read.table("disease.txt", header=F, sep="\t", check.names=F)[,1]]
sample=read.table("sample.txt",sep="\t",header=F,check.names=F,row.names = 1)
data=data[rownames(sample),]
x=as.matrix(data)


afcon=sum(sample[,1]==C)
group=c(rep("0",afcon),rep("1",nrow(data)-afcon))
group=as.matrix(group)
rownames(group)=rownames(data)
y=as.matrix(group[,1])

set.seed(123)
cvfit = cv.glmnet(x, y,family = "binomial", nlambda=100, alpha=1,nfolds = 10) 








fit <- glmnet(x,y,family = "binomial")
cvfit$lambda.min


coef <- coef(fit, s = cvfit$lambda.min)
index <- which(coef != 0)
actCoef <- coef[index]
lassoGene=row.names(coef)[index]
geneCoef=cbind(Gene=lassoGene, Coef=actCoef)
write.table(geneCoef, file="geneCoef.xls", sep="\t", quote=F, row.names=F)
write.table(file="lassoset.txt",lassoGene,sep="\t",quote=F,col.names=F,row.names=F) 


pdf("lasso.pdf",height = 5,width = 7)
layout(matrix(c(1,1,2,2), 2, 2, byrow = F))   

plot(fit,xvar = 'lambda')


plot(cvfit)
abline(v=log(c(cvfit$lambda.min,cvfit$lambda.1se)),lty="dashed")

dev.off()




library(randomForest)
library(limma)
library(ggpubr)
set.seed(123)

inputFile="GSE79962.txt"       
C="C"                        


rt=read.table(inputFile, header=T, sep="\t", check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
data=avereps(data)
data=t(data)
data=data[,read.table("disease.txt", header=F, sep="\t", check.names=F)[,1]]
sample=read.table("sample.txt",sep="\t",header=F,check.names=F,row.names = 1)
data=data[rownames(sample),]
colnames(data)=gsub("-", "afaf", colnames(data))

afcon=sum(sample[,1]==C)
group=c(rep("con",afcon),rep("treat",nrow(data)-afcon))


rf=randomForest(as.factor(group)~., data=data, ntree=500)
pdf(file="forest.pdf", width=6, height=6)
plot(rf, main="Random forest", lwd=2)
dev.off()


optionTrees=which.min(rf$err.rate[,1])
optionTrees
rf2=randomForest(as.factor(group)~., data=data, ntree=optionTrees)



importance=importance(x=rf2)
importance=as.data.frame(importance)



importance$size=rownames(importance)
importance=importance[,c(2,1)]
names(importance)=c("Gene","importance")

af=importance[order(importance$importance,decreasing = T),]
af=af[1:20,]
p=ggdotchart(af, x = "Gene", y = "importance",
             color = "importance", 
             sorting = "descending",                       
             add = "segments",                             
             add.params = list(color = "lightgray", size = 2), 
             dot.size = 6,                        
             font.label = list(color = "white", size = 9,
                               vjust = 0.5),               
             ggtheme = theme_bw()         ,               
             rotate=TRUE                                       )
p1=p+ geom_hline(yintercept = 0, linetype = 2, color = "lightgray")+
  gradient_color(palette =c(ggsci::pal_npg()(2)[2],ggsci::pal_npg()(2)[1])      ) +
  grids()   

pdf(file="importance.pdf", width=6, height=6)
print(p1)
dev.off()

rfGenes=importance[order(importance[,"importance"], decreasing = TRUE),]
write.table(rfGenes, file="rfGenes.xls", sep="\t", quote=F, col.names=T, row.names=F)







library(randomcoloR)
library(venn) 


A="LASSO"
B="RandomForest"
C="SVM-REF"


geneList=list()
rt=read.table(paste0(A,".txt"),header=F,sep="\t",check.names=F)      
geneNames=as.vector(rt[,1])                
geneNames=gsub("^ | $","",geneNames)      
uniqGene=unique(geneNames)                 
geneList[[A]]=uniqGene                   
uniqLength=length(uniqGene)
print(paste("1",uniqLength,sep=" "))
rt=read.table(paste0(B,".txt"),header=F,sep="\t",check.names=F)    
geneNames=as.vector(rt[,1])                
geneNames=gsub("^ | $","",geneNames)       
uniqGene=unique(geneNames)                
geneList[[B]]=uniqGene
uniqLength=length(uniqGene)
print(paste("3",uniqLength,sep=" "))
rt=read.table(paste0(C,".txt"),header=F,sep="\t",check.names=F)    
geneNames=as.vector(rt[,1])                
geneNames=gsub("^ | $","",geneNames)       
uniqGene=unique(geneNames)                
geneList[[C]]=uniqGene
uniqLength=length(uniqGene)
print(paste("3",uniqLength,sep=" "))

mycol <- distinctColorPalette(3)
pdf(file="hub.pdf",width=5,height=5)                                                
venn(geneList,col=mycol[1:length(geneList)],zcolor=mycol[1:length(geneList)],box=F)
dev.off()

intersectGenes=Reduce(intersect,geneList)          
write.table(file="hub.txt",intersectGenes,sep="\t",quote=F,col.names=F,row.names=F) 











install.packages("regplot")
install.packages("timeROC")
install.packages("ggDCA")


library(dplyr)
library(pROC)
library(ggplot2)
library(survival)
library(regplot)
library(rms)
library(ggsci)
library(survminer)
library(timeROC)
library(ggDCA)
library(limma)

rm (list=ls())
setwd("D:/科研/2干试�?生信/ML/GSE79962/5.诊断列线�?)

inputFile="GSE79962.txt"       
hub="hub.txt"        


rt=read.table(inputFile, header=T, sep="\t", check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
data=avereps(data)
data=t(data)
sample=read.table("sample.txt",sep="\t",header=F,check.names=F)
colnames(sample)=c("ID","Type")
data=data[sample$ID,]
aSAH1=data[,read.table(hub, header=F, sep="\t", check.names=F)[,1]]
aSAH=cbind(sample,aSAH1)

aflist=roc(Type~TIMMDC1+MRPS31+FBXO7+LYRM7+BCS1L+PGS1,
           
           data = aSAH)
g3 <- ggroc(aflist, size = 1.2,alpha=.6,)
g5=g3+ggsci::scale_color_lancet()
print(g5)


dd <- datadist(aSAH)
options(datadist="dd")
fit <- lrm(formula = Type ~ TIMMDC1+MRPS31+FBXO7+LYRM7+BCS1L+PGS1, 
           
           data =aSAH)
print(fit)
coef=as.data.frame(fit$coefficients)[-1,,drop=F]
coefout=cbind(ID=rownames(coef),coef)
write.table(coefout,file="coefficients.txt",sep="\t",quote=F,row.names = F)

pdf(file="nomogram.pdf", width=9, height=7.5)
plot(nomogram(fit,fun.at = seq(0.05,0.95,0.05)),funlabel = "nomogram model")
dev.off()







