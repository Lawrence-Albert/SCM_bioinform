## ---------------------------
## Author: Yukun Li
##
## Date Created: 2023-07-10
##
## Copyright (c) Yukun Li, 2023
## Email: lorenli@mail.ccmu.edu.cn

## Purpose of script:
##
## Infiltration of immune cell types compared between SCM and CON; Correlation and pan-tissue analysis between six hub genes and inflammation disorder.


library(limma)
library(reshape2)
library(tidyverse)
library(ggplot2)


immfile="CIBERSORT.txt"         
exp="1.rawexp_GSE54236.txt"       
geness= "keycluster.txt"           
highcol="
midcol="white"                   
lowcol= "

immune=read.table(immfile, header=T, sep="\t", check.names=F, row.names=1)
colnames(immune)=gsub("_CIBERSORT","",colnames(immune))
immune=immune[immune[,"P-value"]<0.05,]
data=as.matrix(immune[,1:(ncol(immune)-3)])
expdata=read.table(exp, header=T, sep="\t", check.names=F, row.names=1)
af=read.table(file =geness,sep = "\t",header = F,check.names = F)
expdata=expdata[intersect(af[,1],rownames(expdata)),]
sameSample=intersect(rownames(data),colnames(expdata))
data=data[sameSample,]
expdata=expdata[,sameSample]
expdata=t(expdata)
data1=data

outTab=data.frame()
for(immune in colnames(data)){
	for(gene in colnames(expdata)){
		x=as.numeric(data[,immune])
		y=as.numeric(expdata[,gene])
		corT=cor.test(x,y,method="spearman")
		cor=corT$estimate
		pvalue=corT$p.value
		text=ifelse(pvalue<0.001,"***",ifelse(pvalue<0.01,"**",ifelse(pvalue<0.05,"*","")))
		outTab=rbind(outTab,cbind(Gene=gene, Immune=immune, cor, text, pvalue))
	}
}
write.table(outTab,quote = F,sep = "\t","cor.xls",row.names = F)
outTab$cor=as.numeric(outTab$cor)

pdf(file="cor.pdf", width=10, height=10)
ggplot(outTab, aes(Gene, Immune)) + 
  geom_tile(aes(fill = cor), colour = "grey", size = 1)+
  scale_fill_gradient2(high=highcol, mid = midcol,low=lowcol) + 
  geom_text(aes(label=text),col ="black",size = 3) +
  theme_minimal() +   
  theme(axis.title.x=element_blank(), axis.ticks.x=element_blank(), axis.title.y=element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 10, face = "bold"),  
        axis.text.y = element_text(size = 10, face = "bold")) +
  labs(fill =paste0("***  p<0.001","\n", "**  p<0.01","\n", " *  p<0.05","\n", "\n","Correlation")) +   
  scale_x_discrete(position = "bottom")+coord_flip()
dev.off()







library(ggsci)
library(randomcoloR)
library(limma)
library(reshape2)
library(ggplot2)
library(ggpubr)
immFile="ssgsea.txt"   
cluFile="sample.txt"             
immune=read.table(immFile, header=T, sep="\t", check.names=F, row.names=1)
immune=immune[immune[,"P-value_CIBERSORT"]<1,]
data=as.matrix(immune[,1:(ncol(immune)-3)])
colnames(immune)=gsub("_CIBERSORT"," ",colnames(immune))
colnames(data)=gsub("_CIBERSORT"," ",colnames(data))
colnames(immune)=gsub("_"," ",colnames(immune))
colnames(data)=gsub("_"," ",colnames(data))
cluster=read.table(cluFile, header=F, sep="\t", check.names=F, row.names=1)
colnames(cluster)="Type"
sameSample=intersect(row.names(data), row.names(cluster))
data=cbind(data[sameSample,,drop=F], cluster[sameSample,,drop=F])
data=data[order(data$Type),]
gaps=c(1, as.vector(cumsum(table(data$Type))))
xlabels=levels(factor(data$Type))

data=melt(data,id.vars=c("Type"))
colnames(data)=c("Type", "Immune", "Expression")

group=levels(factor(data$Type))
data$Type=factor(data$Type, levels=group)
bioCol=pal_jco()(6)
bioCol=bioCol[1:length(group)]
boxplot=ggboxplot(data, x="Immune", y="Expression",  fill="Type",
				  xlab="",
				  ylab="NES",
				  legend.title="Typer", 
				  width=0.8,
				  font.label = list(size = 8, color = "black"),
				  palette=bioCol,add.params = list(size=0.1))
boxplot=boxplot+stat_compare_means(aes(group=Type),symnum.args=list(cutpoints=c(0, 0.001, 0.01, 0.05, 1), symbols=c("***", "**", "*", "ns")), label="p.signif")+theme_bw()+rotate_x_text(50)

pdf(file="immune.diff.pdf", width=7, height=4.2)
print(boxplot)







options(stringsAsFactors=F)
library(corrplot)
library(circlize)
library(limma)
library(PerformanceAnalytics)


rt=read.table("SSGSEA1.txt",sep="\t",header=T,check.names=F,row.names=1)
rt=t(rt)  
M=cor(rt)     
res <- cor.mtest(rt)





























pdf(file="corpot4.pdf",width=7,height=7)
corrplot(M, order = "original", 
         type = "upper", 
         tl.pos = "lt",tl.col="black",	tl.cex = 0.5)  
corrplot(M, add = TRUE, type = "lower", method = "number", order = "AOE",
         diag = FALSE, tl.pos = "n", cl.pos = "n",number.cex=0.4 )     
dev.off()


pdf(file="corpot5.pdf",width=8,height=8)
chart.Correlation(M,method = "spearman")
dev.off()


pdf(file="corpot6.pdf",width=8,height=8)
corrplot(M, method="ellipse",p.mat = res$p, sig.level = 0.2,order = "AOE", type = "upper", tl.pos = "d",tl.cex = 0.45)
corrplot(M, add = TRUE, p.mat = res$p, sig.level = 0.2,type = "lower", method = "number", order = "AOE",
         diag = FALSE, tl.pos = "n", cl.pos = "n",number.cex=0.5)
dev.off()


pdf(file="corpot.pdf", width=8, height=8)
corrplot(M,
         order="original",
         method = "circle",
         type = "upper",
         tl.cex=0.8, pch=T,
         p.mat = res$p,
         insig = "label_sig",
         pch.cex = 1.6,
         sig.level=0.05,
         number.cex = 1,
         col=colorRampPalette(c("blue", "white", "red"))(50),
         tl.col="black")
dev.off()

