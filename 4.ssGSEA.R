## ---------------------------
## Author: Yukun Li
##
## Date Created: 2023-07-10
##
## Copyright (c) Yukun Li, 2023
## Email: lorenli@mail.ccmu.edu.cn

## Purpose of script:
##
## Gene set enrichment of the hub genes






library(ggplot2)
library(limma)
library(pheatmap)
library(ggsci)
lapply(c('clusterProfiler','enrichplot','patchwork'), function(x) {library(x, character.only = T)})
library(org.Hs.eg.db)
library(patchwork)


expFile="GSE79962.txt"          
hub="hub.txt"    

rt=read.table(expFile,sep="\t",header=T,check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
data=avereps(data)

geneaf=read.table(hub,sep="\t",header=F,check.names=F)[,1]


for (genei in geneaf) {
  


group <- ifelse(data[c(genei),]> median(data[c(genei),]), "High", "Low")   
group <- factor(group,levels = c("High","Low"))


design <- model.matrix(~0+group)
colnames(design) <- levels(group)
fit <- lmFit(data,design)
cont.matrix<-makeContrasts(High-Low,levels=design)
fit2 <- contrasts.fit(fit, cont.matrix)
fit2 <- eBayes(fit2)
deg=topTable(fit2,adjust='fdr',number=nrow(data))
Diff=deg






Diff=Diff[order(as.numeric(as.vector(Diff$logFC))),]
diffGene=as.vector(rownames(Diff))
diffLength=length(diffGene)
afGene=c()
if(diffLength>(60)){
  afGene=diffGene[c(1:30,(diffLength-30+1):diffLength)]
}else{
  afGene=diffGene
}
afExp=data[afGene,]

Type1=as.data.frame(group)
Type1=Type1[order(Type1$group,decreasing = T),,drop=F]
Type=Type1[,1]
names(Type)=rownames(Type1)
Type=as.data.frame(Type)

anncolor=list(Type=c(High="red",Low="blue"  ))



logFC_t=0
deg$g=ifelse(deg$P.Value>0.05,'stable',
             ifelse( deg$logFC > logFC_t,'UP',
                     ifelse( deg$logFC < -logFC_t,'DOWN','stable') )
)
table(deg$g)

deg$symbol=rownames(deg)
df <- bitr(unique(deg$symbol), fromType = "SYMBOL",
           toType = c( "ENTREZID"),
           OrgDb = org.Hs.eg.db)
DEG=deg
DEG=merge(DEG,df,by.y='SYMBOL',by.x='symbol')
data_all_sort <- DEG %>% 
  arrange(desc(logFC))

geneList = data_all_sort$logFC 
names(geneList) <- data_all_sort$symbol 
head(geneList)


hallmarks <- read.gmt("h.all.v7.5.1.symbols.gmt")

kk2 <- GSEA(geneList,TERM2GENE =hallmarks)







class(kk2)
colnames(kk2@result)
kegg_result <- as.data.frame(kk2)
rownames(kk2@result)[head(order(kk2@result$enrichmentScore))]
af=as.data.frame(kk2@result)
write.table(af,file=paste0("2.",paste0(genei,"_all_GSEA.xls")),sep="\t",quote=F,col.names=T)


num=5
pdf(paste0("2.",paste0(genei,"_down_GSEA.pdf")),width = 8,height = 8)
af=gseaplot2(kk2, geneSetID = rownames(kk2@result)[head(order(kk2@result$enrichmentScore),num)])
print(af)
dev.off()
pdf(paste0("2.",paste0(genei,"_up_GSEA.pdf")),width = 8,height = 8)
af=gseaplot2(kk2, geneSetID = rownames(kk2@result)[tail(order(kk2@result$enrichmentScore),num)])
print(af)
dev.off()

num=5
pdf(paste0("2.",paste0(genei,"_all_GSEA.pdf")),width = 4,height = 4.5)
af=gseaplot2(kk2, 
             base_size = 6,
             rel_heights = c(1.5, 0.05, 0.05),
             geneSetID = rownames(kk2@result)[c(head(order(kk2@result$enrichmentScore),num),tail(order(kk2@result$enrichmentScore),num))])
print(af)
dev.off()


}






