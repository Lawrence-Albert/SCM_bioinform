## ---------------------------
## Author: Yukun Li
##
## Date Created: 2023-07-10
##
## Copyright (c) Yukun Li, 2023
## Email: lorenli@mail.ccmu.edu.cn

## Purpose of script:
##
## Identification and Functional Analysis of MitoDEGs in SCM.


library(randomcoloR)
library(venn) 


A="WGCNA"
B="DIFF"


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

mycol <- distinctColorPalette(100)

pdf(file="disease.pdf",width=5,height=5)                                                
venn(geneList,col=mycol[1:length(geneList)],zcolor=mycol[1:length(geneList)],box=F)
dev.off()

intersectGenes=Reduce(intersect,geneList)          
write.table(file="disease.txt",intersectGenes,sep="\t",quote=F,col.names=F,row.names=F) 






A="DEG"
B="Mito"


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

mycol <- distinctColorPalette(100)

pdf(file="venn.pdf",width=5,height=5)                                                
venn(geneList,col=mycol[1:length(geneList)],zcolor=mycol[1:length(geneList)],box=F)
dev.off()

intersectGenes=Reduce(intersect,geneList)          
write.table(file="venn.txt",intersectGenes,sep="\t",quote=F,col.names=F,row.names=F) 










library("clusterProfiler")
library("org.Hs.eg.db")
library("enrichplot")
library("ggplot2")
library("pathview")
library("DOSE")

pvalueFilter=0.05        
qvalueFilter=1       
showNum=20

rt=read.table("disease.txt",sep="\t",check.names=F,header=F)      
genes=as.vector(rt[,1])
entrezIDs <- mget(genes, org.Hs.egSYMBOL2EG, ifnotfound=NA)  
entrezIDs <- as.character(entrezIDs)
rt=cbind(rt,entrezID=entrezIDs)
colnames(rt)=c("symbol","entrezID") 
rt=rt[is.na(rt[,"entrezID"])==F,]                        
gene=rt$entrezID
gene=unique(gene)
colorSel="qvalue"
if(qvalueFilter>0.05){
	colorSel="pvalue"
}
kk <- enrichDO(gene =gene, ont = "DO",pvalueCutoff = 1, qvalueCutoff = 1)
KEGG=as.data.frame(kk)
KEGG$geneID=as.character(sapply(KEGG$geneID,function(x)paste(rt$symbol[match(strsplit(x,"/")[[1]],as.character(rt$entrezID))],collapse="/")))
KEGG=KEGG[(KEGG$pvalue<pvalueFilter & KEGG$qvalue<qvalueFilter),]
write.table(KEGG,file="DO.xls",sep="\t",quote=F,row.names = F)
if(nrow(KEGG)<showNum){
	showNum=nrow(KEGG)
}
pdf(file="DO_barplot.pdf",width =9,height = 7)
barplot(kk, drop = TRUE, showCategory = showNum, color = colorSel)+scale_y_discrete(labels=function(x) stringr::str_wrap(x, width=60))
dev.off()
pdf(file="DO_bubble.pdf",width =9,height = 7)
dotplot(kk, showCategory = showNum, orderBy = "GeneRatio",color = colorSel)+scale_y_discrete(labels=function(x) stringr::str_wrap(x, width=60))
dev.off()











library("org.Hs.eg.db")  
library("clusterProfiler")
library("enrichplot")
library("ggplot2")
library("ggnewscale")
library("enrichplot")
library("DOSE")
library(stringr)

pvalueFilter=0.05         
qvalueFilter=1  
showNum=7

rt=read.table("disease.txt",sep="\t",check.names=F,header=F)      
genes=as.vector(rt[,1])
entrezIDs <- mget(genes, org.Hs.egSYMBOL2EG, ifnotfound=NA)  
entrezIDs <- as.character(entrezIDs)
rt=cbind(rt,entrezID=entrezIDs)
colnames(rt)=c("symbol","entrezID") 
rt=rt[is.na(rt[,"entrezID"])==F,]                        
gene=rt$entrezID
gene=unique(gene)

colorSel="qvalue"
if(qvalueFilter>0.05){
  colorSel="pvalue"
}


kk=enrichGO(gene = gene,OrgDb = org.Hs.eg.db, pvalueCutoff =1, qvalueCutoff = 1, ont="all", readable =T)
GO=as.data.frame(kk)
GO=GO[(GO$pvalue<pvalueFilter & GO$qvalue<qvalueFilter),]

write.table(GO,file="GO.xls",sep="\t",quote=F,row.names = F)


if(nrow(GO)<30){
  showNum=nrow(GO)
}

pdf(file="GO_barplot.pdf",width = 9,height =7)
bar=barplot(kk, drop = TRUE, showCategory =showNum,split="ONTOLOGY",color = colorSel) + facet_grid(ONTOLOGY~., scale='free')+scale_y_discrete(labels=function(x) stringr::str_wrap(x, width=60))
print(bar)
dev.off()


pdf(file="GO_bubble.pdf",width = 9,height =7)
bub=dotplot(kk,showCategory = showNum, orderBy = "GeneRatio",split="ONTOLOGY", color = colorSel) + facet_grid(ONTOLOGY~., scale='free')+scale_y_discrete(labels=function(x) stringr::str_wrap(x, width=60))
print(bub)
dev.off()







R.utils::setOption( "clusterProfiler.download.method",'auto' )


library("clusterProfiler")
library("org.Hs.eg.db")
library("enrichplot")
library("ggplot2")
library("pathview")
library("ggnewscale")
library("DOSE")
library(stringr)

pvalueFilter=0.05        
qvalueFilter=1        
showNum=20

rt=read.table("2.DIFF_all.txt",sep="\t",check.names=F,header=T)  
rownames(rt)=rt[,1]
rt=rt[read.table("disease.txt",sep="\t",check.names=F,header=F)[,1] ,c(1,2)]

genes=as.vector(rt[,1])
entrezIDs <- mget(genes, org.Hs.egSYMBOL2EG, ifnotfound=NA)  
entrezIDs <- as.character(entrezIDs)
rt=cbind(rt,entrezID=entrezIDs)
colnames(rt)=c("symbol","logFC","entrezID") 

rt[rt[,"entrezID"]=="NA",]
rt=rt[rt[,"entrezID"]!="NA",]                        
gene=rt$entrezID
gene=unique(gene)
aflogfc=rt$logFC
names(aflogfc)=rt$symbol
colorSel="qvalue"
if(qvalueFilter>0.05){
	colorSel="pvalue"
}
kk <- enrichKEGG(gene = gene, organism = "hsa", pvalueCutoff =1, qvalueCutoff =1)
KEGG=as.data.frame(kk)
KEGG$geneID=as.character(sapply(KEGG$geneID,function(x)paste(rt$symbol[match(strsplit(x,"/")[[1]],as.character(rt$entrezID))],collapse="/")))
KEGG=KEGG[(KEGG$pvalue<pvalueFilter & KEGG$qvalue<qvalueFilter),]
write.table(KEGG,file="KEGG.xls",sep="\t",quote=F,row.names = F)

if(nrow(KEGG)<showNum){
	showNum=nrow(KEGG)
}

pdf(file="KEGG_barplot.pdf",width = 9,height = 7)
barplot(kk, drop = TRUE, showCategory = showNum, color = colorSel) +scale_y_discrete(labels=function(x) stringr::str_wrap(x, width=60))
dev.off()

pdf(file="KEGG_bubble.pdf",width = 9,height = 7)
dotplot(kk, showCategory = showNum, orderBy = "GeneRatio",color = colorSel)+scale_y_discrete(labels=function(x) stringr::str_wrap(x, width=60))
dev.off()

pdf(file="KEGG_cnet.pdf",width = 10,height = 8)
af=setReadable(kk, 'org.Hs.eg.db', 'ENTREZID')
cnetplot(af, showCategory = 5, categorySize="pvalue",circular = TRUE,colorEdge = TRUE,cex_label_category=0.65,cex_label_gene=0.6,foldChange = aflogfc)
dev.off()

pdf(file="KEGG_net.pdf",width = 9,height = 7)
x2 <- pairwise_termsim(kk)
emapplot(x2,showCategory = showNum,cex_label_category=0.65,color = "pvalue",layout ="nicely")
dev.off()  

pdf(file="KEGG_heatplot.pdf",width = 10,height = 7)
kegg=setReadable(kegg, 'org.Hs.eg.db', 'ENTREZID')
heatplot(kegg,foldChange = aflogfc) 
dev.off()

keggId="hsa05171"
geneFC=rt$logFC
names(geneFC)=gene
pv.out=pathview(gene.data = geneFC, pathway.id = keggId, species = "hsa", out.suffix = "pathview")
p <- pathview(gene.data = geneFC, pathway.id = keggId, species = "hsa", kegg.native = F, sign.pos="bottomleft", same.layer = F)






