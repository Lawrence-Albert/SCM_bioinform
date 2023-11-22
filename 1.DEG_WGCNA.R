## ---------------------------
## Author: Yukun Li
##
## Date Created: 2023-07-10
##
## Copyright (c) Yukun Li, 2023
## Email: lorenli@mail.ccmu.edu.cn

## Purpose of script:
##
## Data preprocessing for DEG ;Enrichment levels in genomic weighted gene co-expression network analysis (WGCNA).


if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("impute")

library(ggplot2)
library(limma)
library(pheatmap)
library(ggsci)
library(dplyr)
lapply(c('clusterProfiler','enrichplot','patchwork'), function(x) {library(x, character.only = T)})
library(org.Hs.eg.db)
library(patchwork)
library(WGCNA)
library(GSEABase)

GSE="GSE79962"    
C="C"              
P="P"              
Ccol="blue"        
Pcol="red"         
rt=read.table(paste0(GSE,".txt"),sep="\t",header=T,check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
rt=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
rt=avereps(rt)

sample=read.table("sample.txt",
                  sep="\t",header=F,check.names=F,row.names = 1)
rt=rt[,rownames(sample)]
afcon=sum(sample[,1]==C)
max(rt)
if(max(rt)>50) rt=log2(rt+1)     
rt1=normalizeBetweenArrays(as.matrix(rt))
cols=rainbow(ncol(rt)) 
pdf(file = "1.raw.pdf",width=5,height = 4)
par(cex = 0.7,mar=c(8,8,8,8))
if(ncol(rt)>40) par(cex = 0.5,mar=c(8,8,8,8))  
boxplot(rt,las=2,col =cols ) 
dev.off()

cols=rainbow(ncol(rt1)) 
pdf(file = "1.nor.pdf",width=5,height = 4.5)
par(cex = 0.5,mar=c(8,8,8,8))
if(ncol(rt1)>40) par(cex = 0.5,mar=c(8,8,8,8)) 
boxplot(rt1,las=2,col =cols ) 
dev.off()


rt2=rbind(ID=colnames(rt1),rt1)
write.table(rt2,file=paste0("1.","norexp_",GSE,".txt"),sep="\t",quote=F,col.names = F)


rt3=rbind(ID=colnames(rt),rt)
write.table(rt3,file=paste0("1.","rawexp_",GSE,".txt"),sep="\t",quote=F,col.names = F)


data=rt1


conData=data[,as.vector(colnames(data)[1:afcon])]
aftreat=afcon+1
treatData=data[,as.vector(colnames(data)[aftreat:ncol(data)])]
rt=cbind(conData,treatData)
conNum=ncol(conData)
treatNum=ncol(treatData)


Type=c(rep("con",conNum),rep("treat",treatNum))
design <- model.matrix(~0+factor(Type))
colnames(design) <- c("con","treat")
fit <- lmFit(rt,design)
cont.matrix<-makeContrasts(treat-con,levels=design)
fit2 <- contrasts.fit(fit, cont.matrix)
fit2 <- eBayes(fit2)

Diff=topTable(fit2,adjust='fdr',number=length(rownames(data)))

DIFFOUT=rbind(id=colnames(Diff),Diff)
write.table(DIFFOUT,file=paste0("2.","DIFF_all.xls"),sep="\t",quote=F,col.names=F)



diffSig=Diff[with(Diff, (abs(logFC)>0.1 & adj.P.Val < 0.05 )), ]
diffSigOut=rbind(id=colnames(diffSig),diffSig)
write.table(diffSigOut,file=paste0("2.","DIFF_less.xls"),sep="\t",quote=F,col.names=F)


Diff=Diff[order(as.numeric(as.vector(Diff$logFC))),]
diffGene=as.vector(rownames(Diff))
diffLength=length(diffGene)
afGene=c()
if(diffLength>(60)){                                   
  afGene=diffGene[c(1:30,(diffLength-30+1):diffLength)]
}else{
  afGene=diffGene
}
afExp=rt[afGene,]

Type=c(rep(C,conNum),rep(P,treatNum))
names(Type)=colnames(rt)
Type=as.data.frame(Type)

anncolor=list(Type=c(C=Ccol,P=Pcol))
names(anncolor[[1]])=c(C,P)

pdf(file=paste0("3.", "DIFF_heatmap.pdf"),height=7,width=6)
pheatmap(afExp,                                                                      
         annotation=Type,                                                            
         color = colorRampPalette(c("blue","white","red"))(50),     
         cluster_cols =F,                                                           
         show_colnames = F,                                                         
         scale="row", 
         fontsize = 10,
         fontsize_row=6,
         fontsize_col=8,
         annotation_colors=anncolor
)
dev.off()


adjP=0.05
aflogFC=0.5
Significant=ifelse((Diff$P.Value<adjP & abs(Diff$logFC)>aflogFC), ifelse(Diff$logFC>aflogFC,"Up","Down"), "Not")

p = ggplot(Diff, aes(logFC, -log10(P.Value)))+
  geom_point(aes(col=Significant),size=3)+
  scale_color_manual(values=c(pal_npg()(2)[2], "
  labs(title = " ")+
  theme(plot.title = element_text(size = 16, hjust = 0.5, face = "bold"))+
  geom_hline(aes(yintercept=-log10(adjP)), colour="gray", linetype="twodash",size=1)+
  geom_vline(aes(xintercept=aflogFC), colour="gray", linetype="twodash",size=1)+
  geom_vline(aes(xintercept=-aflogFC), colour="gray", linetype="twodash",size=1)

p

point.Pvalue=0.01
point.logFc=2

Diff$symbol=rownames(Diff)
pdf(paste0("3.", "DIFF_vol.pdf"),width=6.5,height=6)
p=p+theme_bw()
for_label <- Diff %>% 
  filter(abs(logFC) >point.logFc & P.Value< point.Pvalue )
p+geom_point(size = 1.5, shape = 1, data = for_label) +
  ggrepel::geom_label_repel(
    aes(label = symbol),
    data = for_label,
    color="black",
    label.size =0.1
  )
dev.off()




deg=Diff
logFC_t=0
deg$g=ifelse(deg$P.Value>0.05,'stable',
             ifelse( deg$logFC > logFC_t,'UP',
                     ifelse( deg$logFC < -logFC_t,'DOWN','stable') )
)
table(deg$g)

deg$symbol=rownames(deg)
df <- bitr(unique(deg$symbol), fromType = "SYMBOL",
           toType = c( "ENTREZID"),
           OrgDb = org.Hs.eg.db  
           
           )
DEG=deg
DEG=merge(DEG,df,by.y='SYMBOL',by.x='symbol')
data_all_sort <- DEG %>% 
  arrange(desc(logFC))

geneList = data_all_sort$logFC 
names(geneList) <- data_all_sort$ENTREZID 
head(geneList)



if(!requireNamespace("R.utils",quietly = TRUE)) install.packages("R.utils",update = F,ask = F)
R.utils::setOption( "clusterProfiler.download.method",'auto' )

kk2 <- gseKEGG(geneList     = geneList,
               organism     = 'hsa',
              
               nPerm        = 10000,
               minGSSize    = 10,
               maxGSSize    = 200,
               pvalueCutoff = 0.05,
               pAdjustMethod = "none"
               
               )




GSEAOUT=as.data.frame(kk2@result)
write.table(GSEAOUT,file="4.GSEAOUT.xls",sep="\t",quote=F,col.names=T)



num=5
pdf(paste0("4.","down_GSEA.pdf"),width = 8,height = 8)
gseaplot2(kk2, geneSetID = rownames(kk2@result)[head(order(kk2@result$enrichmentScore),num)])
dev.off()
pdf(paste0("4.","up_GSEA.pdf"),width = 8,height = 8)
gseaplot2(kk2, geneSetID = rownames(kk2@result)[tail(order(kk2@result$enrichmentScore),num)])
dev.off()

num=5
pdf(paste0("4.","all_GSEA.pdf"),width = 10,height = 12)
gseaplot2(kk2, geneSetID = rownames(kk2@result)[c(head(order(kk2@result$enrichmentScore),num),tail(order(kk2@result$enrichmentScore),num))])
dev.off()



library(stringr)

pdf(paste0("4.","all_ridgeplot_GSEA.pdf"),width = 6,height = 15)
ridgeplot(kk2)
dev.off()




afdir <- paste0(getwd(),"/5.WGCNA")           
dir.create(afdir)

traitData=sample
traitData[,2]=traitData[,1]
traitData[,1]=ifelse(traitData[,1]==C,1,0)
traitData[,2]=ifelse(traitData[,2]==P,1,0)

colnames(traitData)=c(C,P)




options(stringsAsFactors = FALSE)

fpkm = read.table(paste0("1.rawexp_",GSE,".txt"),header=T,comment.char = "",check.names=F)
rownames(fpkm)=fpkm[,1]
dim(fpkm)
names(fpkm)
datExpr0 = as.data.frame(t(fpkm[,-1]))
names(datExpr0) = fpkm[,1];
rownames(datExpr0) = names(fpkm[,-1])

datExpr0


gsg = goodSamplesGenes(datExpr0, verbose = 3)
gsg$allOK

if (!gsg$allOK)
{
  
  if (sum(!gsg$goodGenes)>0)
    printFlush(paste("Removing genes:", paste(names(datExpr0)[!gsg$goodGenes], collapse = ", ")))
  if (sum(!gsg$goodSamples)>0)
    printFlush(paste("Removing samples:", paste(rownames(datExpr0)[!gsg$goodSamples], collapse = ", ")))
  
  datExpr0 = datExpr0[gsg$goodSamples, gsg$goodGenes]
}


meanFPKM=0.5  
n=nrow(datExpr0)
datExpr0[n+1,]=apply(datExpr0[c(1:nrow(datExpr0)),],2,mean)
datExpr0=datExpr0[1:n,datExpr0[n+1,] > meanFPKM]  


filtered_fpkm=t(datExpr0)
filtered_fpkm=data.frame(rownames(filtered_fpkm),filtered_fpkm)
names(filtered_fpkm)[1]="sample"
head(filtered_fpkm)
write.table(filtered_fpkm, file=paste0(afdir,"/FPKM_filter.xls"),row.names=F, col.names=T,quote=FALSE,sep="\t")

sampleTree = hclust(dist(datExpr0), method = "average")


pdf(file =paste0(afdir,"/1_sampleClustering.pdf"), width = 12, height = 9)
par(cex = 0.6)
par(mar = c(0,4,2,0))
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5,
     cex.axis = 1.5, cex.main = 2)



dev.off()














for (df in colnames(traitData)) {
  traitData[,df]=traitData[,df]/max(traitData[,df])
  print(sd(traitData[,df]))
}
max(traitData)
dim(traitData)
names(traitData)

allTraits = traitData
dim(allTraits)
names(allTraits)


fpkmSamples = rownames(datExpr0)
traitSamples =rownames(allTraits)
traitRows = match(fpkmSamples, traitSamples)
datTraits = allTraits[traitRows,]
rownames(datTraits) 
collectGarbage()


sampleTree2 = hclust(dist(datExpr0), method = "average")

traitColors = numbers2colors(datTraits, signed = FALSE)



pdf(file=paste0(afdir,"/2_Sample dendrogram and trait heatmap.pdf"),width=12,height=11)
plotDendroAndColors(sampleTree2, traitColors,
                    groupLabels = names(datTraits),
                    main = "Sample dendrogram and trait heatmap")
dev.off()












enableWGCNAThreads()

powers = c(1:30)


sft = pickSoftThreshold(datExpr0, powerVector = powers, verbose = 5)



pdf(file=paste0(afdir,"/3_Scale independence.pdf"),width=11,height=5)
par(mfrow = c(1,2))
cex1 = 0.9

plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");

abline(h=0.9,col="red")

plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
dev.off()



softPower =sft$powerEstimate

print(softPower)

adjacency = adjacency(datExpr0, power = softPower)


TOM = TOMsimilarity(adjacency);
dissTOM = 1-TOM


geneTree = hclust(as.dist(dissTOM), method = "average");



pdf(file=paste0(afdir,"/4_Gene clustering on TOM-based dissimilarity.pdf"),width=12,height=9)
plot(geneTree, xlab="", sub="", main = "Gene clustering on TOM-based dissimilarity",
     labels = FALSE, hang = 0.04)
dev.off()



minModuleSize = 40  


dynamicMods = cutreeDynamic(dendro = geneTree, distM = dissTOM,
                            deepSplit = 2, pamRespectsDendro = FALSE,
                            minClusterSize = minModuleSize);
table(dynamicMods)


dynamicColors = labels2colors(dynamicMods)
table(dynamicColors)


pdf(file=paste0(afdir,"/5_Dynamic Tree Cut.pdf"),width=8,height=6)
plotDendroAndColors(geneTree, dynamicColors, "Dynamic Tree Cut",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05,
                    main = "Gene dendrogram and module colors")
dev.off()



MEList = moduleEigengenes(datExpr0, colors = dynamicColors)
MEs = MEList$eigengenes

MEDiss = 1-cor(MEs);

METree = hclust(as.dist(MEDiss), method = "average")


pdf(file=paste0(afdir,"/6_Clustering of module eigengenes.pdf"),width=7,height=6)
plot(METree, main = "Clustering of module eigengenes",
     xlab = "", sub = "")
MEDissThres = 0.25

abline(h=MEDissThres, col = "red")
dev.off()



merge = mergeCloseModules(datExpr0, dynamicColors, cutHeight = MEDissThres, verbose = 3)

mergedColors = merge$colors

mergedMEs = merge$newMEs


pdf(file=paste0(afdir,"/7_merged dynamic.pdf"), width = 9, height = 6.5)
plotDendroAndColors(geneTree, cbind(dynamicColors, mergedColors),
                    c("Dynamic Tree Cut", "Merged dynamic"),
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
dev.off()


moduleColors = mergedColors

colorOrder = c("grey", standardColors(50))
moduleLabels = match(moduleColors, colorOrder)-1
MEs = mergedMEs









nGenes = ncol(datExpr0)
nSamples = nrow(datExpr0)

moduleTraitCor = cor(MEs, datTraits, use = "p")
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples)


pdf(file=paste0(afdir,"/8_Module-trait relationships.pdf"),width=7,height=7.5)

textMatrix = paste(signif(moduleTraitCor, 2), "\n(",
                   signif(moduleTraitPvalue, 1), ")", sep = "")

dim(textMatrix) = dim(moduleTraitCor)
par(mar = c(6, 8.5, 3, 3))


labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = names(datTraits),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.5,
               zlim = c(-1,1),
               main = paste("Module-trait relationships"))
dev.off()








modNames = substring(names(MEs), 3)

geneModuleMembership = as.data.frame(cor(datExpr0, MEs, use = "p"))
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples))

names(geneModuleMembership) = paste("MM", modNames, sep="")
names(MMPvalue) = paste("p.MM", modNames, sep="")


traitNames=names(datTraits)

geneTraitSignificance = as.data.frame(cor(datExpr0, datTraits, use = "p"))
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples))

names(geneTraitSignificance) = paste("GS.", traitNames, sep="")
names(GSPvalue) = paste("p.GS.", traitNames, sep="")

















for (trait in traitNames){
  traitColumn=match(trait,traitNames)
  
  for (module in modNames){
    column = match(module, modNames)
    moduleGenes = moduleColors==module
    
    if (nrow(geneModuleMembership[moduleGenes,]) > 1){
      
      
      pdf(file=paste(afdir,"/9_", trait, "_", module,"_Module membership vs gene significance.pdf",sep=""),width=7,height=7)
      par(mfrow = c(1,1))
      verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                         abs(geneTraitSignificance[moduleGenes, traitColumn]),
                         xlab = paste("Module Membership in", module, "module"),
                         ylab = paste("Gene significance for ",trait),
                         main = paste("Module membership vs. gene significance\n"),
                         cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)
      dev.off()
    }
  }
}


names(datExpr0)
probes = names(datExpr0)




geneInfo0 = data.frame(probes= probes,
                       moduleColor = moduleColors)

for (Tra in 1:ncol(geneTraitSignificance))
{
  oldNames = names(geneInfo0)
  geneInfo0 = data.frame(geneInfo0, geneTraitSignificance[,Tra],
                         GSPvalue[, Tra])
  names(geneInfo0) = c(oldNames,names(geneTraitSignificance)[Tra],
                       names(GSPvalue)[Tra])
}

for (mod in 1:ncol(geneModuleMembership))
{
  oldNames = names(geneInfo0)
  geneInfo0 = data.frame(geneInfo0, geneModuleMembership[,mod],
                         MMPvalue[, mod])
  names(geneInfo0) = c(oldNames,names(geneModuleMembership)[mod],
                       names(MMPvalue)[mod])
}
geneOrder =order(geneInfo0$moduleColor)
geneInfo = geneInfo0[geneOrder, ]

write.table(geneInfo, file = paste0(afdir,"/10_GS_and_MM.xls"),sep="\t",row.names=F)






nGenes = ncol(datExpr0)
nSamples = nrow(datExpr0)



plotTOM = dissTOM^7

diag(plotTOM) = NA




nSelect = 1000

set.seed(10)
select = sample(nGenes, size = nSelect)
selectTOM = dissTOM[select, select]

selectTree = hclust(as.dist(selectTOM), method = "average")
selectColors = moduleColors[select]





plotDiss = selectTOM^7
diag(plotDiss) = NA
library("gplots")

pdf(file=paste0(afdir,"/13_Network heatmap plot_selected genes.pdf"),width=9, height=9)
TOMplot(plotDiss, selectTree, selectColors, main = "Network heatmap plot, selected genes", col=gplots::colorpanel(250,'red',"orange",'lemonchiffon'))
dev.off()



pdf(file=paste0(afdir,"/14_Eigengene dendrogram and Eigengene adjacency heatmap.pdf"), width=5, height=7.5)
par(cex = 0.9)
plotEigengeneNetworks(MEs, "", marDendro = c(0,4,1,2), marHeatmap = c(3,4,1,2), cex.lab = 0.8, xLabelsAngle= 90)
dev.off()
