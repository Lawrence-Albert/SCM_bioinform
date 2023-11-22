## ---------------------------
## Author: Yukun Li
##
## Date Created: 2023-07-10
##
## Copyright (c) Yukun Li, 2023
## Email: lorenli@mail.ccmu.edu.cn

## Purpose of script:
##
## Comparison of single-cell analysis before and after normalization. 
## Overview of the mitochondria/inflammation-related pathway of interest at single-cell resolution.

library(Seurat)
library(scRNAtoolVis)

rm(list = ls())
scdata <- Read10X(data.dir = "S0/s3t0/")


scobj <- CreateSeuratObject(counts = scdata, 
                            project = "s3", 
                            min.cells = 3, 
                            min.features = 200)

metadata = scobj@meta.data
scobj@meta.data$ID = "s3" 
saveRDS(scobj,file = "s3.rds")


data1 <- readRDS(file = "c1.rds")
data2 <- readRDS(file = "c2.rds")

data3 <- readRDS(file = "s1.rds")
data4 <- readRDS(file = "s2.rds")
data5 <- readRDS(file = "s3.rds")

data <- list(data1,data2,data3,data4,data5)

scobj <- merge(x=data[[1]], y = data[-1], project = "scm")

rm(list = ls(pattern="data.*"))
scobj[["percent.mt"]] <- PercentageFeatureSet(scobj, pattern = "^MT-")
VlnPlot(scobj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

scobj <- subset(scobj, subset = nFeature_RNA > 500 & 
                  nFeature_RNA < 2500 & 
                  percent.mt < 20)

scobj <- NormalizeData(scobj)
scobj <- FindVariableFeatures(scobj, selection.method = "vst", nfeatures = 2000)
scobj <- ScaleData(scobj, features = rownames(scobj))
scobj <- RunPCA(scobj, features = VariableFeatures(object = scobj),reduction.name = "pca")


scobj <- RunUMAP(scobj,reduction = "pca", dims = 1:10, reduction.name = "umap_naive")
pdf(file = "rawumap.pdf",width =6,height = 5)
DimPlot(scobj, reduction = "umap_naive",group.by = "ID")
dev.off()



scobj <- RunHarmony(scobj,reduction = "pca",group.by.vars = "ID",reduction.save = "harmony")
scobj <- RunUMAP(scobj, reduction = "harmony", dims = 1:30,reduction.name = "umap")
pdf(file = "HARumap.pdf",width =6,height = 5)
DimPlot(scobj, reduction = "umap",group.by = "ID")
dev.off()

scobj <- FindNeighbors(scobj, reduction = "harmony", dims = 1:30)

scobj <- FindClusters(scobj, resolution = seq(0.2,1.2,0.1))
library(clustree)


pdf(file = "clustertree.pdf",width =10,height =10)
clustree(scobj)
dev.off()

scobj@meta.data$seurat_clusters <- scobj@meta.data$RNA_snn_res.0.3
Idents(scobj) <- "seurat_clusters"

DimPlot(scobj, reduction = "umap", label = T)

scobj@assays$RNA@scale.data <- matrix()
scobj@reductions$umap_naive <- NULL
saveRDS(scobj,file = "hamony_seurat_unannotaion.rds")

rm(list = ls())
library(Seurat)
scobj <- readRDS(file = "hamony_seurat_unannotaion.rds")


pdf(file = "4.UNANumap.pdf",width =6,height = 5)
DimPlot(scobj, reduction = "umap", label = T)
dev.off()


rna.data.average = AverageExpression(scobj,group.by = "seurat_clusters")
rna.data.average = round(rna.data.average$RNA, 2)
write.table(rna.data.average, "CELLiD_input.txt", quote = F, col.names = F, row.names = T, sep="\t")

all_markers <- FindAllMarkers(object = scobj)
saveRDS(all_markers,file = "all_markers.rds")
library(dplyr)
top_markers <- all_markers %>%
  group_by(cluster) %>%
  arrange(desc(avg_log2FC))%>%
  slice_head(n = 20)
write.csv(top_markers,file = "top20.csv")

marker_genes <- c("IGHG4","IGLC2","IGLC3","IGKV3-15")
VlnPlot(scobj, features = marker_genes)
FeaturePlot(scobj, features = marker_genes, order = TRUE,ncol=3)
library(Nebulosa)
plot_density(scobj,features = marker_genes) + plot_layout(ncol = 2)


head(Idents(scobj))
Idents(scobj) <- "seurat_clusters"


scobj <- RenameIdents(scobj,
                      "0"="Naive B",
                      "1"="CD4 Naive T", 
                      "2"="CD14 Mono", 
                      "3"= "NK", 
                      "4"= "CD4 Memory T", 
                      "5"= "Memory B",
                      "6"= "CD16 Mono", 
                      "7"= "Plasma cell", 
                      "8"= "Megakaryocyte",
                      "9"= "DC", 
                      "10"= "RBC", 
                      "11"= "pDC",
                      "12"= "Neutrophil"
)
head(Idents(scobj))

scobj@meta.data$celltype = Idents(scobj)


pdf(file = "5.ANumap.pdf",width =6.5,height = 5)
DimPlot(scobj, reduction = "umap", label = T,
        label.size = 3)
dev.off()

saveRDS(scobj,file = "hamony_seurat_ANnotaion.rds")



rm(list = ls())
scobj <- readRDS(file = "hamony_seurat_ANnotaion.rds")

markers <- c("FBXO7","PGS1","BCS1L","LYRM7","MRPS31","TIMMDC1")

pdf(file = "6.ANumapgroup.pdf",width =2.9,height = 5)
clusterCornerAxes(object = scobj,
                  reduction = 'umap',
                  clusterCol = "celltype",
                  groupFacet = 'group',
                  noSplit = F,
                  nrow=2,
                  cellLabel = T,
                  cornerTextSize=2,
                  cellLabelSize =1.6 ,
                  show.legend = F,
                  aspect.ratio = 1)
dev.off()



library(scRNAtoolVis)
library(scCustomize)

pdf(file = "7.SCMumap.pdf",width =16,height = 2.9)
FeatureCornerAxes(object = subset(scobj,group=="scm"),
                  reduction = 'umap',
                   groupFacet = NULL,
                  relLength = 0.5,relDist = 0.2,
                  minExp = 0,
                  maxExp = 2.5,
                  show.legend = F,
                  features = markers)
dev.off()

pdf(file = "7.CONumap.pdf",width =16,height = 2.9)
FeatureCornerAxes(object = subset(scobj,group=="con"),
                  reduction = 'umap',
                  groupFacet = NULL,
                  relLength = 0.5,relDist = 0.2,
                  minExp = 0,
                  maxExp = 4,
                  show.legend = F,
                  features = markers)
dev.off()


pdf(file = "8.dot.pdf",width =7,height = 6)
DotPlot(scobj, features = markers, 
        cols = c("blue", "red"), 
        dot.scale = 6, 
        split.by = "group")+RotatedAxis()
dev.off()




df <- FetchData(object = scobj, vars = c("FBXO7","celltype","group"))

library(devtools)
library(ggunchained) 
library(ggplot2)
library(ggpubr)

colnames(df) <- c('gene','sample','group')

pdf(file = "8.vlo.pdf",width =14,height = 6)
ggplot(df, aes(x = sample,y = gene, fill = group))+
  geom_split_violin(colour=NA, scale = 'width')+
  scale_fill_manual(values = c("limegreen", "navy"))+
  theme_bw()+
  labs(title = "FBXO7", y="Expression", x = "Celltype")+
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.y = element_text(size = 10, color="black"),
        panel.background = element_blank(),
        axis.text.x = element_text(size = 10, color="black",angle = 90),
        axis.title.y = element_text(size = 12, color="black"))+
  stat_summary(fun = mean,
               fun.min = function(x){quantile(x)[2]},
               fun.max = function(x){quantile(x)[4]},
               geom = "pointrange",
               size=0.3,
               position = position_dodge(width = 0.5),
               color='white')+
  ylim(0,2.5)+
  theme(panel.background = element_rect(fill = "
  stat_compare_means(aes(group = group), label = "p.signif")+
  RotatedAxis()

dev.off()


rm(list = ls())
scobj <- readRDS(file = "RDS/hamony_seurat_ANnotaion.rds")

library(clusterProfiler)
genesets <- read.gmt("human_mitoinflam.gmt")

unique(genesets$term)


signatures <- split(genesets$gene,genesets$term)
signatures <- signatures[-55]

scobjGSEA <- irGSEA.score(object = scobj, assay = "RNA", slot = "data", seeds = 123, ncores = 1,
                      msigdb=F, 
                      custom = T,geneset = signatures, 
                      method =  c("AUCell", "UCell", "singscore", 
                                  "ssgsea"), kcdf = 'Gaussian')

scobjGSEAcon <- irGSEA.score(object = subset(scobj,group=="con"),
                             assay = "RNA", slot = "data", seeds = 123, ncores = 1,
                             msigdb=F, 
                             custom = T,geneset = signatures, 
                             method =  c("AUCell", "UCell", "singscore", 
                                         "ssgsea"), kcdf = 'Gaussian')


scobjGSEAscm <- irGSEA.score(object = subset(scobj,group=="scm"), 
                             assay = "RNA", slot = "data", seeds = 123, ncores = 1,
                             msigdb=F, 
                             custom = T,geneset = signatures, 
                             method =  c("AUCell", "UCell", "singscore", 
                                         "ssgsea"), kcdf = 'Gaussian')

saveRDS(scobjGSEA,file = "RDS/GSEA.rds")
saveRDS(scobjGSEAcon,file = "RDS/GSEAcon.rds")
saveRDS(scobjGSEAscm,file = "RDS/GSEAscm.rds")



scobjGSEA <- readRDS(file = "RDS/GSEA.rds")
scobjGSEAcon <- readRDS(file = "RDS/GSEAcon.rds")
scobjGSEAscm <- readRDS(file = "RDS/GSEAscm.rds")
result.dge.group <- readRDS(file = "RDS/RESULTDGE.rds")

scobjGSEA@meta.data$group26 <- paste(scobjGSEA@meta.data$celltype, 
                                     scobjGSEA@meta.data$group, 
                                     sep = ".")
scobjGSEA@meta.data$group26 <- factor(scobjGSEA@meta.data$group26)


result.dge.group <- irGSEA.integrate(object = scobjGSEA, 
                                     group.by = "group26",
                                     metadata = NULL, col.name = NULL,
                                     method = c("AUCell", "UCell", "singscore", 
                                                "ssgsea"))
saveRDS(result.dge.group, file = "RDS/RESULTDGE.rds")















p2 <- irGSEA.heatmap(object = result.dge.group, 
                     method = "ssgsea",
                     top = 54, 
                     show.geneset = NULL,
                     heatmap.width = 33,
                     heatmap.heigh = 25)

pdf(file = "9.2heatmap.pdf",width =15,height = 10)
p2
dev.off()

result.dge.group2 <- irGSEA.integrate(object = scobjGSEA, 
                                     group.by = "group",
                                     metadata = NULL, col.name = NULL,
                                     method = c("AUCell", "UCell", "singscore", 
                                                "ssgsea"))

p3 <- irGSEA.heatmap(object = result.dge.group2, 
                     method = "ssgsea",
                     top = 54, 
                     show.geneset = NULL,
                     heatmap.width = 10,
                     heatmap.heigh = 25)

pdf(file = "9.3heatmap.pdf",width =6,height = 10)
p3
dev.off()


pdf(file = "9.SCMplot.pdf",width =6,height = 5)
DimPlot(object = subset(scobj,group=="scm"), 
        reduction = "umap", label = T,label.box = T,label.size =2.5)
dev.off()

pdf(file = "9.plot1.pdf",width =6,height = 5)
irGSEA.density.scatterplot(object = scobjGSEAscm,
                                    method = "UCell",
                                    show.geneset = "Oxidative Phosphorylation",
                                    reduction = "umap")
dev.off()


pdf(file = "9.plot2.pdf",width =6,height = 5)
irGSEA.density.scatterplot(object = scobjGSEAscm,
                           method = "UCell",
                           show.geneset = "Reactive Oxygen Species Pathway",
                           reduction = "umap")
dev.off()


pdf(file = "9.plot3.pdf",width =6,height = 5)
irGSEA.density.scatterplot(object = scobjGSEAscm,
                           method = "UCell",
                           show.geneset = "Inflammatory Response",
                           reduction = "umap")
dev.off()

pdf(file = "9.plot3ridge.pdf",width =8,height = 5)
irGSEA.ridgeplot(object = scobjGSEAscm,
                      method = "AUCell",
                      show.geneset = "Inflammatory Response")
dev.off()

pdf(file = "9.plot4.pdf",width =6,height = 5)
irGSEA.density.scatterplot(object = scobjGSEAscm,
                           method = "UCell",
                           show.geneset = "Mitophagy",
                           reduction = "umap")
dev.off()

pdf(file = "9.plot5.pdf",width =6,height = 5)
irGSEA.density.scatterplot(object = scobjGSEAscm,
                           method = "UCell",
                           show.geneset = "Mitochondrial Biogenesis",
                           reduction = "umap")
dev.off()

pdf(file = "9.plot6.pdf",width =6,height = 5)
irGSEA.density.scatterplot(object = scobjGSEAscm,
                           method = "UCell",
                           show.geneset = "Fission",
                           reduction = "umap")
dev.off()


