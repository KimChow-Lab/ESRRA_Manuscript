library(Seurat)
library(ggplot2)
setwd("/Projects/deng/Alzheimer/syn18485175")
data <- Read10X(data.dir = "FilteredData/")
#data was downloaded from syn18485175 using filtered count matrix
Mathys2019 <- CreateSeuratObject(counts = data, project = "Mathys2019")
Info=read.table("/Projects/deng/Alzheimer/syn18485175/Phenotype4Cell.txt",header=TRUE,row.names=1,sep="\t")
Mathys2019$Gender=Info$Gender
Mathys2019$Sex=Info$Sex
Mathys2019$CellType=Info$broad.cell.type
Mathys2019$Cluster=Info$pre.cluster
Mathys2019$Statues=Info$Statues
Mathys2019$OriginalSubCluster=Info$Subcluster
Mathys2019$ProjectID=Info$ProjectID
Mathys2019$BraakStage=paste0("Stage_",Info$BraakStage,sep="")
Mathys2019$Gpath=Info$gpath
Mathys2019$Amyloid=Info$amyloid
Mathys2019$Plaq=Info$plaq_n
Mathys2019$Cogdx=Info$cogdx
Mathys2019$Ceradsc=Info$ceradsc
Mathys2019$NFT=Info$nft
Mathys2019$cogn_global_lv=Info$cogn_global_lv
Mathys2019$gpath_3neocort=Info$gpath_3neocort
Mathys2019$amyloid.group=Info$amyloid.group
Mathys2019$caa_4gp=Info$caa_4gp
Mathys2019$APOEScore=Info$apoeScore
Mathys2019$Age=Info$age_death
Mathys2019$APOEGenotype=Info$APOE
Mathys2019$pathology.group=Info$pathology.group

Mathys2019 <- NormalizeData(Mathys2019, normalization.method = "LogNormalize", scale.factor = 10000)
Mathys2019 <- FindVariableFeatures(Mathys2019, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(Mathys2019)
Mathys2019 <- ScaleData(Mathys2019, features = all.genes)
Mathys2019 <- RunPCA(Mathys2019, features = VariableFeatures(object = Mathys2019))
Mathys2019<- RunTSNE(Mathys2019, reduction = "pca", dims = 1:30)
Mathys2019<- FindNeighbors(Mathys2019, reduction = "pca", dims = 1:30)
Mathys2019<- FindClusters(Mathys2019, resolution = 0.5)
# Visualization


###Figure 1E
colorPalette=colorRampPalette(brewer.pal(8, "Set3"))(length(unique(Mathys2019$CellType))) 
pdf("Human_CellType.pdf",width=5,height=4)
DimPlot(Mathys2019,group.by="CellType",label=F,raster=TRUE,cols =colorPalette)&
theme(plot.title=element_blank(),axis.title.x = element_blank(),axis.title.y = element_blank(),panel.border = element_rect(fill=NA,color="black", size=1.5, linetype="solid"))
dev.off()

DimPlot(Mathys2019, reduction = "umap", group.by="CellType")

##Figure 
pdf("Human_TargetGenes_Expr.pdf",width=5,height=4)
DotPlot(Mathys2019,group.by="CellType",features=c("ESRRA","PPARGC1A"))&RotatedAxis()&
theme(plot.title=element_text(face="italic"),axis.title.x = element_blank(),axis.title.y = element_blank(),panel.border = element_rect(fill=NA,color="black", size=1.5, linetype="solid"))
dev.off()
