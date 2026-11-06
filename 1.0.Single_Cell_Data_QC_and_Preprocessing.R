### Part --[1]----
#' 
setwd("/WData/WorkSpace/xiangjun2024/LiJunHu")

library(Seurat)
library(harmony)
library(ggplot2)
data <- readRDS("/WData/WData1/10XGenomic/Data/Chen_CancerCell_2024_ICI/Exprs/Seurat/scRNAseq.rds")
data <- subset(data, SampleType == "Tumor")
data <- subset(data, TreatmentStage != "On")

data <- NormalizeData(data, normalization.method = "LogNormalize", scale.factor = 10000)
data <- FindVariableFeatures(data, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(data)
data <- ScaleData(data, features = all.genes, vars.to.regress = "percent.mt")
data <- RunPCA(data, features = VariableFeatures(object = data), verbose = FALSE)
data <- JackStraw(data, dims = 50, num.replicate = 100)
data <- ScoreJackStraw(data, ndims = =50)
data <- RunHarmony(data, group.by.vars = "SampleID")
data <- FindNeighbors(data, reduction = "pca", dims = 1:50)
data <- FindClusters(data, resolution = 0.5) ## 0.75.

#' UMAP.
col1 <- c(
    "#6A5ACD",   "#D2B48C",   "#F0E685FF", "#F08080",   "#5DB1DDFF", "#00CC99FF", "#9370DB",
    "#A0522D",   "#6BD76BFF", "#D595A7FF", "#97A1A799", "#87CEEB",   "#DC143C",   "#228B22",   "#7A65A5FF",
    "#BA6338FF", "#FA8072",   "#808000",   "#466983FF", "#5F9EA0",   "#1E90FF",   "#003399FF", "#D2691E",
    "#800080",   "#CDDEB7FF", "#008B8B",   "#20B2AA",   "#9932CC",   "#40E0D0",   "#D60047FF", "#DEB887",
    "#8B008B",   "#802268FF", "#8A2BE2",   "#DC944599", "#8B4513",   "#EE82EE",   "#FF00FF",   "#358DB999",
    "#008099FF", "#CE3D32FF", "#660099FF", "#E9967A",   "#9400D3",   "#99CC00FF", "#749B58FF", "#4682B4",
    "#7B68EE",   "#0099CCFF", "#0000FF",   "#32CD32",   "#5050FFFF", "#FF1493",   "#F0E68C",
    "#FF6347",   "#82581F99", "#0000CD",   "#CC9900FF", "#00D68FFF", "#2F509E99", "#E5614C99", "#98FB98",
    "#FFE4B5"
)

data <- RunUMAP(data, dims = 1:50)
umap_plot <- DimPlot(data, reduction = "umap", label = TRUE, raster = FALSE) + NoLegend()
pdf("./Figure/CellType/2.0.1.cluster.umap.pdf", width = 10, height = 10)
umap_plot + scale_color_manual(values = col1)
dev.off()

#' tSNE.
data <- RunTSNE(data, dims = 1:50)

tsne_plot <- DimPlot(data, reduction = "tsne", label = TRUE, raster = FALSE)
pdf("./Figure/CellType/2.1.1.cluster.tsne.pdf")
tsne_plot + scale_color_manual(values = col1)
dev.off()


celltype <- data.frame(ClusterID = 0:31, celltype = 0:31)

celltype[celltype$ClusterID %in% c(10,11), 2] <- "B cells"
celltype[celltype$ClusterID %in% c(8,19,24,26,28,30), 2] <- "CAFs"
celltype[celltype$ClusterID %in% c(18), 2] <- "Endothelial cells"
celltype[celltype$ClusterID %in% c(5,9,15,23,27), 2] <- "Epithelial cells"
celltype[celltype$ClusterID %in% c(7,20,21,29), 2] <- "Myeloid cells"
celltype[celltype$ClusterID %in% c(6), 2] <- "Plasma cells"
celltype[celltype$ClusterID %in% c(0:4,12,13,14,16,17,22,25,31), 2] <- "T cells"

data@meta.data$CellType <- "NA"

for (i in 1:nrow(celltype)) { ## 进行细胞注释时要注意选择的分辨率.
    data@meta.data[which(data@meta.data$RNA_snn_res.0.5 == celltype$ClusterID[i]), "CellType"] <- celltype$celltype[i]
}
Idents(data) <- data$CellType

col2 <- c("#F08080", "#D2B48C", "#F0E685FF", "#6A5ACD", "#5DB1DDFF", "#00CC99FF", "#808000")

a <- c("T cells", "B cells", "Plasma cells", "Myeloid cells", "CAFs", "Endothelial cells", "Epithelial cells")
data$CellType <- factor(data$CellType, levels = a)
Idents(data) <- data$CellType

pdf("./Figure/CellType/2.9.celltype.umap.pdf", width = 10, height = 10)
DimPlot(data, reduction = "umap", cols = col2, label = TRUE, pt.size = 1, raster = FALSE) + NoLegend()
dev.off()

pdf("./Figure/CellType/2.10.celltype.tsne.pdf", width = 10, height = 10)
DimPlot(data, reduction = "tsne", cols = col2, label = TRUE, pt.size = 1, raster = FALSE) # + NoLegend()
dev.off()
