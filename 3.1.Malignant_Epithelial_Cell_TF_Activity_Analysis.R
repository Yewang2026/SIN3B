### Part --[1]----
# conda create -n pyscenic_3.7 python=3.7
conda activate pyscenic_3.7

setwd("/WData/WorkSpace/xiangjun2024/LiJunHu/Data/SCENIC")

library(Seurat)
library(tidyverse)
library(patchwork)
library(SCENIC)
library(openxlsx)
library(SCopeLoomR)
library(AUCell)
library(pheatmap)
library(circlize)

source("/WData/WData1/Biodata/Code/scRNA/pySCENIC.R")

seurat.obj <- readRDS("/WData/WorkSpace/xiangjun2024/LiJunHu/Data/Exprs/Malignant.processed.rds")

Idents(seurat.obj) <- "CellCluster"
res.scenic <- pySCENIC.Seurat(seurat.obj, layer = "count", genome_version = "hg38v10",
  output.loom = "/WData/WorkSpace/xiangjun2024/LiJunHu/Data/SCENIC/sample_SCENIC.loom"
)
saveRDS(res.scenic, file = "/WData/WorkSpace/xiangjun2024/LiJunHu/Data/SCENIC/res.scenic.rds")


### Part --[3]----

#' Point 1.

loom <- open_loom("sample_SCENIC.loom")
regulons_incidMat <- get_regulons(loom, column.attr.name = "Regulons")
regulons <- regulonsToGeneLists(regulons_incidMat)
class(regulons)

regulonAUC <- get_regulons_AUC(loom, column.attr.name = "RegulonsAUC")
regulonAucThresholds <- get_regulon_thresholds(loom)

#' Point 2.

cellinfo <- seurat.obj@meta.data[, c("CellCluster", "nFeature_RNA", "nCount_RNA")] # 细胞meta信息
colnames(cellinfo) <- c("celltype", "nGene", "nUMI")
cellTypes <- as.data.frame(subset(cellinfo, select = "celltype"))
selectedResolution <- "celltype"
sub_regulonAUC <- regulonAUC
rss <- calcRSS(AUC = getAUC(sub_regulonAUC), cellAnnotation = cellTypes[colnames(sub_regulonAUC), selectedResolution])
rss <- na.omit(rss)

#' Point 3.

a <- c(names(table(cellTypes[, 1])))
a <- match(a, colnames(rss))
rss <- rss[, a]

rss.zscore <- t(apply(rss, 1, scale)) 
colnames(rss.zscore) <- colnames(rss)
rss.zscore <- data.frame(rss.zscore)

colanno <- data.frame(
  CellID = c("cluster1", "cluster2", "cluster3", "cluster4","cluster5"),
  CellCluster = c("cluster1", "cluster2", "cluster3", "cluster4", "cluster5")
)
rownames(colanno) <- colanno$CellID
colanno$CellID <- NULL
colanno$CellCluster <- factor(colanno$CellCluster, levels = c("cluster1", "cluster2", "cluster3", "cluster4", "cluster5"))

CellCluster <- c("#C3D0E3", "#A4E1E8", "#79D5C3", "#F9CBAA", "#F2B5B5")
names(CellCluster) <- c("cluster1", "cluster2", "cluster3", "cluster4", "cluster5")
ann_colors <- list(CellCluster = CellCluster) 

col_fun <- colorRamp2(
  c(-1, -0.5, 0, 0.5, 1),
  rev(c("#E07A5F", "#F4A261", "#F6D186", "#C5E1A5", "#81C784")) 
)

pdf("/WData/WorkSpace/xiangjun2024/LiJunHu/Figure/5.Epi/6.1.上皮恶性细胞TF活性.pdf", width = 5, height = 6)
pheatmap(as.matrix(rss.zscore[1:10,]),
  scale = "row", 
  show_colnames = F, show_rownames = T,
  #color = col_fun(seq(-1, 1, length.out = 100)),
  color = colorRampPalette(c("#f0f0f5", "white", "darkred"))(1000),
  annotation_col = colanno, annotation_colors = ann_colors,
  cluster_cols = F, cluster_rows = F,
  border_color = NA, 
  fontsize_col = 8, fontsize_row = 12,
  main = "TF activity"
)
dev.off()

