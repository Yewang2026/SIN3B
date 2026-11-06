### Part --[1]----

rm(list = ls())
library(monocle)
library(Seurat)
library(Matrix)
library(Seurat)
library(stringr)
library(openxlsx)
library(gridExtra)
setwd("/WData/WorkSpace/xiangjun2024/LiJunHu")

data <- readRDS("./Data/Exprs/Malignant.processed.rds")
a <- as.character(names(table(data$CellCluster)))
data <- subset(data, subset = CellCluster %in% a)
counts <- GetAssayData(data, assay = "RNA", layer = "counts")

sample.info <- data@meta.data
sample.info$SampleID <- as.factor(sample.info$SampleID)
gene_anno <- data.frame(gene_short_name = rownames(data))
rownames(gene_anno) <- gene_anno$gene_short_name


### Part --[2]----


### Point 1

pd <- new("AnnotatedDataFrame", data = sample.info)
fd <- new("AnnotatedDataFrame", data = gene_anno)
cds <- newCellDataSet(counts, phenoData = pd, featureData = fd, expressionFamily = negbinomial.size())

### Point 2


cds <- estimateSizeFactors(cds)
cds <- estimateDispersions(cds)

### Point 3

cds <- detectGenes(cds, min_expr = 0.1)

expressed_genes <- row.names(subset(fData(cds), num_cells_expressed >= 10))
cds <- cds[expressed_genes, ]
disp_table <- dispersionTable(cds)
disp.genes <- subset(disp_table, mean_expression >= 0.1 & dispersion_empirical >= 1 * dispersion_fit)
ordered_genes <- disp.genes$gene_id

cds <- setOrderingFilter(cds, ordered_genes)
cds <- reduceDimension(cds, max_components = 2, method = "DDRTree", verbose = TRUE, norm_method = "log", pseudo_expr = 1, num_dim = 6)
cds <- orderCells(cds)

cds <- orderCells(cds, root_state = 4)

### Point 4.
#' 提取结果.
monocle_meta <- data.frame(t(cds@reducedDimS), cds$Pseudotime, as.character(cds$State), as.character(cds$CellCluster))
colnames(monocle_meta) <- c("C1", "C2", "Pseudotime", "State", "CellCluster")


### Part --[3]----

#' Point 1.

a <- c("cluster1", "cluster2", "cluster3", "cluster4", "cluster5")

col2 <- c("#C3D0E3", "#A4E1E8", "#79D5C3", "#F9CBAA", "#F2B5B5")
names(col2) <- a

pdf("./Figure/5.Epi/6.3.Epi拟时序Pseudotime.pdf", width = 10, height = 10)
plot_cell_trajectory(cds, color_by = "Pseudotime", size = 1, show_backbone = TRUE)
dev.off()

pdf("./Figure/5.Epi/6.4.Epi拟时序State.pdf", width = 10, height = 10)
plot_cell_trajectory(cds, color_by = "State")
dev.off()

pdf("./Figure/5.Epi/6.5.Epi拟时序CellCluster.pdf", width = 10, height = 10)
plot_cell_trajectory(cds, color_by = "CellCluster", size = 1, show_backbone = TRUE) +
    scale_color_manual(values = col2)
dev.off()

pdf("./Figure/5.Epi/6.6.Epi拟时序CellCluster1.pdf", width = 10, height = 10)
plot_cell_trajectory(cds, color_by = "CellCluster", size = 1, show_backbone = TRUE) +
    facet_wrap(~CellCluster, nrow = 2) +
    scale_color_manual(values = col2)
dev.off()


#' Point 2.

cytotrace.meta <- readRDS("./Data/聚类/CytoTRACE/Malignant.cytotrace.meta.rds")

b <- match(rownames(monocle_meta), rownames(cytotrace.meta))
cytotrace.meta <- cytotrace.meta[b, ]
monocle_meta$CytoTRACE2_Relative <- cytotrace.meta$CytoTRACE2_Relative

pdf("./Figure/5.Epi/6.7.Epi拟时序-CytoTRACE.pdf", width = 10, height = 10)
ggplot(monocle_meta, aes(x = C1, y = C2, color = CytoTRACE2_Relative)) +
    geom_point(size = 1, alpha = 1) +
    scale_color_gradientn(
        colors = c("#2A0A4A", "#FF1493", "#FF7F50", "#FFDAB9"),
        name = "Relative Score"
    ) +
    labs(title = "", x = "UMAP1", y = "UMAP2") +
    theme_minimal() +
    theme(
        plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
        axis.title = element_text(size = 12),
        axis.text = element_text(size = 10),
        axis.ticks = element_line(linewidth = 0.5),
        legend.title = element_text(size = 10),
        legend.text = element_text(size = 8),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(color = "black")
    )
dev.off()
