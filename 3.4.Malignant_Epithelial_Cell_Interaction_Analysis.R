library(CellChat)
library(patchwork)
library(ggplot2)
library(dplyr)
library(Seurat)
library(stringr)

setwd("/WData/WorkSpace/xiangjun2024/LiJunHu")


### Part --[1]----

Tcell <- readRDS("./Data/Exprs/data.processed.rds")
Malignant <- readRDS("./Data/Exprs/Malignant.processed.rds")

Tcell <- subset(Tcell, CellType != "Epithelial cells")
Tcells.counts <- GetAssayData(Tcell, assay = "RNA", layer = "data") ## 需要产生标准化后的数据,也就是data里面的数据.
Malignant.counts <- GetAssayData(Malignant, assay = "RNA", layer = "data") ## 需要产生标准化后的数据,也就是data里面的数据.

Tcells.meta.data <- Tcell@meta.data
Tcells.meta.mat <- data.frame(cbind(cell = Tcells.meta.data$CellID, cell_type = as.character(Tcells.meta.data$CellType)))

Malignant.meta.data <- Malignant@meta.data
Malignant.meta.mat <- data.frame(cbind(cell = Malignant.meta.data$CellID, cell_type = as.character(Malignant.meta.data$CellCluster)))

#' Point 2.

counts.mat <- cbind(Malignant.counts, Tcells.counts)

meta.mat <- rbind(Malignant.meta.mat, Tcells.meta.mat)
meta.mat[, 2] <- str_replace_all(meta.mat[, 2], "\\.", "_")
rownames(meta.mat) <- meta.mat[, 1]


### Part --[2]----

cellchat <- createCellChat(
    object = as.matrix(counts.mat),
    meta = meta.mat,
    group.by = "cell_type"
)

CellChatDB <- CellChatDB.human

cellchat@DB <- CellChatDB

cellchat <- subsetData(cellchat)
cellchat <- identifyOverExpressedGenes(cellchat) #
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- computeCommunProb(cellchat,
    population.size = TRUE,
    type = "truncatedMean",
    trim = 0.1
)
cellchat <- filterCommunication(cellchat, min.cells = 5)
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)

saveRDS(cellchat, file = "./Data/聚类/CCI/cellchat.rds")


### Part --[3]----

rownames(cellchat@net$count)
cell.order <- c(
    "cluster1", "cluster2", "cluster3", "cluster4", "cluster5",
    "T cells", "B cells", "Plasma cells", "Myeloid cells", "CAFs", "Endothelial cells"
)
cellchat@idents <- factor(cellchat@idents, levels = cell.order)

idx <- match(cell.order, rownames(cellchat@net$count))
cellchat@net$count <- cellchat@net$count[idx, idx]
cellchat@net$weight <- cellchat@net$weight[idx, idx]
cellchat@net$pval <- cellchat@net$pval[idx, idx]

cellchat@netP$netP <- lapply(cellchat@netP$netP, function(mat) {
    mat[idx, idx]
})

groupSize <- as.numeric(table(cellchat@idents))
col1 <- c(
    "#C3D0E3", "#A4E1E8", "#79D5C3", "#F9CBAA", "#F2B5B5",
    "#F08080", "#D2B48C", "#F0E685FF", "#6A5ACD", "#5DB1DDFF", "#00CC99FF"
)

pdf("./Figure/5.Epi/7.1.细胞与上皮细胞互作数量图.pdf", width = 8, height = 8)
netVisual_circle(cellchat@net$count,
    vertex.weight = groupSize,
    weight.scale = TRUE,
    label.edge = FALSE,
    color.use = col1,
    title.name = "Number of Interactions"
)
dev.off()

mat <- cellchat@net$count
num_rows <- nrow(mat)
num_cols <- 2
num_plot_rows <- ceiling(num_rows / num_cols)
png("./Figure/5.Epi/7.1.细胞与上皮细胞互作数量图.png", width = 800, height = num_plot_rows * 400)
par(mfrow = c(num_plot_rows, num_cols), xpd = TRUE)
for (i in 1:num_rows) {
    mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
    mat2[i, ] <- mat[i, ]
    netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, edge.weight.max = max(mat), title.name = rownames(mat)[i])
}
dev.off()

pdf("./Figure/5.Epi/7.2.细胞与上皮细胞互作强度图.pdf", width = 8, height = 8)
netVisual_circle(cellchat@net$weight,
    vertex.weight = groupSize,
    weight.scale = TRUE,
    label.edge = FALSE,
    color.use = col1,
    title.name = "Interaction Strength"
)
dev.off()

mat <- cellchat@net$weight
num_rows <- nrow(mat)
num_cols <- 2
num_plot_rows <- ceiling(num_rows / num_cols)
png("./Figure/5.Epi/7.2.细胞与上皮细胞互作强度图.png", width = 800, height = num_plot_rows * 400)
par(mfrow = c(num_plot_rows, num_cols), xpd = TRUE)
for (i in 1:num_rows) {
    mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
    mat2[i, ] <- mat[i, ]
    netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, edge.weight.max = max(mat), title.name = rownames(mat)[i])
}
dev.off()

pdf("./Figure/5.Epi/7.3.细胞与上皮细胞整体通讯热图.pdf", width = 12, height = 8)
netVisual_heatmap(cellchat, color.heatmap = "Reds", color.use = col1, row.show = cell.order, col.show = cell.order)
dev.off()

cellchat@netP$pathways

pathways.all <- cellchat@netP$pathways

base.dir <- "./Figure/5.Epi/CellChat"

library(circlize)

for (pw in pathways.all) {
    pathways.show <- pw
    pw.safe <- gsub("[^A-Za-z0-9]+", "_", pathways.show)
    out.dir <- file.path(base.dir, pw.safe)
    if (!dir.exists(out.dir)) dir.create(out.dir, recursive = TRUE)

    pdf(file.path(out.dir, paste0("7.4.", pw.safe, ".circle.pdf")), width = 8, height = 8)
    netVisual_aggregate(cellchat, signaling = pathways.show, color.use = col1, layout = "circle")
    dev.off()

    tryCatch(
        {
            pdf(file.path(out.dir, paste0("7.5.", pw.safe, ".chord.pdf")), width = 16, height = 8)
            circos.clear()
            netVisual_individual(cellchat, signaling = pathways.show, layout = "chord", color.use = col1, nCol = 2)
            dev.off()
        },
        error = function(e) {
            message("  [跳过 chord 图] ", pathways.show, " : ", e$message)
            circos.clear()
            try(dev.off(), silent = TRUE)
        }
    )

    pdf(file.path(out.dir, paste0("7.6.", pw.safe, ".hierarchy.pdf")), width = 12, height = 8)
    netVisual_aggregate(cellchat, signaling = pathways.show, layout = "hierarchy", color.use = col1, vertex.receiver = c(1, 2))
    dev.off()
}


for (pw in pathways.all) {
    pathways.show <- pw
    pw.safe <- gsub("[^A-Za-z0-9]+", "_", pathways.show)
    out.dir <- file.path(base.dir, pw.safe)

    ComplexHeatmap::ht_opt(RESET = TRUE)

    pdf(file.path(out.dir, paste0("7.7.", pw.safe, ".heatmap.pdf")), width = 12, height = 8)

    ht <- netVisual_heatmap(
        object        = cellchat,
        signaling     = pathways.show,
        color.use     = col1,
        row.show      = cell.order,
        col.show      = cell.order,
        color.heatmap = c("#2166AC", "#B2182B")
    )

    ComplexHeatmap::draw(
        ht,
        heatmap_legend_side    = "right",
        annotation_legend_side = "right"
    )
    dev.off()
}
