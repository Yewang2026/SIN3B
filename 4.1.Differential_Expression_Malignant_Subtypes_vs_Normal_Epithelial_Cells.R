### Part --[1]----
#'

setwd("/WData/WorkSpace/xiangjun2024/LiJunHu")

rm(list = ls())
library(dplyr)
library(openxlsx)
library(tidyr)
library(Seurat)
library(limma)

load("./Data/gtf.protein.coding.RData")
load("./Data/inferCNV/inferCNV.TumorCell.RData")
data <- readRDS("./Data/Exprs/data.processed.rds")
data <- subset(data, CellType == "Epithelial cells")
Malignant <- readRDS("./Data/Exprs/Malignant.processed.rds")


### Part --[2]----

normal.cells <- inferCNV.TumorCell[inferCNV.TumorCell$Cell.index == "Unresolved", ]

normal.cells <- intersect(rownames(normal.cells), colnames(data))
NormalEpi <- subset(data, cells = normal.cells)

cat("Number of normal epithelial cells:", ncol(NormalEpi), "\n")
cat("Number of malignant cells:", ncol(Malignant), "\n")


### Part --[3]----

normal.expr <- GetAssayData(NormalEpi, assay = "RNA", layer = "data")
malignant.expr <- GetAssayData(Malignant, assay = "RNA", layer = "data")

common.genes <- intersect(rownames(normal.expr), rownames(malignant.expr))
normal.expr <- normal.expr[common.genes, ]
malignant.expr <- malignant.expr[common.genes, ]


### Part --[4]----
meta.malignant <- Malignant@meta.data
meta.malignant$CellCluster <- as.character(meta.malignant$CellCluster)

subtypes <- c("cluster1", "cluster2", "cluster3", "cluster4", "cluster5")


### Part --[5]----
for (a in subtypes) {
    cat("Now processing:", a, "vs Normal\n")

    subtype.cells <- rownames(meta.malignant)[meta.malignant$CellCluster == a]
    subtype.cells <- intersect(subtype.cells, colnames(malignant.expr))

    if (length(subtype.cells) < 3) {
        cat("Skip", a, ": too few cells\n")
        next
    }

    expr.use <- cbind(
        malignant.expr[, subtype.cells, drop = FALSE],
        normal.expr[, , drop = FALSE]
    )

    group_list <- c(
        rep(a, length(subtype.cells)),
        rep("Normal", ncol(normal.expr))
    )
    group_list <- factor(group_list, levels = c("Normal", a))

    design <- model.matrix(~ 0 + group_list)
    colnames(design) <- levels(group_list)
    rownames(design) <- colnames(expr.use)

    contrast.matrix <- makeContrasts(contrasts = paste0(a, "-Normal"), levels = design)

    fit <- lmFit(expr.use, design)
    fit2 <- contrasts.fit(fit, contrast.matrix)
    fit2 <- eBayes(fit2)

    allDiff <- topTable(fit2, adjust = "fdr", number = 200000)
    allDiff <- na.omit(allDiff)

    test_p <- allDiff$adj.P.Val <= 0.05
    test_up <- allDiff$logFC >= 0.25
    test_down <- allDiff$logFC <= -0.25

    allDiff$change <- ifelse(test_p & test_up, "UP",
        ifelse(test_p & test_down, "DOWN", "NOT")
    )

    print(table(allDiff$change))

    allDiff <- allDiff[order(allDiff$logFC, decreasing = TRUE), ]
    allDiff <- tibble::rownames_to_column(allDiff, "Symbol")

    if (exists("gtf.protein.coding")) {
        DEG <- allDiff[allDiff$Symbol %in% gtf.protein.coding, ]
    } else {
        DEG <- allDiff
    }

    save(allDiff, DEG,
        file = paste0("./Data/DEG/", a, "_vs_Normal.DEG.RData")
    )
}
