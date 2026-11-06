### Part --[1]----


library(Seurat)
library(infercnv)
library(dplyr)
library(GenomicRanges)

setwd("/WData/WorkSpace/xiangjun2024/LiJunHu")
load("./Data/Exprs/data.processed.RData")

data <- readRDS("./Data/Exprs/data.processed.rds")
load("/WData/WData1/Biodata/infercnv/ENSEMBL.Protein.Coding.hg38.RData")
cell.cluster <- data@meta.data$CellType
data <- AddMetaData(data, cell.cluster, col.name = "cell.cluster")


### Part --[2]----

source("/WData/WData1/Biodata/Code/scRNA/1.inferCNV-2.2.R")

inferCNV.dir <- "/WData/WorkSpace/xiangjun2024/LiJunHu/Data/inferCNV" # inferCNV.subcluster

infercnv.results <- list()

#'
for (sample_id in sample_ids) {
    cat("Running inferCNV for SampleID:", sample_id, "\n")

    sample_seurat_obj <- subset(data, subset = Response == sample_id)

    ref.group.name <- unique(sample_seurat_obj$cell.cluster)
    ref.group.name <- ref.group.name[ref.group.name != "Epithelial cells"]


    out_dir_sample <- file.path(inferCNV.dir, sample_id)
    dir.create(out_dir_sample, recursive = TRUE, showWarnings = FALSE)

    #'
    infercnv_result <- inferCNV(
        seurat.obj = sample_seurat_obj, genome.version = "hg38", ref.group.name = ref.group.name,
        out.dir = out_dir_sample, cutoff = 0.1,
        scale_data = FALSE, denoise = TRUE, HMM = FALSE, cluster_by_groups = TRUE, write_phylo = TRUE, write_expr_matrix = TRUE,
        analysis_mode = "subclusters", ## samples subclusters
        tumor_subcluster_partition_method = "leiden", num_threads = 10, save_rds = FALSE
    )

    #'
    infercnv.results[[sample_id]] <- infercnv_result

    #'
    rm(sample_seurat_obj, infercnv_result)
    gc()
}


### Part --[3]----

source("/WData/WData1/Biodata/Code/scRNA/根据inferCNV结果识别肿瘤细胞.R")

file.name <- list.files("./Data/inferCNV", full.names = FALSE, recursive = FALSE)

out.dir <- "/WData/WorkSpace/xiangjun2024/LiJunHu/Figure/4.inferCNV/"

res.df <- list()

#' 循环遍历每个inferCNV目录
for (i in 1:length(file.name)) {
    inferCNV.dir <- paste("./Data/inferCNV/", file.name[i], sep = "")
    output.file <- paste(out.dir, file.name[i], ".tumor.cell.classification.pdf", sep = "")
    res.df[[i]] <- inferCNV.TumorCell(inferCNV.dir, isDisplayNormal = TRUE, output.file)
}

inferCNV.TumorCell <- do.call(rbind, res.df)
inferCNV.TumorCell <- distinct(inferCNV.TumorCell, Cell.names, .keep_all = TRUE)
save(inferCNV.TumorCell, file = "./Data/inferCNV/inferCNV.TumorCell.RData")


### Part --[4]----

library(ggplot2)
library(RColorBrewer)
load("./Data/inferCNV/inferCNV.TumorCell.RData")
dataToPlot <- inferCNV.TumorCell

pdf("./Figure/4.inferCNV/1.0.Total.Tumor.classification.pdf", width = 7, height = 6)
ggplot(dataToPlot, aes(x = CNV.signal, y = CNV.correlation)) +
    geom_point(aes(color = Cell.index), size = 1, alpha = .8, shape = 16) +
    scale_color_manual(values = c("Malignant cells" = "#FF6347", "Unresolved" = "#b6aaaa", "Normal cells" = "#6BD76BFF")) +
    xlab("CNV signal") +
    ylab("CNV correlation") +
    theme_classic() +
    xlim(0, 0.03) +
    theme(axis.title.x = element_text(size = 16), axis.text.x = element_text(size = 12)) +
    theme(axis.title.y = element_text(size = 16), axis.text.y = element_text(size = 12))
dev.off()
