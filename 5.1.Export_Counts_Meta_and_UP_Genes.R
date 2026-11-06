setwd("/WData/WorkSpace/xiangjun2024/LiJunHu")

library(Matrix)
library(Seurat)

dir.create("./Data/Geneformer", showWarnings = FALSE, recursive = TRUE)
dir.create("./Data/Geneformer/up_genes", showWarnings = FALSE, recursive = TRUE)

# =========================
# Part 1.
# =========================
data <- readRDS("./Data/Exprs/Malignant.processed.rds")

counts <- GetAssayData(data, layer = "counts", assay = "RNA")
counts <- as(counts, "dgCMatrix")

Matrix::writeMM(counts, file = "./Data/Geneformer/malignant_counts.mtx")

write.csv(
    data.frame(gene_name = rownames(counts)),
    file = "./Data/Geneformer/malignant_genes.csv",
    row.names = FALSE,
    quote = FALSE
)

meta <- data@meta.data[colnames(counts), , drop = FALSE]
write.csv(
    meta,
    file = "./Data/Geneformer/malignant_meta.csv",
    row.names = TRUE
)

cat("Matrix / genes / meta exported.\n")

# =========================
# Part 2.
# =========================
deg_dir <- "./Data/DEG"
out_dir <- "./Data/Geneformer/up_genes"
groups <- c("cluster1", "cluster2", "cluster3", "cluster4", "cluster5")

for (g in groups) {
    f <- file.path(deg_dir, paste0(g, "_vs_Normal.DEG.RData"))
    load(f)

    deg_df <- DEG

    if (!all(c("Symbol", "change") %in% colnames(deg_df))) {
        stop(paste("File", f, "does not contain Symbol/change columns"))
    }

    up_deg <- deg_df[deg_df$change == "UP", , drop = FALSE]
    up_deg <- up_deg[!duplicated(up_deg$Symbol), , drop = FALSE]
    up_deg <- up_deg[order(up_deg$logFC, decreasing = TRUE), ]
    up_deg <- head(up_deg, 10)

    write.csv(
        up_deg,
        file = file.path(out_dir, paste0(g, "_UP_DEG.csv")),
        row.names = FALSE,
        quote = FALSE
    )

    write.csv(
        data.frame(Symbol = up_deg$Symbol),
        file = file.path(out_dir, paste0(g, "_UP_genes.csv")),
        row.names = FALSE,
        quote = FALSE
    )

    cat(g, "UP genes:", nrow(up_deg), "\n")
}

cat("UP genes exported for cluster1-cluster5.\n")
