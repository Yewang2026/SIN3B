### Part -[1]-
library(dplyr)
library(logger)
library(httr)
library(jsonlite)
library(readr)
setwd("E:/其他/LiJunHu")

source("D:/BaiduSyncdisk/Jun/Central/Code/scRNAseq/scTherapy/R/predict_compounds.R")

gene_list <- "D:/BaiduSyncdisk/Jun/Central/Code/scRNAseq/scTherapy/geneinfo_beta_input.txt"
gene_info <- data.table::fread(gene_list) %>% as.data.frame()


### Part -[2]-
load("./Data/DEG/cluster1_vs_Normal.DEG.RData")
rownames(DEG) <- DEG$Symbol
malignant_cells_DEG <- DEG[, c(2, 6)]
colnames(malignant_cells_DEG) <- c("avg_log2FC", "p_val_adj")

# filter DEG
DEG_malignant <- malignant_cells_DEG %>%
    mutate(gene_symbol = rownames(.)) %>%
    inner_join(gene_info, by = "gene_symbol") %>%
    filter((p_val_adj <= 0.05 & (avg_log2FC > 1 | avg_log2FC < -1)) | (avg_log2FC > -0.1 & avg_log2FC < 0.1))
DEG_malignant_list <- setNames(as.list(DEG_malignant$avg_log2FC), DEG_malignant$gene_symbol)

# predict monotherapies for malignant cluster
monotherapy_drugs <- predict_drugs(DEG_malignant_list)
cluster1.drugs <- monotherapy_drugs
saveRDS(cluster1.drugs, file = "./Data/scTherapy/cluster1.drugs.rds")


### Part -[3]-
load("./Data/DEG/cluster2_vs_Normal.DEG.RData")
rownames(DEG) <- DEG$Symbol
malignant_cells_DEG <- DEG[, c(2, 6)]
colnames(malignant_cells_DEG) <- c("avg_log2FC", "p_val_adj")

# filter DEG
DEG_malignant <- malignant_cells_DEG %>%
    mutate(gene_symbol = rownames(.)) %>%
    inner_join(gene_info, by = "gene_symbol") %>%
    filter((p_val_adj <= 0.05 & (avg_log2FC > 1 | avg_log2FC < -1)) | (avg_log2FC > -0.1 & avg_log2FC < 0.1))
DEG_malignant_list <- setNames(as.list(DEG_malignant$avg_log2FC), DEG_malignant$gene_symbol)

# predict monotherapies for malignant cluster
monotherapy_drugs <- predict_drugs(DEG_malignant_list)
cluster2.drugs <- monotherapy_drugs
saveRDS(cluster2.drugs, file = "./Data/scTherapy/cluster2.drugs.rds")


### Part -[4]-
load("./Data/DEG/cluster3_vs_Normal.DEG.RData")
rownames(DEG) <- DEG$Symbol
malignant_cells_DEG <- DEG[, c(2, 6)]
colnames(malignant_cells_DEG) <- c("avg_log2FC", "p_val_adj")

# filter DEG
DEG_malignant <- malignant_cells_DEG %>%
    mutate(gene_symbol = rownames(.)) %>%
    inner_join(gene_info, by = "gene_symbol") %>%
    filter((p_val_adj <= 0.05 & (avg_log2FC > 1 | avg_log2FC < -1)) | (avg_log2FC > -0.1 & avg_log2FC < 0.1))
DEG_malignant_list <- setNames(as.list(DEG_malignant$avg_log2FC), DEG_malignant$gene_symbol)

# predict monotherapies for malignant cluster
monotherapy_drugs <- predict_drugs(DEG_malignant_list)
cluster3.drugs <- monotherapy_drugs
saveRDS(cluster3.drugs, file = "./Data/scTherapy/cluster3.drugs.rds")


### Part -[5]-
load("./Data/DEG/cluster4_vs_Normal.DEG.RData")
rownames(DEG) <- DEG$Symbol
malignant_cells_DEG <- DEG[, c(2, 6)]
colnames(malignant_cells_DEG) <- c("avg_log2FC", "p_val_adj")

# filter DEG
DEG_malignant <- malignant_cells_DEG %>%
    mutate(gene_symbol = rownames(.)) %>%
    inner_join(gene_info, by = "gene_symbol") %>%
    filter((p_val_adj <= 0.05 & (avg_log2FC > 1 | avg_log2FC < -1)) | (avg_log2FC > -0.1 & avg_log2FC < 0.1))
DEG_malignant_list <- setNames(as.list(DEG_malignant$avg_log2FC), DEG_malignant$gene_symbol)

# predict monotherapies for malignant cluster
monotherapy_drugs <- predict_drugs(DEG_malignant_list)
cluster4.drugs <- monotherapy_drugs
saveRDS(cluster4.drugs, file = "./Data/scTherapy/cluster4.drugs.rds")


### Part -[6]-
load("./Data/DEG/cluster5_vs_Normal.DEG.RData")
rownames(DEG) <- DEG$Symbol
malignant_cells_DEG <- DEG[, c(2, 6)]
colnames(malignant_cells_DEG) <- c("avg_log2FC", "p_val_adj")

# filter DEG
DEG_malignant <- malignant_cells_DEG %>%
    mutate(gene_symbol = rownames(.)) %>%
    inner_join(gene_info, by = "gene_symbol") %>%
    filter((p_val_adj <= 0.05 & (avg_log2FC > 1 | avg_log2FC < -1)) | (avg_log2FC > -0.1 & avg_log2FC < 0.1))
DEG_malignant_list <- setNames(as.list(DEG_malignant$avg_log2FC), DEG_malignant$gene_symbol)

# predict monotherapies for malignant cluster
monotherapy_drugs <- predict_drugs(DEG_malignant_list)
cluster5.drugs <- monotherapy_drugs
saveRDS(cluster5.drugs, file = "./Data/scTherapy/cluster5.drugs.rds")


### Part -[7]-
cluster1.drugs$Subtype <- "cluster1"
cluster2.drugs$Subtype <- "cluster2"
cluster3.drugs$Subtype <- "cluster3"
cluster4.drugs$Subtype <- "cluster4"
cluster5.drugs$Subtype <- "cluster5"

scTherapy.drugs <- rbind(cluster1.drugs, cluster2.drugs, cluster3.drugs, cluster4.drugs, cluster5.drugs)
saveRDS(scTherapy.drugs, file = "./Data/scTherapy/scTherapy.drugs.rds")
