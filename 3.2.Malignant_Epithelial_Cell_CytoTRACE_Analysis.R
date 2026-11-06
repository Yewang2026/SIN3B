### Part --[1]----

setwd("/WData/WorkSpace/xiangjun2024/LiJunHu")

library(Seurat)
library(CytoTRACE2)
library(dplyr)
library(ggpubr)
rm(list = ls())

data <- readRDS("./Data/Exprs/Malignant.processed.rds")
annotation <- data.frame(phenotype = as.character(data$CellCluster)) %>%
    set_rownames(., colnames(data))


### Part --[2]----

result <- cytotrace2(data, is_seurat = TRUE, slot_type = "counts", species = "human", seed = 1234)
cytotrace.meta <- result@meta.data
cytotrace.meta$CytoTRACE2_Relative1 <- cytotrace.meta$CytoTRACE2_Relative


### Part --[3]----

col2 <- c("#C3D0E3", "#A4E1E8", "#79D5C3", "#F9CBAA", "#F2B5B5")

cytotrace.meta$CellCluster <- factor(cytotrace.meta$CellCluster)

pdf("./Figure/5.Epi/6.2.恶性细胞CytoTRACE.Boxplot.pdf", width = 5, height = 5)
ggplot(cytotrace.meta, aes(x = CellCluster, y = CytoTRACE2_Relative, fill = CellCluster)) +
    geom_boxplot(width = 0.3, position = position_dodge(width = 0.95)) +
    stat_compare_means() +
    scale_fill_manual(values = col2) +
    scale_y_continuous(limits = c(0, 1)) +
    theme_classic() +
    theme(
        legend.position = "none",
        plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 10)
    )
dev.off()
