### Part --[1]----

setwd("/WData/WorkSpace/xiangjun2024/LiJunHu")
library(openxlsx)
library(dplyr)
library(ggplot2)
library(tidyr)

Knockdown.gene <- read.xlsx("./Data/Geneformer/Knockdown.gene.xlsx")
Knockdown.gene <- Knockdown.gene[Knockdown.gene$sum == 5, ]
Knockdown.gene <- Knockdown.gene$GeneSymbol

Geneeffect <- read.csv("/WData/WData2/xiangjun/ChangZeWen/Data/CRISPRGeneEffect.csv", check.names = FALSE, row.names = 1)
colnames(Geneeffect) <- gsub(" *\\(.*\\)$", "", colnames(Geneeffect))
Geneeffect <- Geneeffect[, colnames(Geneeffect) %in% Knockdown.gene]

Genedependency <- read.csv("/WData/WData2/xiangjun/ChangZeWen/Data/CRISPRGeneDependency.csv", check.names = FALSE, row.names = 1)
colnames(Genedependency) <- gsub(" *\\(.*\\)$", "", colnames(Genedependency))
Genedependency <- Genedependency[, colnames(Genedependency) %in% Knockdown.gene]


### Part --[2]----
#' gene effect.

plot_df <- Geneeffect %>%
    tibble::rownames_to_column("CellLine") %>%
    pivot_longer(
        cols = -CellLine,
        names_to = "Gene",
        values_to = "Score"
    )

#' 按中位数排序.
plot_df$Gene <- factor(
    plot_df$Gene,
    levels = plot_df %>%
        group_by(Gene) %>%
        summarise(med = median(Score, na.rm = TRUE)) %>%
        arrange(med) %>%
        pull(Gene)
)

#' 画图.
pdf("./Figure/8.Gene.Knockdown/2.3.DepMap.Geneeffect图.pdf", width = 8, height = 5)
ggplot(plot_df, aes(x = Gene, y = Score, fill = Gene)) +
    geom_boxplot(width = 0.65, outlier.size = 1.8, outlier.alpha = 0.8) +
    geom_hline(yintercept = -0.5, linetype = "dashed", color = "red", linewidth = 1) +
    theme_classic(base_size = 14) +
    theme(
        axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"),
        axis.text.y = element_text(face = "bold"),
        axis.title = element_text(face = "bold"),
        legend.position = "none"
    ) +
    labs(x = NULL, y = "CRISPR Gene Effect")
dev.off()


### Part --[3]----
#' gene dependency

#' 宽表转长表.
plot_df1 <- Genedependency %>%
    tibble::rownames_to_column("CellLine") %>%
    pivot_longer(
        cols = -CellLine,
        names_to = "Gene",
        values_to = "Score"
    )

#'
pdf("./Figure/8.Gene.Knockdown/2.4.DepMap.Genedependency图.pdf", width = 8, height = 5)
ggplot(plot_df1, aes(x = Gene, y = Score, fill = Gene)) +
    geom_boxplot(width = 0.65, outlier.size = 1.8, outlier.alpha = 0.8) +
    geom_hline(yintercept = -0.5, linetype = "dashed", color = "red", linewidth = 1) +
    theme_classic(base_size = 14) +
    theme(
        axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"),
        axis.text.y = element_text(face = "bold"),
        axis.title = element_text(face = "bold"),
        legend.position = "none"
    ) +
    labs(x = NULL, y = "CRISPR Gene Effect")
dev.off()
