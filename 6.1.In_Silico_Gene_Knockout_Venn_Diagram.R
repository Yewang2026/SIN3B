### Part --[1]----

setwd("/WData/WorkSpace/xiangjun2024/LiJunHu")
library(ggvenn) # 读入包
library(UpSetR)
library(openxlsx)
library(ggVennDiagram)
library(dplyr)
library(ggplot2)

cluster1 <- read.csv("./Data/Geneformer/geneformer_significant_affected_genes/cluster1_all_significant_affected_genes.csv")
cluster2 <- read.csv("./Data/Geneformer/geneformer_significant_affected_genes/cluster2_all_significant_affected_genes.csv")
cluster3 <- read.csv("./Data/Geneformer/geneformer_significant_affected_genes/cluster3_all_significant_affected_genes.csv")
cluster4 <- read.csv("./Data/Geneformer/geneformer_significant_affected_genes/cluster4_all_significant_affected_genes.csv")
cluster5 <- read.csv("./Data/Geneformer/geneformer_significant_affected_genes/cluster5_all_significant_affected_genes.csv")

cluster1 <- unique(cluster1$Affected_gene)
cluster2 <- unique(cluster2$Affected_gene)
cluster3 <- unique(cluster3$Affected_gene)
cluster4 <- unique(cluster4$Affected_gene)
cluster5 <- unique(cluster5$Affected_gene)

### Part --[2]----

VENN.LIST <- list(
    cluster1 = as.character(cluster1),
    cluster2 = as.character(cluster2),
    cluster3 = as.character(cluster3),
    cluster4 = as.character(cluster4),
    cluster5 = as.character(cluster5)
)


pdf("./Figure/8.Gene.Knockdown/2.1.敲除基因韦恩图.pdf", width = 8, height = 8)
ggvenn(VENN.LIST, c( # "All.Samples",
    "cluster1", "cluster2", "cluster3", "cluster4", "cluster5"
),
fill_color = c("#C3D0E3", "#A4E1E8", "#79D5C3", "#F9CBAA", "#F2B5B5"), # 填充颜色
fill_alpha = 0.5, # 改变图形透明度
stroke_color = "black", # 边界线的颜色
stroke_alpha = 0.5, # 边界线的透明度
stroke_linetype = "dashed", # 将边界线转化为虚线
text_size = 4
) # 调节字体大小
dev.off()


### Part --[3]----


pdf("./Figure/8.Gene.Knockdown/2.2.敲除基因upset图.pdf", width = 6, height = 5)
upset(fromList(VENN.LIST), order.by = "freq")
dev.off()


### Part --[4]----


x <- VENN.LIST
names(x)
id <- unique(unlist(x))
id
y <- do.call(cbind, lapply(x, function(i) {
    sign(id %in% i)
}))
rownames(y) <- id
y <- data.frame(y)
y$sum <- rowSums(y)
Knockdown.gene <- y[order(y$sum, decreasing = T), ]
Knockdown.gene <- tibble::rownames_to_column(Knockdown.gene, "GeneSymbol")
write.xlsx(Knockdown.gene, file = "./Data/Geneformer/Knockdown.gene.xlsx")
