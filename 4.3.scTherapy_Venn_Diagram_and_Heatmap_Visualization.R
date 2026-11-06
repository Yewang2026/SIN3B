### Part --[1]----
setwd("/WData/WorkSpace/xiangjun2024/LiJunHu")
library(ggvenn)
library(UpSetR)
library(openxlsx)
library(ggVennDiagram)
library(dplyr)
library(ggplot2)

scTherapy.drugs <- readRDS("./Data/scTherapy/scTherapy.drugs.rds")

cluster1 <- scTherapy.drugs[scTherapy.drugs$Subtype == "cluster1", ]
cluster1 <- cluster1$Drug_Name
cluster2 <- scTherapy.drugs[scTherapy.drugs$Subtype == "cluster2", ]
cluster2 <- cluster2$Drug_Name
cluster3 <- scTherapy.drugs[scTherapy.drugs$Subtype == "cluster3", ]
cluster3 <- cluster3$Drug_Name
cluster4 <- scTherapy.drugs[scTherapy.drugs$Subtype == "cluster4", ]
cluster4 <- cluster4$Drug_Name
cluster5 <- scTherapy.drugs[scTherapy.drugs$Subtype == "cluster5", ]
cluster5 <- cluster5$Drug_Name


### Part --[2]----
VENN.LIST <- list(
  cluster1 = as.character(cluster1),
  cluster2 = as.character(cluster2),
  cluster3 = as.character(cluster3),
  cluster4 = as.character(cluster4),
  cluster5 = as.character(cluster5)
)

pdf("./Figure/6.Bulk分析/4.1.scTherapy韦恩图.pdf", width = 8, height = 8)
ggvenn(VENN.LIST, c( # "All.Samples",
  "cluster1", "cluster2", "cluster3", "cluster4", "cluster5"
),
fill_color = c("#C3D0E3", "#A4E1E8", "#79D5C3", "#F9CBAA", "#F2B5B5"), # 填充颜色
fill_alpha = 0.5, #
stroke_color = "black", #
stroke_alpha = 0.5, #
stroke_linetype = "dashed", #
text_size = 4
)
dev.off()


### Part --[3]----
pdf("./Figure/6.Bulk分析/4.2.scTherapy.upset图.pdf", width = 6, height = 5)
upset(fromList(VENN.LIST), order.by = "freq")
dev.off()


### Part --[4]----

drug.df <- scTherapy.drugs %>%
  dplyr::select(Drug_Name, Subtype) %>%
  distinct()

drug.df$Subtype <- factor(drug.df$Subtype, levels = c("cluster1", "cluster2", "cluster3", "cluster4", "cluster5"))

drug.stat <- dplyr::summarise(
  dplyr::group_by(drug.df, Drug_Name),
  nSubtype = dplyr::n_distinct(Subtype),
  subtype_set = paste(sort(unique(as.character(Subtype))), collapse = "_"),
  .groups = "drop"
)

drug.plot <- drug.df %>%
  left_join(drug.stat, by = "Drug_Name")


drug.order <- drug.plot %>%
  arrange(desc(nSubtype), subtype_set, Drug_Name) %>%
  pull(Drug_Name) %>%
  unique()

drug.plot$Drug_Name <- factor(drug.plot$Drug_Name, levels = rev(drug.order))

col.use <- c("#C3D0E3", "#A4E1E8", "#79D5C3", "#F9CBAA", "#F2B5B5")
names(col.use) <- c("cluster1", "cluster2", "cluster3", "cluster4", "cluster5")

pdf("./Figure/6.Bulk分析/4.3.scTherapy共享热图.pdf", width = 5, height = 10)
ggplot(drug.plot, aes(x = Subtype, y = Drug_Name, color = Subtype)) +
  geom_point(size = 4.5) +
  scale_color_manual(values = col.use) +
  scale_x_discrete(labels = c(
    "cluster1" = "cluster1",
    "cluster2" = "cluster2",
    "cluster3" = "cluster3",
    "cluster4" = "cluster4",
    "cluster5" = "cluster5"
  )) +
  labs(x = NULL, y = "Drug") +
  theme_bw() +
  theme(
    panel.grid.major = element_line(color = "#D9D9D9", linewidth = 0.4),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 12, face = "bold"),
    axis.text.y = element_text(size = 10, face = "bold"),
    axis.title.y = element_text(size = 14, face = "bold"),
    legend.position = "none",
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8)
  )
dev.off()
