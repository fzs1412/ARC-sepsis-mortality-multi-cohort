# --- 依赖包 ---
library(data.table)
library(dplyr)
library(stringr)
library(tibble)
library(DESeq2)
library(clusterProfiler)
library(org.Hs.eg.db)
library(fgsea)        # 快速GSEA分析包（此脚本未直接用到，但保留）
library(ggplot2)
library(Cairo)
library(EnhancedVolcano)
library(pheatmap)
library(patchwork)
library(grid)

# --- Step 2: 加载并准备表达矩阵 ---
setwd('d:/OneDrive/R/2024Sepsis_ARC/DataAnalysis')

# 来自同一批次的计数文件
count_lc <- fread("D:/OneDrive/R/01sepsis/OMIX006457/OMIX006457-02.csv", data.table = FALSE)
count_singlorn <- fread("D:/OneDrive/R/01sepsis/OMIX006457/OMIX006457-07.csv", data.table = FALSE)

# 将第一列设为行名
rownames(count_lc) <- count_lc[, 1];      count_lc <- count_lc[, -1]
rownames(count_singlorn) <- count_singlorn[, 1]; count_singlorn <- count_singlorn[, -1]

# --- Step 3: 标准化基因ID（同时作用于两份矩阵）---
# 从行名中提取 ENSG
rownames(count_lc)       <- str_split_i(rownames(count_lc),       "\\|", 1)
rownames(count_singlorn) <- str_split_i(rownames(count_singlorn), "\\|", 1)
cat("已将两份计数矩阵的基因名标准化为纯 ENSG。\n")

# --- Step 4: 合并原始计数矩阵 ---
common_genes <- intersect(rownames(count_lc), rownames(count_singlorn))
cat("共有的基因数量: ", length(common_genes), "\n")

count_lc_filtered       <- count_lc[common_genes, ]
count_singlorn_filtered <- count_singlorn[common_genes, ]

combined_counts <- cbind(count_lc_filtered, count_singlorn_filtered)
cat("合并后的原始计数矩阵维度: ", nrow(combined_counts), "个基因, ", ncol(combined_counts), "个样本\n")

# 计数为整数
combined_counts <- round(combined_counts)

# (可选) 保存合并后的原始计数矩阵
fwrite(as.data.frame(combined_counts) %>% rownames_to_column("GeneID"),
       "D:/OneDrive/R/01sepsis/OMIX006457/combined_raw_counts.csv")
cat("合并后的原始计数矩阵已保存。\n")

# --- Step 5: 加载并准备临床数据 (coldata) ---
coldata_raw <- fread("D:/OneDrive/R/2024Sepsis_ARC/DataAnalysis/CMAISEday1.csv") %>%
  select(-no) %>%
  dplyr::rename(class = ARC) %>%
  mutate(class = factor(class, levels = c("0", "1"))) %>%
  column_to_rownames("SampleName")

# --- Step 6: 匹配表达矩阵和临床数据 ---
# 筛选 Day 1 的样本（列名不以 d3/d5 结尾）
day1_samples <- colnames(combined_counts)[!grepl("d3$|d5$", colnames(combined_counts))]
countdata_d1 <- combined_counts[, day1_samples]
cat("筛选出 Day 1 的样本数: ", ncol(countdata_d1), "\n")

common_samples <- intersect(colnames(countdata_d1), rownames(coldata_raw))
cat("临床数据与表达谱共有的样本数: ", length(common_samples), "\n")

countdata_final <- countdata_d1[, common_samples]
coldata_final   <- coldata_raw[common_samples, ]

stopifnot(all(colnames(countdata_final) == rownames(coldata_final)))

# --- Step 7: 协变量校正后的 DESeq2 差异表达 ---
# 移除协变量 NA
needed_cov <- c("age", "sex", "SOFA", "lac", "class")
complete_covariates <- complete.cases(coldata_final[, needed_cov])
coldata_final <- coldata_final[complete_covariates, ]
countdata_final <- countdata_final[, rownames(coldata_final)]
cat("移除协变量NA后，最终样本数: ", nrow(coldata_final), "\n")

dds <- DESeqDataSetFromMatrix(countData = countdata_final,
                              colData   = coldata_final,
                              design    = ~ age + sex + SOFA + lac + class)

# 过滤低表达基因
keep <- rowSums(counts(dds)) >= 10
dds  <- dds[keep,]
cat("过滤低表达基因后剩余: ", nrow(dds), "个基因\n")

dds_analyzed <- DESeq(dds)
cat("DESeq2分析完成。\n")

res <- results(dds_analyzed, contrast = c("class", "1", "0"))
summary(res)

# 保存结果（保留保存；不再立即读回）
fwrite(as.data.frame(res) %>% rownames_to_column("GeneID"),
       "DESeq2_results_ARC_vs_NoARC_Adjusted.csv")
saveRDS(res,          file = "DESeq2_results_object.rds")
saveRDS(dds_analyzed, file = "dds_analyzed_object.rds")

# 与 dds 保持一致的 coldata（用于作图注释）
coldata_final <- as.data.frame(colData(dds_analyzed))

# --- 步骤 8: 创建并清理唯一的基因ID-SYMBOL映射表（去重后仅保留一次） ---
# 用原始文件中最完整的 "ENSEMBL|SYMBOL" 信息构建映射
count_lc_map_df     <- fread("D:/OneDrive/R/01sepsis/OMIX006457/OMIX006457-02.csv", data.table = FALSE)
original_gene_names <- count_lc_map_df[, 1]

initial_gene_map <- data.frame(
  ENSEMBL = str_split_i(original_gene_names, "\\|", 1),
  SYMBOL  = str_split_i(original_gene_names, "\\|", 2),
  stringsAsFactors = FALSE
) %>% filter(!is.na(ENSEMBL) & ENSEMBL != "")

# 合并映射并解决多对一（每个SYMBOL保留 padj 最小的 ENSEMBL）
res_df_for_map <- as.data.frame(res) %>%
  rownames_to_column("ENSEMBL") %>%
  left_join(initial_gene_map, by = "ENSEMBL")

gene_map_cleaned <- res_df_for_map %>%
  mutate(SYMBOL = ifelse(is.na(SYMBOL) | SYMBOL == "", ENSEMBL, SYMBOL)) %>%
  group_by(SYMBOL) %>%
  slice_min(order_by = padj, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  select(ENSEMBL, SYMBOL)

gene_map_cleaned$SYMBOL <- make.unique(gene_map_cleaned$SYMBOL)
cat("创建了包含", nrow(gene_map_cleaned), "个基因的唯一 ID-SYMBOL 映射表。\n")

# --- 步骤 9: 火山图 (Figure 5A) ---
 
res_df_for_volcano <- as.data.frame(res) %>%
  rownames_to_column("ENSEMBL") %>%
  left_join(gene_map_cleaned, by = "ENSEMBL")

genes_to_label <- res_df_for_volcano %>%
  filter(padj < 0.05 & abs(log2FoldChange) > 1) %>%
  arrange(padj) %>%
  head(20) %>%
  pull(SYMBOL)

p_volcano <- EnhancedVolcano(
  res_df_for_volcano,
  lab = res_df_for_volcano$SYMBOL,
  selectLab = genes_to_label,
  x = 'log2FoldChange',
  y = 'padj',
  xlab = expression(Log[2]~'Fold Change'),
  ylab = expression(-Log[10]~'Adjusted P-value'),
  axisLabSize = 14,
  pCutoff = 0.05,
  FCcutoff = 1.0,
  pointSize = 2.0,
  labSize = 4.0,
  labFace = 'bold',
  title = "A. Volcano plot of DEGs (ARC vs Non-ARC)",
  subtitle = "Adjusted for age, sex, SOFA score, and lactate",
  col = c('grey30', '#2166AC', '#2166AC', '#B2182B'),
  colAlpha = 0.8,
  legendPosition = "bottom",
  legendLabels = c('NS', 'Log2 FC', 'Adjusted P', 'Adjusted P & Log2 FC'),
  legendLabSize = 12,
  legendIconSize = 5.0,
  drawConnectors = TRUE,
  widthConnectors = 0.5,
  colConnectors = 'grey50',
  gridlines.major = FALSE,
  gridlines.minor = FALSE,
  border = 'full',
  borderWidth = 1.5,
  borderColour = 'black'
) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5))


# --- 步骤 10: 热图 (Figure 5B) ---
# --- Step 10: Heatmap (Figure 5B, English) ---
# vst transform
vsd <- vst(dds_analyzed, blind = FALSE)
vst_matrix <- assay(vsd)

# Select DEGs for heatmap
degs_for_heatmap <- res_df_for_volcano %>%
  filter(padj < 0.05 & abs(log2FoldChange) > 1.5)

cat("Number of significant genes for heatmap: ", nrow(degs_for_heatmap), "\n")

heatmap_matrix <- vst_matrix[degs_for_heatmap$ENSEMBL, ]
rownames(heatmap_matrix) <- degs_for_heatmap$SYMBOL

# Column annotation in English
annotation_col <- data.frame(
  Group = factor(coldata_final$class, levels = c("0", "1"),
                 labels = c("Non-ARC", "ARC")),
  row.names = colnames(heatmap_matrix)
)

pheatmap_plot <- pheatmap(
  heatmap_matrix,
  scale = "row",
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  show_colnames = FALSE,
  show_rownames = TRUE,
  fontsize_row = 8,
  annotation_col = annotation_col,
  annotation_colors = list(
    Group = c("Non-ARC" = "#2166AC", "ARC" = "#B2182B")
  ),
  border_color = 'white',
  main = "B. Heatmap of Differentially Expressed Genes",
  color = colorRampPalette(c("navy", "white", "firebrick3"))(50)
)

# --- 步骤 11: 组合并保存 ---
p_heatmap_wrapped <- wrap_elements(pheatmap_plot$gtable)

final_plot <- p_volcano + p_heatmap_wrapped + plot_layout(widths = c(1, 1.2))

ggsave(
  "../figure/Figure5_Volcano_Heatmap_Combined.pdf",
  plot = final_plot,
  width = 16, height = 9, units = "in",
  device = cairo_pdf
)
ggsave(
  "../figure/Figure5_Volcano_Heatmap_Combined.tiff",
  plot = final_plot,
  width = 16, height = 9, units = "in",
  dpi = 300, compression = "lzw"
)

cat("Figure 5 (火山图与热图) 已成功保存为 PDF 和 TIFF。\n")
