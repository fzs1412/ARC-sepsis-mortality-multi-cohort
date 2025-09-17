
library(dplyr)
library(tibble)
library(clusterProfiler)
library(org.Hs.eg.db)
library(fgsea) # 快速GSEA分析包
library(ggplot2)
library(Cairo)
library(msigdbr) 
library(data.table)   # collapsePathways/plotGseaTable 依赖 data.table 语法
options(expressions = 5e5) 

# **【网络修复】**: 延长网络下载超时时间至300秒 (5分钟)
# 这可以解决因网络慢导致msigdbr数据下载失败的问题
options(timeout = 300)


# --- 步骤 2: 准备用于GSEA的基因排序列表 ---

# 1. 加载您之前保存的DESeq2结果对象
#    确保 "DESeq2_results_object.rds" 文件在您的工作目录中
if (!file.exists("DESeq2_results_object.rds")) {
  stop("错误: 找不到已保存的 'DESeq2_results_object.rds' 文件。请先运行 `normalize_omics_data.R` 脚本。")
}
res <- readRDS("DESeq2_results_object.rds")

# 2. 准备数据框，并移除NA值
res_df <- as.data.frame(res) %>%
  rownames_to_column("ENSEMBL") %>%
  filter(!is.na(stat)) # GSEA使用 'stat' 列进行排序，它比log2FC更稳健

# 3. 将ENSEMBL ID转换为ENTREZ ID，这是GSEA进行通路匹配需要的
ensembl_to_entrez <- bitr(res_df$ENSEMBL,
                          fromType = "ENSEMBL",
                          toType = "ENTREZID",
                          OrgDb = "org.Hs.eg.db")

# 4. 将ENTREZ ID合并回结果，并创建最终的排序列表
#    这个列表是一个以ENTREZ ID为名，以 'stat' 值为内容的命名向量
ranks <- res_df %>%
  inner_join(ensembl_to_entrez, by = "ENSEMBL") %>%
  # 处理一个ENTREZ ID对应多个ENSEMBL ID的情况，保留 'stat' 值绝对值最大的那个
  group_by(ENTREZID) %>%
  slice_max(order_by = abs(stat), n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  { setNames(.$stat, .$ENTREZID) } %>%
  sort(decreasing = TRUE) # 降序排列

cat("已成功创建包含", length(ranks), "个基因的排序列表用于GSEA分析。\n")


# --- 步骤 3: 【核心修复】准备通路基因集 ---
msigdbr_go_bp <- msigdbr(
  species = "Homo sapiens",
  category = "C5",
  subcategory = "GO:BP"
)

# 关键：把 entrez_gene 转为字符；否则与 ranks 名称（字符）不匹配
msig_t2g <- msigdbr_go_bp %>%
  dplyr::select(gs_name, entrez_gene) %>%
  dplyr::distinct() %>%
  dplyr::mutate(ENTREZID = as.character(entrez_gene)) %>%
  dplyr::select(gs_name, ENTREZID)

go_pathways_entrez <- split(msig_t2g$ENTREZID, msig_t2g$gs_name)

cat("已加载 ", length(go_pathways_entrez), " 个GO:BP通路；覆盖 ",
    length(unique(unlist(go_pathways_entrez))), " 个ENTREZ基因（字符型）。\n", sep = "")

# --- 步骤 4: 运行 GSEA 分析 ---
fgsea_res <- fgsea(
  pathways = go_pathways_entrez,
  stats = ranks,
  minSize = 15,
  maxSize = 500,
  nPermSimple = 10000
)

fgsea_res_sorted <- fgsea_res %>%
  as_tibble() %>%
  arrange(desc(NES))

print("--- GSEA分析中富集得分最高的通路 ---")
print(head(fgsea_res_sorted %>% select(-leadingEdge, -ES), 10))

top_pathway_up <- fgsea_res_sorted %>%
  filter(NES > 0) %>%
  slice_min(padj, n = 1, with_ties = FALSE) %>%
  pull(pathway)

top_pathway_down <- fgsea_res_sorted %>%
  filter(NES < 0) %>%
  slice_min(padj, n = 1, with_ties = FALSE) %>%
  pull(pathway)

p_enrich_up <- if (length(top_pathway_up) > 0)
  plotEnrichment(go_pathways_entrez[[top_pathway_up]], ranks) + labs(title = top_pathway_up) else NULL
p_enrich_down <- if (length(top_pathway_down) > 0)
  plotEnrichment(go_pathways_entrez[[top_pathway_down]], ranks) + labs(title = top_pathway_down) else NULL

# --- 步骤5（建议放到库加载处）---

# ======= 最简输出：Figure 6（两联图）+ 两个补充表（PDF） =======

# 需要一个轻量依赖：gridExtra 用来拼图 & 画表格
if (!requireNamespace("gridExtra", quietly = TRUE)) install.packages("gridExtra")
library(gridExtra)
library(grid)

# 1) 幂等清理：基因集与 ranks 类型统一
go_pathways_entrez <- lapply(go_pathways_entrez, function(x) unique(as.character(x)))
ranks_ok <- ranks[is.finite(ranks)]
names(ranks_ok) <- as.character(names(ranks_ok))
ranks_ok <- ranks_ok[!duplicated(names(ranks_ok))]
ranks_ok <- sort(ranks_ok, decreasing = TRUE)

# 2) 选“最显著上/下调”通路（按 FDR/padj 最小）
fg_df <- as.data.frame(fgsea_res)
fg_df <- fg_df[!is.na(fg_df$padj), ]

# 最显著上调
stopifnot(any(fg_df$NES > 0))
up_df <- fg_df[fg_df$NES > 0, ]
top_up <- up_df[which.min(up_df$padj), ]
pw_up  <- top_up$pathway
gs_up  <- go_pathways_entrez[[pw_up]]

# 最显著下调（若没有则只画上调）
has_down <- any(fg_df$NES < 0)
if (has_down) {
  dn_df <- fg_df[fg_df$NES < 0, ]
  top_dn <- dn_df[which.min(dn_df$padj), ]
  pw_dn  <- top_dn$pathway
  gs_dn  <- go_pathways_entrez[[pw_dn]]
}

# 3) 画 enrichment 曲线（A/B）
p_up <- fgsea::plotEnrichment(gs_up, ranks_ok) +
  ggplot2::labs(
    title    = paste0("A. ", pw_up),
    subtitle = sprintf("NES = %.2f | FDR = %.2g | size = %d",
                       top_up$NES, top_up$padj, top_up$size)
  ) +
  theme_minimal(base_size = 12) +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major.y = element_blank(),
        axis.title = element_blank(),
        plot.title = element_text(face = "bold"))

if (has_down) {
  p_dn <- fgsea::plotEnrichment(gs_dn, ranks_ok) +
    ggplot2::labs(
      title    = paste0("B. ", pw_dn),
      subtitle = sprintf("NES = %.2f | FDR = %.2g | size = %d",
                         top_dn$NES, top_dn$padj, top_dn$size)
    ) +
    theme_minimal(base_size = 12) +
    theme(panel.grid.minor = element_blank(),
          panel.grid.major.y = element_blank(),
          axis.title = element_blank(),
          plot.title = element_text(face = "bold"))
}

# 4) 导出 Figure 6（PDF，竖排两联图；若无下调则仅一页一图）
dir.create("../figure", showWarnings = FALSE, recursive = TRUE)
Cairo::CairoPDF("../figure/Figure6.pdf", width = 8, height = if (has_down) 10 else 5)
if (has_down) {
  grid.arrange(p_up, p_dn, ncol = 1)
} else {
  grid.arrange(p_up)
}
dev.off()

# ===== 补充表合并为一个 Word（Top10 上/下调）=====
if (!requireNamespace("officer", quietly = TRUE)) install.packages("officer")
if (!requireNamespace("flextable", quietly = TRUE)) install.packages("flextable")
library(officer); library(flextable)

dir.create("../table", showWarnings = FALSE, recursive = TRUE)

fg_df <- as.data.frame(fgsea_res)
fg_df <- fg_df[!is.na(fg_df$padj), c("pathway","NES","pval","padj","size")]

up10 <- subset(fg_df, NES > 0)
up10 <- head(up10[order(up10$padj), ], 10)

dn10 <- subset(fg_df, NES < 0)
dn10 <- head(dn10[order(dn10$padj), ], 10)

for (v in c("NES","pval","padj")) {
  if (nrow(up10) > 0) up10[[v]] <- signif(up10[[v]], 3)
  if (nrow(dn10) > 0) dn10[[v]] <- signif(dn10[[v]], 3)
}

doc <- read_docx()
doc <- body_add_par(doc, "补充表：GSEA Top10（上调/下调）", style = "heading 1")
if (nrow(up10) > 0) {
  doc <- body_add_par(doc, "上调 Top10 通路", style = "heading 2")
  doc <- body_add_flextable(doc, autofit(flextable(up10)))
}
if (nrow(dn10) > 0) {
  doc <- body_add_par(doc, "下调 Top10 通路", style = "heading 2")
  doc <- body_add_flextable(doc, autofit(flextable(dn10)))
}

print(doc, target = "../table/Supp_Table_GSEA_Top10.docx")
cat("已生成：../table/Supp_Table_GSEA_Top10.docx\n")

