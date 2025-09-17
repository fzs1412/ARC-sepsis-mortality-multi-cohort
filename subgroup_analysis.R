# --- Step 1: Load Libraries and Prepare Combined Data ---
# 确保所有必要的包都已安装和加载
# install.packages(c("data.table", "dplyr", "meta", "Cairo", "forestploter", "lmtest"))
library(data.table)
library(dplyr)
library(meta)
library(Cairo)
library(forestploter)
library(grid)
library(lmtest) # 用于更稳健地计算交互作用P值

# 加载您的主数据集
dat <- fread("d:/OneDrive/R/2024Sepsis_ARC/DataAnalysis/RICU_sepsis_com3.csv")

# **关键步骤: 创建一个用于分析的完整、合并的数据集**

# 1. 为CMAISE创建代理charlson评分
dat_CMAISE <- dat[dat$source == "CMAISE1.5v", ] %>%
  mutate(
    charlson_proxy = Hypertension + Diabete + CHF + MI + COPD,
    charlson = charlson_proxy
  )
# 将CMAISE与其他数据合并
dat_other <- dat[dat$source != "CMAISE1.5v", ]
dat_combined <- rbind(dat_other, dat_CMAISE, fill = TRUE)

# 2. 创建亚组变量
dat_combined <- dat_combined %>%
  mutate(
    age_group = ifelse(age < 65, "< 65 years", ">= 65 years"),
    sofa_group = ifelse(sofa < 8, "< 8", ">= 8"),
    sex = factor(sex, levels = c(0, 1), labels = c("Female", "Male")),
    source_group = factor(case_when(
      source == "miiv" ~ "MIMIC-IV",
      source == "mimic" ~ "MIMIC-III",
      source == "eicu" ~ "eICU",
      source == "aumc" ~ "AUMC",
      source == "CMAISE1.5v" ~ "CMAISE",
      TRUE ~ as.character(source)
    ))
  )

# --- Step 2: 创建一个用于亚组分析的函数 ---

run_subgroup_analysis <- function(subgroup_var, subgroup_label, data) {
  
  results <- list()
  
  full_covariates <- "age + sex + charlson + vaso_day1 + sofa + wbc + lact + ph + bili + crea + pf + hr_max + temp_min"
  reduced_covariates <- "age + sex + sofa + lact" 
  
  subgroup_levels <- levels(factor(data[[subgroup_var]]))
  for (level in subgroup_levels) {
    sub_data <- data[data[[subgroup_var]] == level, ]
    
    if(nrow(sub_data) < 100 || length(unique(sub_data$ARC)) < 2 || sum(sub_data$death, na.rm=TRUE) < 10) next
    
    num_events <- sum(sub_data$death, na.rm = TRUE)
    base_covariates <- if (num_events < 50) reduced_covariates else full_covariates
    
    current_covariates <- base_covariates
    if (length(unique(na.omit(sub_data$sex))) < 2) {
      current_covariates <- gsub("\\+ sex|sex \\+", "", current_covariates)
    }
    
    formula_sub <- as.formula(paste("death ~ ARC +", current_covariates))
    model_sub <- glm(formula_sub, data = sub_data, family = "binomial")
    
    if ("ARC" %in% names(coef(model_sub))) {
      or <- exp(coef(model_sub)["ARC"])
      ci <- exp(confint.default(model_sub)["ARC", ])
      
      results[[as.character(level)]] <- data.frame(
        Subgroup = paste("  ", level),
        N = nrow(sub_data),
        Events = num_events,
        OR = or,
        Lower = ci[1],
        Upper = ci[2]
      )
    }
  }
  
  interaction_covariates <- full_covariates
  if (length(unique(na.omit(data$sex))) < 2) {
    interaction_covariates <- gsub("\\+ sex|sex \\+", "", interaction_covariates)
  }
  
  formula_interaction <- as.formula(paste("death ~ ARC *", subgroup_var, "+", interaction_covariates))
  model_interaction <- glm(formula_interaction, data = data, family = "binomial")
  
  formula_main <- as.formula(paste("death ~ ARC +", subgroup_var, "+", interaction_covariates))
  model_main <- glm(formula_main, data=data, family="binomial")
  
  lrt <- lrtest(model_main, model_interaction)
  p_interaction_value <- lrt$`Pr(>Chisq)`[2]
  
  final_df <- do.call(rbind, results)
  if(is.null(final_df)) return(NULL)
  
  header_row <- data.frame(Subgroup = subgroup_label, N = NA, Events = NA, OR = NA, Lower = NA, Upper = NA, P_interaction_text = sprintf("%.3f", p_interaction_value))
  final_df$P_interaction_text <- ""
  
  return(rbind(header_row, final_df))
}


# --- Step 3: 运行所有分析 ---

run_binary_subgroup <- function(var_name, label, data) {
  data_temp <- data %>% mutate(sub_var = factor(ifelse(!!sym(var_name) == 1, "Yes", "No")))
  return(run_subgroup_analysis("sub_var", label, data_temp))
}

subgroup_results <- rbind(
  run_subgroup_analysis("source_group", "Database", dat_combined),
  run_subgroup_analysis("age_group", "Age Group", dat_combined),
  run_subgroup_analysis("sex", "Sex", dat_combined),
  run_subgroup_analysis("sofa_group", "SOFA Score", dat_combined),
  run_binary_subgroup("Hypertension", "Hypertension", dat_combined),
  run_binary_subgroup("Diabete", "Diabetes", dat_combined)
)

overall_covariates <- "age + sex + charlson + vaso_day1 + sofa + wbc + lact + ph + bili + crea + pf + hr_max + temp_min"
if (length(unique(na.omit(dat_combined$sex))) < 2) {
  overall_covariates <- gsub("\\+ sex|sex \\+", "", overall_covariates)
}

full_model <- glm(as.formula(paste("death ~ ARC +", overall_covariates)), data = dat_combined, family = "binomial")
overall_or <- exp(coef(full_model)["ARC"])
overall_ci <- exp(confint.default(full_model)["ARC", ])
overall_df <- data.frame(
  Subgroup = "Overall", 
  N = nrow(dat_combined), 
  Events = sum(dat_combined$death, na.rm = TRUE), 
  OR = overall_or, 
  Lower = overall_ci[1], 
  Upper = overall_ci[2],
  P_interaction_text = ""
)

final_table <- rbind(overall_df, subgroup_results)


# --- Step 4: 格式化表格并生成森林图 ---

# **【解决标题问题】**: 创建一个有正确列名的数据框用于绘图
plot_data <- final_table %>%
  mutate(
    N = ifelse(!is.na(N), format(N, big.mark = ","), ""),
    Events = ifelse(!is.na(Events), as.character(Events), ""),
    # **【解决重叠】**: 将文本和图形分开
    `OR (95% CI)` = ifelse(
      !is.na(OR),
      sprintf("%.2f (%.2f-%.2f)", OR, Lower, Upper),
      ""
    ),
    # 空白列用于容纳森林图的图形部分
    ` ` = paste(rep(" ", 20), collapse = ""),
    # **【解决标题问题】**: 这是最后一列，将自动获得标题
    `P for Interaction` = P_interaction_text,
    # 标记标题行以便加粗
    is_header = !grepl("  ", Subgroup) & Subgroup != "Overall"
  )

# **【最终美化】: 定义一个更专业、更宽敞的森林图主题**
tm <- forest_theme(
  base_size = 10,
  # **【解决重叠】**: 增加行高
  row_height = unit(0.9, "cm"),
  ci_pch = 16, ci_col = "#377eb8", ci_lty = 1, ci_lwd = 1.5,
  refline_lwd = 1, refline_lty = "dashed", refline_col = "grey20",
  summary_pch = 18, summary_col = "#e41a1c",
  footnote_cex = 0.8,
  # 对齐列标题
  colhead_just = c("left", "right", "right", "right", "center", "right")
)

# 绘制森林图
p <- forest(
  data = plot_data[, c("Subgroup", "N", "Events", "OR (95% CI)", " ", "P for Interaction")],
  est = plot_data$OR,
  lower = plot_data$Lower,
  upper = plot_data$Upper,
  sizes = 0.8,
  is_summary = plot_data$Subgroup == "Overall",
  # **【解决重叠】**: ci_column现在指向空白列
  ci_column = 5, 
  ref_line = 1,
  arrow_lab = c("ARC Protective", "ARC Harmful"),
  xlim = c(0.4, 1.6),
  # **【解决重叠】**: 精确设置列宽
  widths = unit(c(3.5, 1.5, 1.5, 2.5, 3.5, 1.5), "cm"), # Subgroup, N, Events, OR(text), Plot, P-int
  theme = tm
)

# 为标题行加粗
p <- edit_plot(p, row = which(plot_data$is_header), gp = gpar(fontface = "bold"))

# --- Step 5: 保存图表 ---
ggsave(
  "../figure/Figure4_Subgroup_Analysis.pdf",
  plot = p, device = "pdf", width = 12, height = 11, units = "in" # 增加尺寸
)

ggsave(
  "../figure/Figure4_Subgroup_Analysis.tiff",
  plot = p, device = "tiff", width = 12, height = 11, units = "in", dpi = 300, compression = "lzw"
)

print("Figure 4 (Subgroup Analysis Forest Plot) has been saved as PDF and TIFF.")

