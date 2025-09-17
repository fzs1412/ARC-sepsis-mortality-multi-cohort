# --- Step 1: Load Libraries and Prepare Data ---
# 确保所有必要的包都已安装和加载
# install.packages(c("data.table", "dplyr", "survival", "cmprsk", "forestploter", "grid", "Cairo"))
library(data.table)
library(dplyr)
library(survival)
library(cmprsk)
library(forestploter)
library(grid)
library(Cairo)

# 加载您的主数据集
dat_original <- fread("d:/OneDrive/R/2024Sepsis_ARC/DataAnalysis/RICU_sepsis_com3.csv")

# --- Step 2: Prepare Data for Analysis ---

# 1. 为CMAISE创建代理charlson评分
dat_CMAISE <- dat_original[dat_original$source == "CMAISE1.5v", ] %>%
  mutate(
    charlson_proxy = Hypertension + Diabete + CHF + MI + COPD,
    charlson = charlson_proxy
  )
# 将CMAISE与其他数据合并
dat_other <- dat_original[dat_original$source != "CMAISE1.5v", ]
dat_combined <- rbind(dat_other, dat_CMAISE, fill = TRUE)

# 2. 创建多状态结局变量
dat_combined <- dat_combined %>%
  mutate(
    event_status_cr = factor(ifelse(death == 1, 1, 2)),
    time_to_event = los_hosp
  )

# 拆分数据
dat_list <- split(dat_combined, dat_combined$source)
# Reorder for desired plot appearance
desired_order <- c("MIMIC-IV", "MIMIC-III", "eICU", "AUMC", "CMAISE")
names(dat_list) <- c("AUMC", "CMAISE", "eICU", "MIMIC-III", "MIMIC-IV")
dat_list <- dat_list[desired_order]


# --- Step 3: 创建一个函数来运行两种模型并返回适合绘图的数据 ---

run_both_models <- function(df, study_name) {
  
  full_covariates <- c("age", "sex", "charlson", "vaso_day1", "sofa", "wbc", "lact", "ph", "bili", "crea", "pf", "hr_max", "temp_min")
  reduced_covariates <- c("age", "sex", "sofa", "lact")
  
  covariates <- if (study_name == "CMAISE") reduced_covariates else full_covariates
  
  df_complete <- df %>% select(all_of(c("time_to_event", "event_status_cr", "death", "ARC", covariates))) %>% na.omit()
  
  if(nrow(df_complete) < 100) return(NULL)
  
  # --- Cox Model ---
  cox_model <- coxph(as.formula(paste("Surv(time_to_event, death) ~ ARC +", paste(covariates, collapse = " + "))), data = df_complete)
  hr <- exp(coef(cox_model)["ARC"])
  hr_ci <- exp(confint.default(cox_model)["ARC", ])
  
  cox_result <- data.frame(
    Subgroup = "  Standard Cox Model",
    Model = "Cox",
    OR = hr, Lower = hr_ci[1], Upper = hr_ci[2]
  )
  
  # --- Competing Risk Model ---
  formula_str <- paste("~", paste(c("ARC", covariates), collapse = " + "))
  cov_matrix <- model.matrix(as.formula(formula_str), data = df_complete)[, -1]
  
  cr_model <- crr(ftime = df_complete$time_to_event, fstatus = df_complete$event_status_cr, 
                  cov1 = cov_matrix, failcode = 1, cencode = 2)
  
  summary_cr <- summary(cr_model)
  shr <- summary_cr$coef["ARC", "exp(coef)"]
  ci_low <- summary_cr$conf.int["ARC", "2.5%"]
  ci_high <- summary_cr$conf.int["ARC", "97.5%"]
  
  cr_result <- data.frame(
    Subgroup = "  Competing Risk Model",
    Model = "CRR",
    OR = shr, Lower = ci_low, Upper = ci_high
  )
  
  # 添加数据库标题行
  header_row <- data.frame(Subgroup = study_name, Model = "Header", OR = NA, Lower = NA, Upper = NA)
  
  return(rbind(header_row, cox_result, cr_result))
}

# --- Step 4: 运行分析并整合所有结果 ---
all_results <- lapply(names(dat_list), function(name) {
  run_both_models(dat_list[[name]], name)
})

plot_data <- do.call(rbind, all_results)

# --- Step 5: 准备用于forestploter的最终数据框 ---
plot_data <- plot_data %>%
  mutate(
    `Effect Size (95% CI)` = ifelse( # More generic header
      !is.na(OR),
      sprintf("%.2f (%.2f-%.2f)", OR, Lower, Upper),
      ""
    ),
    ` ` = paste(rep(" ", 20), collapse = ""), # 空白列用于绘图
    is_header = Model == "Header"
  )

# --- Step 6: 【最终优化】生成出版级质量的森林图 ---
tm <- forest_theme(
  base_size = 10,
  row_height = unit(0.8, "cm"),
  # **【错误修复】**: 定义图例和对应的美学
  legend_name = "Model Type",
  legend_value = c("Standard Cox (HR)", "Competing Risk (SHR)"),
  ci_pch = c(16, 15), # Circle for Cox, Square for CRR
  ci_col = c("#2166AC", "#B2182B"), # Blue for Cox, Red for CRR
  refline_lwd = 1, refline_lty = "dashed", refline_col = "grey20",
  colhead_just = c("left", "center", "center")
)

p <- forest(
  data = plot_data[, c("Subgroup", "Effect Size (95% CI)", " ")],
  est = plot_data$OR,
  lower = plot_data$Lower,
  upper = plot_data$Upper,
  # **【错误修复】**: 添加 legend 参数来激活分组
  legend = plot_data$Model,
  ci_column = 3,
  ref_line = 1,
  arrow_lab = c("ARC Protective", "ARC Harmful"),
  xlim = c(0.4, 3.0), # Adjusted for better visualization
  # 设置列宽
  widths = unit(c(5, 4, 5), "cm"),
  title = "Sensitivity Analysis: Standard Cox vs. Competing Risk Models",
  theme = tm
)

# 为标题行加粗
p <- edit_plot(p, row = which(plot_data$is_header), gp = gpar(fontface = "bold"))
# 移除标题行的绘图元素
p <- edit_plot(p, row = which(plot_data$is_header), col = 3, which = "ci", gp = gpar(col = "white"))

# --- Step 7: 保存图表 ---
ggsave(
  "../figure/FigureS_Sensitivity_CompetingRisks_Forest.pdf",
  plot = p, device = "pdf", width = 9, height = 11
)
ggsave(
  "../figure/FigureS_Sensitivity_CompetingRisks_Forest.tiff",
  plot = p, device = "tiff", width = 9, height = 11, units = "in", dpi = 300, compression = "lzw"
)

print("Competing risks sensitivity analysis forest plot has been saved to the '../figure/' directory.")

