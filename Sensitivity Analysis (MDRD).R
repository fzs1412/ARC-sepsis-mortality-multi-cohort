# --- Step 1: Load Libraries and Prepare Data ---
# 确保所有必要的包都已安装和加载
#  (Random Effects Model)
library(data.table)
library(dplyr)
library(meta)
library(Cairo)
library(forestploter)
library(grid)

# 加载您的主数据集
# Main data loading section remains the same as your provided code
dat <- fread("d:/OneDrive/R/2024Sepsis_ARC/DataAnalysis/RICU_sepsis_com3.csv") %>%
  select(source, age, sex, charlson,
         Hypertension, Diabete, CHF, MI, COPD, sofa, vaso_day1, 
         wbc, plt, ph, pf, lact, bili, crea, CrCl, bun,
         hr_max, hr_min, map_max, map_min,
         resp_max, resp_min, temp_max, temp_min,
         death, ARC_MDRD, MV_time, los_icu, los_hosp)


# --- Step 2: Adaptive Regression Function ---
# This section (run_regression_analysis function) remains the same as your code

run_regression_analysis <- function(df, study_name) {
  
  full_formula <- death ~ ARC_MDRD + age + sex + charlson + vaso_day1 + sofa + wbc + lact + ph + bili + crea + pf + hr_max + temp_min
  reduced_formula <- death ~ ARC_MDRD + age + sex + sofa + lact
  
  if (study_name == "CMAISE") {
    formula <- reduced_formula
    message(paste("Using REDUCED model for", study_name))
  } else {
    formula <- full_formula
    message(paste("Using FULL model for", study_name))
  }
  
  # 在运行模型前处理缺失值
  model_vars <- all.vars(formula)
  df_complete <- df %>% select(all_of(model_vars)) %>% na.omit()
  
  model <- tryCatch({
    glm(formula, data = df_complete, family = "binomial")
  }, error = function(e) {
    message(paste("Error in", study_name, ":", e$message)); return(NULL)
  })
  
  if (is.null(model)) return(NULL)
  
  if ("ARC_MDRD" %in% names(coef(model))) {
    summary_model <- summary(model)
    conf_intervals <- tryCatch(confint(model), error = function(e) confint.default(model))
    
    or <- exp(coef(model)["ARC_MDRD"])
    ci_low <- exp(conf_intervals["ARC_MDRD", 1])
    ci_high <- exp(conf_intervals["ARC_MDRD", 2])
    
    # 从完整数据中获取事件计数以进行报告
    events_arc <- sum(df$death[df$ARC_MDRD == 1], na.rm = TRUE)
    total_arc <- sum(df$ARC_MDRD == 1, na.rm = TRUE)
    events_no_arc <- sum(df$death[df$ARC_MDRD == 0], na.rm = TRUE)
    total_no_arc <- sum(df$ARC_MDRD == 0, na.rm = TRUE)
    
    return(data.frame(
      study = study_name,
      N = nrow(df), # 添加总样本量
      Events = events_arc + events_no_arc, # 添加总事件数
      or = or,
      ci.lo = ci_low,
      ci.hi = ci_high
    ))
  } else {
    return(NULL)
  }
}


# --- Step 3: Run Analysis and Prepare Plotting Data ---
# This section remains largely the same, but we will refine the final plot_data

# Splitting dataframes
dat_aumc <- dat[dat$source == "aumc", ]
dat_eicu <- dat[dat$source == "eicu", ]
dat_miiv <- dat[dat$source == "miiv", ]
dat_mimic <- dat[dat$source == "mimic", ]

# CMAISE data prep
dat_CMAISE <- dat[dat$source == "CMAISE1.5v", ] %>%
  mutate(
    charlson_proxy = Hypertension + Diabete + CHF + MI + COPD,
    charlson = charlson_proxy
  ) %>%
  select(
    source, age, sex, charlson, vaso_day1, sofa, wbc, lact, ph, bili, crea, pf,
    hr_max, temp_min, death, ARC_MDRD
  )

results_list_glm <- list(
  run_regression_analysis(dat_miiv, "MIMIC-IV"),
  run_regression_analysis(dat_mimic, "MIMIC-III"),
  run_regression_analysis(dat_eicu, "eICU"),
  run_regression_analysis(dat_aumc, "AUMC"),
  run_regression_analysis(dat_CMAISE, "CMAISE")
)

final_results_glm <- do.call(rbind, results_list_glm)

meta_analysis_glm <- metagen(
  TE = log(or),
  lower = log(ci.lo),
  upper = log(ci.hi),
  studlab = study,
  data = final_results_glm,
  sm = "OR",
  comb.random = TRUE
)

# --- Step 4: 【最终优化】 Prepare Dataframe for `forestploter` ---

# 添加总结行的数据
summary_data <- data.frame(
  study = "Overall ", # **【恢复模型标签】**
  N = sum(final_results_glm$N),
  Events = sum(final_results_glm$Events),
  or = exp(meta_analysis_glm$TE.random),
  ci.lo = exp(meta_analysis_glm$lower.random),
  ci.hi = exp(meta_analysis_glm$upper.random)
)

# 合并研究数据和总结数据
plot_data <- rbind(final_results_glm, summary_data)

# 格式化用于显示的文本列
plot_data <- plot_data %>%
  mutate(
    N = format(N, big.mark = ","),
    Events = format(Events, big.mark = ","),
    # **【解决重叠】**: 将文本和图形分开
    `OR (95% CI)` = sprintf("%.2f (%.2f-%.2f)", or, ci.lo, ci.hi),
    # 空白列用于容纳森林图的图形部分
    ` ` = paste(rep(" ", 20), collapse = "")
  )

# --- Step 5: 【最终优化】 Generate Publication-Quality Forest Plot ---
# 定义主题
tm <- forest_theme(
  base_size = 10,
  row_height = unit(0.9, "cm"),
  ci_pch = 16, ci_col = "#377eb8", ci_lty = 1, ci_lwd = 1.5,
  refline_lwd = 1, refline_lty = "dashed", refline_col = "grey20",
  summary_pch = 18, summary_col = "#e41a1c",
  # 对齐列标题
  colhead_just = c("left", "right", "right", "right", "center")
)

# 绘制森林图
p <- forest(
  data = plot_data[, c("study", "N", "Events", "OR (95% CI)", " ")],
  est = plot_data$or,
  lower = plot_data$ci.lo,
  upper = plot_data$ci.hi,
  sizes = 0.8,
  is_summary = (1:nrow(plot_data)) == nrow(plot_data),
  # **【解决重叠】**: ci_column现在指向空白列
  ci_column = 5,
  ref_line = 1,
  arrow_lab = c("ARC_MDRD Protective", "ARC_MDRD Harmful"),
  xlim = c(0.2, 2.5), # 调整X轴范围以适应CMAISE的宽CI
  # **【解决重叠】**: 精确设置列宽
  widths = unit(c(4.5, 1.5, 1.5, 3, 3.5), "cm"), # Study, N, Events, OR(text), Plot
  #title = "Association of ARC_MDRD with Hospital Mortality (Multivariable Model)",
  theme = tm
)

# 为总结行加粗
p <- edit_plot(p, row = nrow(plot_data), gp = gpar(fontface = "bold"))

# --- Step 6: 保存图表 ---
ggsave(
  "../figure/FigureS_ForestPlot_MDRD.pdf",
  plot = p, device = "pdf", width = 12, height = 7
)

ggsave(
  "../figure/FigureS_ForestPlot_MDRD.tiff",
  plot = p, device = "tiff", width = 12, height = 7, units = "in", dpi = 300, compression = "lzw"
)

print("Figure 2 (using adaptive models) has been saved to the '../figure/' directory.")
