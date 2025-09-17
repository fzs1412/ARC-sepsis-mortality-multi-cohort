# --- Step 1: Load Libraries and Prepare Data ---
# 确保所有必要的包都已安装和加载
# install.packages(c("data.table", "dplyr", "MatchIt", "meta", "Cairo"))
library(data.table)
library(dplyr)
library(MatchIt)
library(meta)
library(Cairo)

# 加载您的主数据集
dat <- fread("d:/OneDrive/R/2024Sepsis_ARC/DataAnalysis/RICU_sepsis_com3.csv")

# 将数据拆分为单独的数据框
dat_aumc <- dat[dat$source == "aumc", ]
dat_eicu <- dat[dat$source == "eicu", ]
dat_miiv <- dat[dat$source == "miiv", ]
dat_mimic <- dat[dat$source == "mimic", ]
dat_CMAISE <- dat[dat$source == "CMAISE1.5v", ]

# --- Step 2: 对每个数据库执行倾向性评分匹配(PSM) ---

# --- CMAISE Cohort ---
# 由于CMAISE缺少charlson，我们使用单个合并症进行匹配
message("--- Starting PSM for CMAISE ---")
m.out_cmaise <- matchit(ARC ~ age + sex + vaso_day1 + sofa + wbc + lact + ph + bili + crea + pf +
                          hr_max + temp_min + Hypertension + Diabete + CHF + MI + COPD,
                        data = dat_CMAISE, method = "nearest", ratio = 1)
summary(m.out_cmaise)
# plot(m.out_cmaise, type = "jitter") # 可选：检查匹配质量
PSM_CMAISE <- match.data(m.out_cmaise)
message(paste("CMAISE matched pairs:", nrow(PSM_CMAISE) / 2))
# fwrite(PSM_CMAISE, "PSM_CMAISE.csv") # 可选：保存匹配后的数据

# --- MIMIC-III Cohort ---
message("--- Starting PSM for MIMIC-III ---")
m.out_mimic <- matchit(ARC ~ age + sex + vaso_day1 + sofa + wbc + lact + ph + bili + crea + pf +
                         hr_max + temp_min + charlson,
                       data = dat_mimic, method = "nearest", ratio = 1)
summary(m.out_mimic)
PSM_mimic <- match.data(m.out_mimic)
message(paste("MIMIC-III matched pairs:", nrow(PSM_mimic) / 2))

# --- MIMIC-IV Cohort ---
message("--- Starting PSM for MIMIC-IV ---")
m.out_miiv <- matchit(ARC ~ age + sex + vaso_day1 + sofa + wbc + lact + ph + bili + crea + pf +
                        hr_max + temp_min + charlson,
                      data = dat_miiv, method = "nearest", ratio = 1)
summary(m.out_miiv)
PSM_miiv <- match.data(m.out_miiv)
message(paste("MIMIC-IV matched pairs:", nrow(PSM_miiv) / 2))

# --- eICU Cohort (续写部分) ---
message("--- Starting PSM for eICU ---")
m.out_eicu <- matchit(ARC ~ age + sex + vaso_day1 + sofa + wbc + lact + ph + bili + crea + pf +
                        hr_max + temp_min + charlson,
                      data = dat_eicu, method = "nearest", ratio = 1)
summary(m.out_eicu)
PSM_eicu <- match.data(m.out_eicu)
message(paste("eICU matched pairs:", nrow(PSM_eicu) / 2))

# --- AUMC Cohort (续写部分) ---
message("--- Starting PSM for AUMC ---")
m.out_aumc <- matchit(ARC ~ age + sex + vaso_day1 + sofa + wbc + lact + ph + bili + crea + pf +
                        hr_max + temp_min + charlson,
                      data = dat_aumc, method = "nearest", ratio = 1)
summary(m.out_aumc)
PSM_aumc <- match.data(m.out_aumc)
message(paste("AUMC matched pairs:", nrow(PSM_aumc) / 2))


# --- Step 3: 在每个匹配后的队列中分析结局 ---

# 创建一个函数来简化结局分析
analyze_matched_outcome <- function(matched_df, study_name) {
  # 检查数据是否有效
  if (is.null(matched_df) || nrow(matched_df) == 0) {
    message(paste("No matched data for", study_name))
    return(NULL)
  }
  
  # 在匹配后的数据上运行逻辑回归
  # 'death' 是您的住院死亡率变量
  model <- glm(death ~ ARC, data = matched_df, family = "binomial")
  
  # 提取结果
  summary_model <- summary(model)
  conf_intervals <- confint.default(model) # 使用默认CI以避免潜在错误
  
  or <- exp(coef(model)["ARC"])
  ci_low <- exp(conf_intervals["ARC", 1])
  ci_high <- exp(conf_intervals["ARC", 2])
  
  # 获取事件计数
  events_arc <- sum(matched_df$death[matched_df$ARC == 1], na.rm = TRUE)
  total_arc <- sum(matched_df$ARC == 1, na.rm = TRUE)
  events_no_arc <- sum(matched_df$death[matched_df$ARC == 0], na.rm = TRUE)
  total_no_arc <- sum(matched_df$ARC == 0, na.rm = TRUE)
  
  return(data.frame(
    study = study_name,
    or = or,
    ci.lo = ci_low,
    ci.hi = ci_high,
    events.arc = events_arc,
    total.arc = total_arc,
    events.noarc = events_no_arc,
    total.noarc = total_no_arc
  ))
}

# 为每个匹配后的数据集运行结局分析
results_list <- list(
  analyze_matched_outcome(PSM_miiv, "MIMIC-IV"),
  analyze_matched_outcome(PSM_mimic, "MIMIC-III"),
  analyze_matched_outcome(PSM_eicu, "eICU"),
  analyze_matched_outcome(PSM_aumc, "AUMC"),
  analyze_matched_outcome(PSM_CMAISE, "CMAISE")
)

# 合并所有结果
final_results <- do.call(rbind, results_list)

# --- Step 4: 执行Meta分析并生成出版级森林图 ---

if (!is.null(final_results) && nrow(final_results) > 0) {
  
  # 添加用于在图表中显示的列
  final_results <- final_results %>%
    mutate(
      n_arc_display = paste0(events.arc, "/", total.arc),
      n_noarc_display = paste0(events.noarc, "/", total.noarc)
    )
  
  print("--- 来自PSM分析的最终结果 ---")
  print(final_results)
  
  # 执行Meta分析
  meta_analysis <- metagen(
    TE = log(or),
    lower = log(ci.lo),
    upper = log(ci.hi),
    studlab = study,
    data = final_results,
    sm = "OR",
    comb.random = TRUE, # 推荐使用随机效应模型
    title = "Association of ARC with Hospital Mortality after PSM"
  )
  
  # 创建绘图函数
  create_publication_forest_plot <- function(meta_obj, data_for_plot) {
    forest(
      meta_obj,
      leftcols = c("studlab", "n_arc_display", "n_noarc_display", "effect.ci", "w.random"),
      leftlabs = c("Study (PSM Cohort)", "Events/N (ARC)", "Events/N (No ARC)", "Odds Ratio (95% CI)", "Weight (%)"),
      rightcols = FALSE,
      col.square = "black",
      col.diamond = "black",
      prediction = FALSE,
      xlab = "OR for Hospital Mortality (ARC vs. No ARC)",
      smlab = "Overall Effect",
      spacing = 1.2,
      print.I2 = TRUE,
      print.tau2 = TRUE
    )
  }
  
  # --- 生成PDF和TIFF格式的图片 ---
  pdf("../figure/FigureS_ForestPlot_PSM_HospitalDeath_Final.pdf", width = 11, height = 7)
  create_publication_forest_plot(meta_analysis, final_results)
  dev.off()
  
  tiff("../figure/FigureS_ForestPlot_PSM_HospitalDeath_Final.tiff", width = 11, height = 7, units = "in", res = 300, compression = "lzw")
  create_publication_forest_plot(meta_analysis, final_results)
  dev.off()
  
  print("来自PSM的补充森林图(PDF和TIFF格式)已保存到 '../figure/' 目录。")
} else {
  print("没有有效的匹配结果来生成Meta分析图表。")
}

