# --- Step 1: Load Libraries ---
# 确保所有必要的包都已安装和加载
# install.packages(c("data.table", "dplyr", "MatchIt", "survival", "survminer", "patchwork", "Cairo"))
library(data.table)
library(dplyr)
library(survival)
library(survminer)
library(patchwork) # 用于组合图表
library(Cairo)

# --- Step 2: Load Data and Assume PSM Dataframes Exist ---
# 这部分代码假设您已经运行了PSM脚本，并且环境中存在以下数据框：
# PSM_miiv, PSM_mimic, PSM_eicu, PSM_aumc, PSM_CMAISE
#
# 如果没有，请先运行您之前的PSM代码来生成它们。
# 示例：
# dat <- fread("d:/OneDrive/R/2024Sepsis_ARC/DataAnalysis/RICU_sepsis_com3.csv")
# ... (此处省略您所有的matchit代码) ...
# PSM_miiv <- match.data(m.out_miiv)
# ... etc ...

# --- Step 3: Create a Professional, Reusable KM Plot Function ---
# 这个函数是核心，它统一了所有图的风格，并创建了正确的生存对象

create_km_plot <- function(df, title, show_risk_table = TRUE, time_limit_days = 28) {
  
  # **关键修正: 创建正确的28天生存对象**
  df <- df %>%
    mutate(
      # 生存时间：取住院天数和28天中的较小值
      time_to_event = pmin(los_hosp, time_limit_days, na.rm = TRUE),
      # 结局事件：仅当在28天内死亡时，事件为1
      event_status = ifelse(death == 1 & los_hosp <= time_limit_days, 1, 0)
    )
  
  # 创建生存拟合对象
  fit <- survfit(Surv(time_to_event, event_status) ~ ARC, data = df)
  
  # 绘制KM图
  p <- ggsurvplot(
    fit,
    data = df,
    # --- 美学风格 ---
    palette = c("#B2182B", "#2166AC"), # 专业、色觉友好的颜色 (暗红, 深蓝)
    ggtheme = theme_classic(base_size = 12),
    size = 1.1, # **美化**: 略微加粗线条
    
    # **【新增加】**: 限制纵轴范围为0.5到1
    ylim = c(0.5, 1),
    
    # --- 标签和图例 ---
    title = title,
    xlab = paste0("Time in Days (Censored at ", time_limit_days, " Days)"),
    ylab = "Overall Survival Probability",
    legend.title = "Group",
    legend.labs = c("No ARC", "ARC"),
    
    # --- P值 ---
    pval = TRUE,
    pval.method = TRUE, # **美化**: 显示检验方法和更规范的P值
    pval.coord = c(1, 0.55), # **【新调整】**: 将P值移到更醒目的位置
    pval.size = 4,
    
    # --- 置信区间 ---
    conf.int = TRUE,
    conf.int.style = "step",
    
    # --- 风险表 (根据参数决定是否显示) ---
    risk.table = show_risk_table,
    risk.table.y.text = FALSE,  
    risk.table.title = "Number at Risk",
    risk.table.fontsize = 3.5,
    tables.theme = theme_survminer(font.main = 10),
    
    # --- 其他细节 ---
    censor.size = 3,
    censor.shape = "|" # **美化**: 使用更简洁的竖线作为删失符号
  )
  
  # 如果显示风险表，调整图和表的相对高度
  if (show_risk_table) {
    p$plot <- p$plot + theme(plot.title = element_text(hjust = 0.5, face = "bold"),
                             legend.position = "top") + 
      # **【新增加】**: 在图的主体部分也应用ylim
      coord_cartesian(ylim = c(0.5, 1))
    p <- ggarrange(p$plot, p$table, heights = c(0.7, 0.3), ncol = 1, nrow = 2)
  } else {
    # 如果不显示风险表，p本身就是ggplot对象
    p <- p$plot + theme(plot.title = element_text(hjust = 0.5, face = "bold"),
                        legend.position = "top") +
      # **【新增加】**: 在图的主体部分也应用ylim
      coord_cartesian(ylim = c(0.5, 1))
  }
  
  return(p)
}

# --- Step 4: Generate the Main Figure for the Manuscript (Figure 3) ---

# 主图: MIMIC-IV, 带有风险表
p_main <- create_km_plot(PSM_miiv, title = "Survival in the Primary Cohort (MIMIC-IV)", show_risk_table = TRUE)

# --- Step 5: Save the Main Figure ---
# PDF (矢量图，用于排版)
ggsave(
  "../figure/Figure3_KM_Primary_Cohort.pdf",
  plot = p_main,
  device = "pdf",
  width = 8,
  height = 7,
  units = "in",
  useDingbats = FALSE
)

# TIFF (高分辨率位图，用于提交)
ggsave(
  "../figure/Figure3_KM_Primary_Cohort.tiff",
  plot = p_main,
  device = "tiff",
  width = 8,
  height = 7,
  units = "in",
  dpi = 300,
  compression = "lzw"
)

print("Figure 3 (Primary Cohort KM Plot) has been saved as PDF and TIFF.")


# --- Step 6: Generate and Optimize the Supplementary Figure for Validation Cohorts ---

# 为4个验证队列生成单独的图，不带风险表
p_mimic3 <- create_km_plot(PSM_mimic, title = "MIMIC-III", show_risk_table = FALSE)
p_eicu   <- create_km_plot(PSM_eicu, title = "eICU", show_risk_table = FALSE)
p_aumc   <- create_km_plot(PSM_aumc, title = "AUMC", show_risk_table = FALSE)
p_cmaise <- create_km_plot(PSM_CMAISE, title = "CMAISE", show_risk_table = FALSE)

# **【美化核心】: 优化组合图，共享坐标轴，采用2x2布局**

# 移除顶部图的X轴
p_mimic3 <- p_mimic3 + theme(axis.title.x = element_blank(), axis.text.x = element_blank(), legend.position = "none")
p_eicu   <- p_eicu   + theme(axis.title.x = element_blank(), axis.text.x = element_blank(), legend.position = "none")

# 移除右侧图的Y轴
p_eicu   <- p_eicu   + theme(axis.title.y = element_blank(), axis.text.y = element_blank())
p_cmaise <- p_cmaise + theme(axis.title.y = element_blank(), axis.text.y = element_blank(), legend.position = "none")

# 移除左下角图的图例
p_aumc <- p_aumc + theme(legend.position = "none")


# 使用 patchwork 组合成2x2网格
validation_plot <- (p_mimic3 | p_eicu) / (p_aumc | p_cmaise)

# 添加一个总标题和共享的图例
validation_plot_final <- validation_plot + 
  plot_layout(guides = "collect") + # 收集所有图例并只显示一个
  plot_annotation(
    title = "Survival in the Validation Cohorts",
    caption = "Kaplan-Meier survival curves censored at 28 days for each validation cohort after propensity score matching.",
    theme = theme(plot.title = element_text(hjust = 0.5, face = "bold", size=16))
  ) & theme(legend.position = 'bottom') # 将统一的图例放在底部

# --- Step 7: Save the Supplementary Figure ---
ggsave(
  "../figure/FigureS_KM_Validation_Cohorts.pdf",
  plot = validation_plot_final,
  device = "pdf",
  width = 10, # 调整为更适合2x2布局的尺寸
  height = 8,
  units = "in",
  useDingbats = FALSE
)

ggsave(
  "../figure/FigureS_KM_Validation_Cohorts.tiff",
  plot = validation_plot_final,
  device = "tiff",
  width = 10,
  height = 8,
  units = "in",
  dpi = 300,
  compression = "lzw"
)

print("Supplementary Figure (Validation Cohorts KM Plot) has been saved as PDF and TIFF.")