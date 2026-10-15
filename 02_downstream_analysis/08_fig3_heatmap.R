# ==============================================================================
# 🏆 SCI Figure 3: 核心植物科跨方法检出频率全景热图 (Nature暖红 + 零值隐藏)
# ==============================================================================
suppressMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(patchwork)
})

DATA_DIR <- Sys.getenv("MUNTJAC_DATA", unset=".")
setwd(file.path(DATA_DIR, "02_abundance_matrices"))

# 固定的 7 维度分析策略
methods <- c("Traditional", "DIR_Loose", "DIR_Medium", "DIR_Strict", "ASM_90", "ASM_95", "ASM_99")
methods_label <- c("Traditional"="Literature-based",
                   "DIR_Loose"="DIR_Loose","DIR_Medium"="DIR_Medium","DIR_Strict"="DIR_Strict",
                   "ASM_90"="ASM_90","ASM_95"="ASM_95","ASM_99"="ASM_99")

cat("=================================================================\n")
cat("🚀 启动 Figure 3 核心植物全景热图绘制引擎 (暖红定制版)...\n")
cat("=================================================================\n")

# 提取数据的函数
get_occurrence_data <- function(host_cn, top_n) {
  csv_file <- ifelse(host_cn=="小麂", "Matrix_detections_Reeves.csv", "Matrix_detections_crinifrons.csv")
  if (!file.exists(csv_file)) stop(sprintf("❌ 找不到文件: %s", csv_file))
  
  df <- read.csv(csv_file, check.names = FALSE)
  df$Total_Abundance <- rowSums(df[, methods], na.rm = TRUE)
  
  df_top <- df %>% slice_max(order_by = Total_Abundance, n = top_n, with_ties = FALSE)
  
  df_long <- df_top %>% 
    select(Family_Latin, all_of(methods)) %>%
    pivot_longer(cols = all_of(methods), names_to = "Method", values_to = "Frequency")
  
  df_long$Method <- factor(methods_label[as.character(df_long$Method)], levels = methods_label[methods])
  df_long$Family_Latin <- factor(df_long$Family_Latin, levels = rev(df_top$Family_Latin))
  
  return(df_long)
}

# 核心绘图函数
plot_occurrence_heatmap <- function(df_long, global_max) {
  
  # 计算颜色反转阈值 (高热度区域使用白字，低热度区域使用黑字/深灰)
  text_color_threshold <- global_max * 0.5
  
  p <- ggplot(df_long, aes(x = Method, y = Family_Latin, fill = Frequency)) +
    geom_tile(color = "white", linewidth = 0.6) +
    # 🎨 【核心修改1】零值留白 (Frequency > 0 才显示)，并根据底色深浅自动反转文字颜色
    geom_text(aes(label = ifelse(Frequency > 0, round(Frequency, 0), ""), 
                  color = Frequency > text_color_threshold), 
              size = 2.5, fontface = "bold") +
    scale_color_manual(values = c("TRUE" = "white", "FALSE" = "grey20"), guide = "none") +
    # 🎨 【核心修改2】采用你指定的同款 Nature 暖红渐变配色
    scale_fill_gradient(low = "#F5F5F5", high = "#E64B35", 
                        name = "Occurrence\nFrequency", 
                        limits = c(0, global_max), 
                        na.value = "white") +
    coord_fixed(ratio = 0.53) +
    labs(x = "", y = "", title = NULL) + 
    theme_minimal(base_family = "sans") +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, color = "black", face = "bold", size = 8),
      axis.text.y = element_text(color = "black", face = "italic", size = 7),
      panel.grid = element_blank(),
      plot.margin = margin(3, 6, 3, 6)
    )
  return(p)
}

# ==============================================================================
# 🧩 批量生成并组装 (Top 25 和 Top 20 独立出图)
# ==============================================================================
top_levels <- c(25, 20)

for (n in top_levels) {
  cat(sprintf("⏳ 正在生成 Top %d 的小麂与黑麂对比组图...\n", n))
  
  # 1. 提取小麂和黑麂的数据
  df_reevesi    <- get_occurrence_data("小麂", n)
  df_crinifrons <- get_occurrence_data("黑麂", n)
  
  # 2. 计算全局最大值，确保两张图使用绝对一致的颜色标尺
  global_max <- max(c(df_reevesi$Frequency, df_crinifrons$Frequency), na.rm = TRUE)
  
  # 3. 绘制单图
  p_reevesi    <- plot_occurrence_heatmap(df_reevesi, global_max)
  p_crinifrons <- plot_occurrence_heatmap(df_crinifrons, global_max)
  
  # 4. Patchwork 左右拼图 (左小麂，右黑麂)
  final_plot <- p_reevesi | p_crinifrons
  
  # 5. 剥离独立图例，统一放置，并加上 A B 标签
  final_plot <- final_plot + 
    plot_annotation(tag_levels = 'A') + 
    plot_layout(guides = "collect") & 
    theme(
      plot.tag = element_text(size = 18, face = "bold"),
      legend.position = "right",
      legend.title = element_text(face = "bold", size = 10),
      legend.key.height = unit(1.5, "cm")
    )
  
  # 6. 导出 PDF
  out_pdf <- file.path("..", "figures", sprintf("Fig3_Heatmap_Top%d.pdf", n)); dir.create(dirname(out_pdf), showWarnings=FALSE, recursive=TRUE)
  plot_height <- ifelse(n == 25, 9, 8) 
  ggsave(out_pdf, final_plot, width = 210, height = 146, units = "mm", dpi = 600)
  svglite::svglite(file.path("..", "figures", sprintf("Fig3_Heatmap_Top%d.svg", n)), width = 210/25.4, height = 146/25.4); print(final_plot); dev.off()
  
  cat(sprintf("  ├─ ✅ 生成成功: %s (全局最大标尺: %d)\n", out_pdf, round(global_max)))
}

cat("\n🎉 [大功告成] Nature 暖红配色、零值清爽留白的 Top 25/20 全景热图已生成！\n")