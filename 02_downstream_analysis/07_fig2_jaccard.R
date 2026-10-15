# ==============================================================================
# 🏆 SCI 巨幅排版：多层级食谱一致性矩阵图版 (Patchwork + 极致细节优化版)
# ==============================================================================
suppressMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(patchwork) 
})

DATA_DIR <- Sys.getenv("MUNTJAC_DATA", unset=".")
setwd(file.path(DATA_DIR, "02_abundance_matrices"))

# 定义分析参数
methods <- c("Traditional", "DIR_Loose", "DIR_Medium", "DIR_Strict", "ASM_90", "ASM_95", "ASM_99")
# 显示标签：Traditional 显示为 literature-based（读数据仍用原列名）
methods_label <- c("Traditional"="Literature-based",
                   "DIR_Loose"="DIR_Loose","DIR_Medium"="DIR_Medium","DIR_Strict"="DIR_Strict",
                   "ASM_90"="ASM_90","ASM_95"="ASM_95","ASM_99"="ASM_99")
hosts_cn <- c("小麂", "黑麂")
hosts_latin <- c("M. reevesi", "M. crinifrons")

cat("=================================================================\n")
cat("🚀 正在构建 A-F 面板的 Jaccard 一致性 SCI 全景图版...\n")
cat("=================================================================\n")

# 核心绘图函数
generate_jaccard_panel <- function(host_cn, host_latin, top_n, panel_label) {
  csv_file <- ifelse(host_cn=="小麂", "Matrix_detections_Reeves.csv", "Matrix_detections_crinifrons.csv")
  if (!file.exists(csv_file)) stop(sprintf("找不到文件: %s", csv_file))
  
  df <- read.csv(csv_file, check.names = FALSE)
  df$Total_Abundance <- rowSums(df[, methods], na.rm = TRUE)
  
  if (top_n == "All") {
    df_top <- df
    title_suffix <- "All detected families"
  } else {
    df_top <- df %>% slice_max(order_by = Total_Abundance, n = as.numeric(top_n), with_ties = FALSE)
    title_suffix <- paste0("Top ", top_n, " families")
  }
  
  df_bin <- df_top %>% mutate(across(all_of(methods), ~ifelse(. > 0, 1, 0)))
  jaccard_matrix <- matrix(nrow = length(methods), ncol = length(methods), dimnames = list(methods, methods))
  
  for (i in 1:length(methods)) {
    for (j in 1:length(methods)) {
      set1 <- df_bin$Family_Latin[df_bin[[methods[i]]] == 1]
      set2 <- df_bin$Family_Latin[df_bin[[methods[j]]] == 1]
      shared <- length(intersect(set1, set2))
      union_set <- length(union(set1, set2))
      jaccard_matrix[i, j] <- ifelse(union_set == 0, 0.00, round((shared / union_set) * 100, 2))
    }
  }
  
  df_jaccard <- as.data.frame(as.table(jaccard_matrix))
  colnames(df_jaccard) <- c("Method1", "Method2", "Jaccard")
  
  # 提取下三角（直角在左下，对角线从左上到右下）
  df_jaccard$idx1 <- match(df_jaccard$Method1, methods)
  df_jaccard$idx2 <- match(df_jaccard$Method2, methods)
  df_jaccard <- df_jaccard %>% filter(idx1 <= idx2) 
  
  df_jaccard$Method1 <- factor(methods_label[as.character(df_jaccard$Method1)], levels = methods_label[methods])
  df_jaccard$Method2 <- factor(methods_label[as.character(df_jaccard$Method2)], levels = rev(methods_label[methods]))
  
  # 🎨 SCI 细节调优
  p <- ggplot(df_jaccard, aes(x = Method1, y = Method2, fill = Jaccard)) +
    geom_tile(color = "white", linewidth = 0.6) +
    # 【优化 1】字体大小设为 3，阈值设为 55%，浅色底配深灰(grey20)，深色底配纯白
    geom_text(aes(label = sprintf("%.2f%%", Jaccard), color = Jaccard > 55), 
              size = 1.8, fontface = "bold") +
    scale_color_manual(values = c("TRUE" = "white", "FALSE" = "grey20"), guide = "none") +
    scale_fill_gradient(low = "#F5F5F5", high = "#1565C0", name = "Jaccard\nSimilarity (%)", limits = c(0, 100), na.value = "white") +
    coord_fixed() + 
    labs(
      title = substitute(italic(host) ~ " (" ~ suffix ~ ")", list(host = host_latin, suffix = title_suffix)),
      x = "", y = ""
    ) +
    theme_minimal(base_family = "sans") +
    theme(
      # 【优化 2】X轴 45度斜排，同时 vjust=1 确保文字右上角完美贴合刻度线
      axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, color = "black", face = "bold", size = 7),
      # 【优化 3】Y轴 绝对水平，vjust=0.5 确保垂直居中对齐方块
      axis.text.y = element_text(angle = 0, hjust = 1, vjust = 0.5, color = "black", face = "bold", size = 7),
      plot.title = element_text(hjust = 0.5, size = 9, face = "bold"),
      panel.grid = element_blank(),
      plot.margin = margin(2, 2, 2, 2)
    )
  
  return(p)
}

# 独立生成 6 个面板 (A - F)
pA <- generate_jaccard_panel("小麂", "M. reevesi", "All", "A")
pB <- generate_jaccard_panel("黑麂", "M. crinifrons", "All", "B")
pC <- generate_jaccard_panel("小麂", "M. reevesi", "20", "C")
pD <- generate_jaccard_panel("黑麂", "M. crinifrons", "20", "D")
pE <- generate_jaccard_panel("小麂", "M. reevesi", "10", "E")
pF <- generate_jaccard_panel("黑麂", "M. crinifrons", "10", "F")

# 组装并渲染最终图版
final_plot <- (pA | pB) / (pC | pD) / (pE | pF) +
  plot_annotation(
    title = "Consistency of dietary reconstruction across thresholds and analytical strategies",
    tag_levels = 'A', 
    theme = theme(
      plot.title = element_text(size = 12, face = "bold", hjust = 0.5, margin = margin(b = 8)),
      plot.tag = element_text(size = 12, face = "bold") 
    )
  ) +
  plot_layout(guides = "collect") & 
  theme(legend.position = "right", legend.key.height = unit(0.8, "cm"))

out_pdf <- file.path("..", "figures", "Fig2_Jaccard.pdf"); dir.create(dirname(out_pdf), showWarnings=FALSE, recursive=TRUE)
ggsave(out_pdf, final_plot, width = 210, height = 297, units = "mm", dpi = 600)
svglite::svglite(file.path("..", "figures", "Fig2_Jaccard.svg"), width = 210/25.4, height = 297/25.4); print(final_plot); dev.off()

cat(sprintf("\n🎉 [完美组图] 细节全部校准，SCI 级多面板图版已生成！请查看文件: %s\n", out_pdf))