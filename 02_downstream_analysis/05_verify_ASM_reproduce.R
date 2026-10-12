# =============================================================================
# 验证脚本：从 ASM 科水平矩阵复现 mRRA/FOO/wPOO + Alpha/Beta/Niche
# 覆盖 ASM_90 与 ASM_95 两个阈值
# 输入文件（与本脚本放在同一目录，或改 setwd 到该目录）：
#   小麂_独立OTU矩阵_ASM90.csv / _ASM95.csv
#   黑麂_独立OTU矩阵_ASM90.csv / _ASM95.csv
# 依赖：vegan（PERMANOVA/Bray-Curtis/PCoA）。EcoSimR 可选。
# =============================================================================

## ---- 0. 环境 ----
if (!requireNamespace("vegan", quietly = TRUE)) install.packages("vegan")
library(vegan)

setwd(file.path(Sys.getenv("MUNTJAC_DATA", unset="."), "02_abundance_matrices"))
HOST_R <- "M. reevesi"    # 小麂
HOST_B <- "M. crinifrons" # 黑麂

## ---- 1. 读入并合并某一阈值的科水平矩阵 ----
load_matrix <- function(thr) {
  fr <- read.csv(sprintf("小麂_独立OTU矩阵_ASM%s.csv", thr), check.names = FALSE)
  fb <- read.csv(sprintf("黑麂_独立OTU矩阵_ASM%s.csv", thr), check.names = FALSE)
  rownames(fr) <- fr$SampleID; fr$SampleID <- NULL
  rownames(fb) <- fb$SampleID; fb$SampleID <- NULL
  fams <- union(colnames(fr), colnames(fb))
  fr[setdiff(fams, colnames(fr))] <- 0
  fb[setdiff(fams, colnames(fb))] <- 0
  fr <- fr[, fams]; fb <- fb[, fams]
  host <- c(rep(HOST_R, nrow(fr)), rep(HOST_B, nrow(fb)))
  m <- rbind(fr, fb)
  m[is.na(m)] <- 0
  list(tpm = as.matrix(m), host = host)   # tpm: 原始 TPM 丰度（样本×科）
}

## ---- 2. mRRA / FOO / %FOO / POO / wPOO（基于样本内 RRA 比例）----
diet_metrics <- function(tpm, host, which_host) {
  sub <- tpm[host == which_host, , drop = FALSE]
  sub <- sub[rowSums(sub) > 0, , drop = FALSE]   # S = 有有效植物序列的样本
  S <- nrow(sub)
  RRA <- sweep(sub, 1, rowSums(sub), "/")        # 每个样本归一化为比例
  mRRA <- colMeans(RRA) * 100                    # 平均相对丰度 %
  FOO  <- colSums(sub > 0)                        # 出现样本数
  pFOO <- FOO / S * 100                           # 出现频率 %
  POO  <- FOO / sum(FOO) * 100
  wraw <- (mRRA / 100) * FOO
  wPOO <- wraw / sum(wraw) * 100                  # 加权出现频率 %
  out <- data.frame(Family = names(mRRA), mRRA_pct = round(mRRA, 3),
                    FOO_count = FOO, FOO_pct = round(pFOO, 2),
                    POO_pct = round(POO, 2), wPOO_pct = round(wPOO, 3),
                    S_valid = S, row.names = NULL)
  out <- out[out$FOO_count > 0, ]
  out[order(-out$wPOO_pct), ]
}

## ---- 3. Alpha 多样性（样本内 RRA 比例上计算）----
alpha_div <- function(tpm, host) {
  keep <- rowSums(tpm) > 0
  sub <- tpm[keep, , drop = FALSE]; h <- host[keep]
  RRA <- sweep(sub, 1, rowSums(sub), "/")
  rich <- rowSums(sub > 0)
  shannon <- apply(RRA, 1, function(p) { p <- p[p > 0]; -sum(p * log(p)) })  # 自然对数
  even <- ifelse(rich > 1, shannon / log(rich), NA)                          # Pielou J
  data.frame(SampleID = rownames(sub), Host = h,
             Richness = rich, Shannon = shannon, Evenness = even, row.names = NULL)
}

## ---- 4. Beta 多样性：Bray-Curtis 在【原始 TPM】上，PERMANOVA + PCoA ----
beta_div <- function(tpm, host, perm = 999) {
  keep <- rowSums(tpm) > 0
  sub <- tpm[keep, , drop = FALSE]; h <- factor(host[keep])
  D <- vegdist(sub, method = "bray")            # 注意：用原始 TPM，不做样本内归一化
  set.seed(2024)
  ad <- adonis2(D ~ h, permutations = perm)
  pc <- cmdscale(D, k = 2, eig = TRUE)
  pos <- pc$eig[pc$eig > 0]
  list(F = ad$F[1], R2 = ad$R2[1], p = ad$`Pr(>F)`[1],
       PCoA1 = pc$eig[1] / sum(pos) * 100, PCoA2 = pc$eig[2] / sum(pos) * 100)
}

## ---- 5. Pianka 生态位重叠 + Gotelli RA3 零模型 ----
##       宿主利用向量 = 该宿主【池化 TPM 列和】归一化
pianka <- function(a, b) sum(a * b) / sqrt(sum(a^2) * sum(b^2))

niche_overlap <- function(tpm, host, iter = 1000, seed = 2024) {
  util <- function(hh) {
    sub <- tpm[host == hh, , drop = FALSE]
    sub <- sub[rowSums(sub) > 0, , drop = FALSE]
    cs <- colSums(sub); cs / sum(cs)
  }
  u1 <- util(HOST_R); u2 <- util(HOST_B)
  O <- pianka(u1, u2)
  set.seed(seed)
  null <- replicate(iter, pianka(sample(u1), sample(u2)))  # RA3: 各自打乱资源位置
  E <- mean(null); SD <- sd(null)
  list(O = O, E = E, SD = SD, SES = (O - E) / SD,
       p = (sum(null >= O) + 1) / (iter + 1))
}

## ---- 6. 一键运行某阈值 + 某样本集 ----
run_scenario <- function(thr, drop = NULL, tag = "") {
  d <- load_matrix(thr)
  if (!is.null(drop)) {
    idx <- !(rownames(d$tpm) %in% drop)
    d$tpm <- d$tpm[idx, , drop = FALSE]; d$host <- d$host[idx]
  }
  cat(sprintf("\n############ ASM_%s  [%s]  rows=%d (Reeves=%d, black=%d)\n",
              thr, tag, nrow(d$tpm), sum(d$host == HOST_R), sum(d$host == HOST_B)))

  cat("\n--- Diet metrics: Reeves (小麂) top8 ---\n")
  print(head(diet_metrics(d$tpm, d$host, HOST_R), 8), row.names = FALSE)
  cat("\n--- Diet metrics: black (黑麂) top8 ---\n")
  print(head(diet_metrics(d$tpm, d$host, HOST_B), 8), row.names = FALSE)

  a <- alpha_div(d$tpm, d$host)
  cat("\n--- Alpha (mean±SD by host) ---\n")
  for (hh in c(HOST_R, HOST_B)) {
    s <- a[a$Host == hh, ]
    cat(sprintf("  %-14s n=%d: Richness %.2f±%.2f | Shannon %.3f±%.3f | Evenness %.3f±%.3f\n",
        hh, nrow(s), mean(s$Richness), sd(s$Richness),
        mean(s$Shannon), sd(s$Shannon),
        mean(s$Evenness, na.rm = TRUE), sd(s$Evenness, na.rm = TRUE)))
  }
  for (col in c("Richness", "Shannon", "Evenness")) {
    x <- a[a$Host == HOST_R, col]; y <- a[a$Host == HOST_B, col]
    w <- suppressWarnings(wilcox.test(x, y))
    cat(sprintf("  Wilcoxon %-9s W=%.1f, p=%.4f\n", col, w$statistic, w$p.value))
  }

  b <- beta_div(d$tpm, d$host)
  cat(sprintf("\n--- Beta: PERMANOVA pseudo-F=%.3f, R2=%.3f, p=%.4f | PCoA1=%.1f%%, PCoA2=%.1f%%\n",
              b$F, b$R2, b$p, b$PCoA1, b$PCoA2))

  n <- niche_overlap(d$tpm, d$host)
  cat(sprintf("--- Niche: Pianka O=%.3f | RA3 null E=%.3f±%.3f, SES=%.3f, p=%.4f\n",
              n$O, n$E, n$SD, n$SES, n$p))
}

## ---- 7. 跑全部 4 个场景 ----
for (thr in c("90", "95")) {
  run_scenario(thr, NULL,     "FULL 66")
  run_scenario(thr, "final",  "final dataset (65)")
}
