# =============================================================================
# compute_downstream_stats.R
# CORRECT downstream indicators — strictly per Methods 2.7/2.8 — from the 
# (65-sample) family-level RRA matrices & Diet tables built by build_downstream_data.py.
#
# Currency (Methods 2.7): within-sample RRA; mRRA = mean RRA across samples.
# Computed here (all vegan / closed-form, no external network packages needed):
#   * Beta  : Bray-Curtis PCoA, PERMANOVA (adonis2), PERMDISP (betadisper),
#             SIMPER  — all on the family-level RRA matrix (Methods 2.8).
#   * IndVal: Dufrene-Legendre indicator value (specificity x fidelity x 100)
#             + 999-permutation P.  (labdsv/indicspecies equivalent.)
#   * Niche : Levins' standardized niche breadth (Ba, common pool);
#             Pianka O on mRRA + Gotelli RA3 null (1999 perms, seed=1);
#             Bluthgen raw d' and H2' (bipartite-equivalent, closed form).
#   * Alpha : richness/Shannon/Pielou on rarefied raw reads (Counts_rarefied_*,
#             rarefied to min total sequencing depth by build_downstream_data.py).
#
# Run:  "D:/R-4.5.1/bin/Rscript.exe" compute_downstream_stats.R
# =============================================================================
suppressPackageStartupMessages({ library(vegan); library(dplyr) })
set.seed(1)

BASE <- Sys.getenv("MUNTJAC_DATA", unset = ".")  # 指向 Zenodo 数据包根目录
RES  <- file.path(BASE, "03_diet_composition", "rarefied_counts")   # RRA + rarefied counts
DIET <- file.path(BASE, "03_diet_composition")
dir.create(RES, showWarnings = FALSE, recursive = TRUE)

rd <- function(p) tryCatch(read.csv(p, fileEncoding="UTF-8-BOM", check.names=FALSE, stringsAsFactors=FALSE),
                           error=function(e) read.csv(p, check.names=FALSE, stringsAsFactors=FALSE))
wr <- function(df, name) write.csv(df, file.path(RES, name), row.names=FALSE)

SP <- c("M. reevesi","M. crinifrons")

load_rra <- function(thr) {
  d <- rd(file.path(RES, sprintf("RRA_ASM%s.csv", thr)))
  m <- as.matrix(d[, -(1:2)]); rownames(m) <- d$SampleID
  host <- factor(d$Host, levels = SP)
  list(m = m, host = host)
}

# ---- Alpha diversity on rarefied raw reads (Methods 2.8) --------------------
# Methods 2.8: alpha (richness / Shannon / Pielou) on RAW read counts rarefied to
# the minimum per-sample sequencing depth. The rarefaction is performed on the full
# (all-contig) raw read vector to the minimum TOTAL mapped-read depth (build_downstream_data.py,
# Counts_rarefied_*), NOT to the minimum plant depth (which is ~5 reads and degenerate).
# Rarefied plant-family counts -> richness/Shannon/Pielou here.
alpha_diversity <- function(thr) {
  d <- rd(file.path(RES, sprintf("Counts_rarefied_ASM%s.csv", thr)))
  m <- as.matrix(d[, -(1:2)]); rownames(m) <- d$SampleID
  host <- d$Host
  S <- rowSums(m > 0)                            # family richness (rarefied)
  p <- m / rowSums(m)                            # read proportions
  H <- -rowSums(p * log(p), na.rm = TRUE)        # Shannon (natural log); 0*log(0)->NaN removed
  J <- H / log(S); J[S <= 1] <- NA_real_         # Pielou evenness (undefined for S<=1)
  data.frame(SampleID = rownames(m), Host = host, Richness = S,
             Shannon = round(H, 6), Pielou = round(J, 6))
}

# ---- Pianka O + Gotelli RA3 null -------------------------------------------
pianka_ra3 <- function(fr, fb, nperm = 1999) {
  fams <- union(fr$Family, fb$Family)
  px <- setNames(rep(0, length(fams)), fams); py <- px
  px[fr$Family] <- fr$mRRA_pct/100; py[fb$Family] <- fb$mRRA_pct/100
  pk <- function(a,b) sum(a*b)/sqrt(sum(a^2)*sum(b^2))
  O <- pk(px, py)
  ra3 <- function(p){ nz <- which(p>0); p[nz] <- sample(p[nz]); p }  # Gotelli RA3: reshuffle non-zero within row
  set.seed(1)
  nul <- replicate(nperm, pk(ra3(px), ra3(py)))
  data.frame(O=O, E=mean(nul), SES=(O-mean(nul))/sd(nul), p=mean(nul >= O))
}

# ---- Dufrene-Legendre IndVal ----------------------------------------------
indval <- function(m, host, nperm=999) {
  lv <- levels(host)
  out <- lapply(colnames(m), function(f){
    v <- m[, f]
    av <- sapply(lv, function(g) mean(v[host==g]))
    A <- av/sum(av)                                   # specificity
    B <- sapply(lv, function(g) mean(v[host==g] > 0)) # fidelity
    iv <- A*B*100
    k <- which.max(iv)
    obs <- max(iv)
    # permutation: shuffle host labels, recompute max IndVal
    null <- replicate(nperm, {
      h <- sample(host)
      avn <- sapply(lv, function(g) mean(v[h==g])); An <- avn/sum(avn)
      Bn <- sapply(lv, function(g) mean(v[h==g] > 0))
      max(An*Bn*100)
    })
    data.frame(Family=f, IndVal=unname(obs), Indicator_of=lv[k], p_value=mean(null >= obs))
  })
  do.call(rbind, out)
}

# ---- Levins Ba (common pool) -----------------------------------------------
levins_ba <- function(fr, fb) {
  fams <- union(fr$Family, fb$Family)
  px <- setNames(rep(0,length(fams)),fams); py <- px
  px[fr$Family] <- fr$mRRA_pct/100; py[fb$Family] <- fb$mRRA_pct/100
  ba <- function(p){ B <- 1/sum(p^2); (B-1)/(length(p)-1) }
  data.frame(Host=SP, Levins_Ba_commonpool=c(ba(px), ba(py)))
}

# ---- Bluthgen d' & H2' (bipartite-equivalent, closed form) -----------------
# Interaction matrix = mRRA (species x family), the common currency (Methods 2.7).
#   H2' = 1 - H2/H2max ;  H2 = -sum p_ij ln p_ij over realised cells (x ln x -> 0 at 0);
#   H2max = entropy of the marginal-product matrix (summed over ALL cells).
#   d' (raw) = sum_j p_ij ln(p_ij/q_j); d' (std) = raw / -ln(min realised q_j).
bipartite <- function(fr, fb) {
  fams <- union(fr$Family, fb$Family)
  A <- matrix(0, 2, length(fams), dimnames=list(SP, fams))
  A["M. reevesi", fr$Family] <- fr$mRRA_pct/100
  A["M. crinifrons", fb$Family] <- fb$mRRA_pct/100
  p <- A/sum(A)                                     # p_ij (fraction of total)
  H2  <- -sum(p[p>0]*log(p[p>0]))
  pmarg <- outer(rowSums(p), colSums(p))            # q_i * q_j  (all > 0)
  H2max <- -sum(pmarg*log(pmarg))                   # over ALL cells
  H2prime <- 1 - H2/H2max
  q <- colSums(A)/sum(A)                            # partner marginals
  dprime <- sapply(SP, function(sp){
    pp <- A[sp,]/sum(A[sp,])                        # p_ij for species sp (sums to 1)
    d_raw <- sum(pp[pp>0]*log(pp[pp>0]/q[pp>0]))
    d_std <- d_raw / -log(min(q[pp>0]))             # standardised d' (0-1)
    c(raw=d_raw, std=d_std)
  })
  data.frame(H2prime=H2prime,
             dprime_raw_Reeves=unname(dprime["raw","M. reevesi"]),
             dprime_raw_black=unname(dprime["raw","M. crinifrons"]),
             dprime_std_Reeves=unname(dprime["std","M. reevesi"]),
             dprime_std_black=unname(dprime["std","M. crinifrons"]))
}

# =============================================================================
for (thr in c("95","90")) {
  key <- paste0("ASM", thr)
  D <- load_rra(thr)
  m <- D$m; host <- D$host
  fr <- rd(file.path(DIET, sprintf("Diet_ASM%s_Reeves.csv", thr)))
  fb <- rd(file.path(DIET, sprintf("Diet_ASM%s_black.csv", thr)))

  # -- Alpha on rarefied raw reads (Methods 2.8; see alpha_diversity note) ----
  alpha <- alpha_diversity(thr)
  wr(alpha, sprintf("Alpha_%s.csv", thr))
  aR <- alpha[alpha$Host=="M. reevesi",]; aB <- alpha[alpha$Host=="M. crinifrons",]
  cat(sprintf("  Alpha(rarefied raw reads): Reeves richness %.2f Shannon %.3f Pielou %.3f | black richness %.2f Shannon %.3f Pielou %.3f\n",
              mean(aR$Richness), mean(aR$Shannon, na.rm=TRUE), mean(aR$Pielou, na.rm=TRUE),
              mean(aB$Richness), mean(aB$Shannon, na.rm=TRUE), mean(aB$Pielou, na.rm=TRUE)))

  # -- Beta: Bray-Curtis + PERMANOVA + PERMDISP + SIMPER ---------------------
  bc  <- vegdist(m, method="bray")
  per <- adonis2(bc ~ host, permutations=999)
  perm <- data.frame(scenario=key, F=per$F[1], R2=per$R2[1], p=per$`Pr(>F)`[1])
  pd  <- betadisper(bc, host); pdtest <- permutest(pd, permutations=999)
  permd <- data.frame(scenario=key, F=pdtest$statistic, p=pdtest$tab$`Pr(>F)`[1])
  sim <- simper(m, host, permutations=999)
  dsim <- summary(sim)[[1]]
  fams <- rownames(dsim)
  contrib <- c(dsim$cumsum[1], diff(dsim$cumsum)) * 100
  ava <- colMeans(m[host==SP[1],,drop=FALSE]); avb <- colMeans(m[host==SP[2],,drop=FALSE])
  higher <- ifelse(ava[fams] >= avb[fams], SP[1], SP[2])
  simout <- data.frame(Family=fams, higher_in=higher, contrib_pct=unname(contrib),
                       stringsAsFactors=FALSE)
  simout <- simout[order(-simout$contrib_pct), ]

  # -- IndVal ----------------------------------------------------------------
  iv <- indval(m, host)

  # -- Niche -----------------------------------------------------------------
  nb <- levins_ba(fr, fb); nb$scenario <- key
  bp <- bipartite(fr, fb); bp$scenario <- key
  pk <- pianka_ra3(fr, fb); pk$scenario <- key

  wr(perm,  sprintf("PERMANOVA_%s.csv", key))
  wr(permd, sprintf("PERMDISP_%s.csv", key))
  wr(simout, sprintf("SIMPER_%s.csv", key))
  wr(iv,     sprintf("IndVal_%s.csv", key))
  wr(nb,     sprintf("NicheBreadth_Levins_%s.csv", key))
  wr(bp,     sprintf("Bipartite_specialization_%s.csv", key))
  wr(pk,     sprintf("Pianka_%s.csv", key))

  cat(sprintf("== %s () ==============\n", key))
  cat(sprintf("  PERMANOVA  F=%.3f R2=%.4f p=%.3f\n", perm$F, perm$R2, perm$p))
  cat(sprintf("  PERMDISP   F=%.3f p=%.3f\n", permd$F, permd$p))
  cat(sprintf("  Pianka O=%.3f (RA3 p=%.3f, SES=%.3f) | H2'=%.3f | d'std(R)=%.3f d'std(B)=%.3f\n",
              pk$O, pk$p, pk$SES, bp$H2prime, bp$dprime_std_Reeves, bp$dprime_std_black))
  cat(sprintf("  Levins Ba  Reeves=%.3f black=%.3f\n", nb$Levins_Ba_commonpool[1], nb$Levins_Ba_commonpool[2]))
  cat(sprintf("  SIMPER top: %s\n", paste(sprintf("%s(%.1f%%→%s)", simout$Family[1:5], simout$contrib_pct[1:5], simout$higher_in[1:5]), collapse=", ")))
  cat(sprintf("  IndVal top: %s\n", paste(sprintf("%s(%.1f, %s, p=%.3f)", iv$Family[order(-iv$IndVal)][1:5], iv$IndVal[order(-iv$IndVal)][1:5], iv$Indicator_of[order(-iv$IndVal)][1:5], iv$p_value[order(-iv$IndVal)][1:5]), collapse=", ")))
}
cat("== DONE -> ", RES, " ==\n")
