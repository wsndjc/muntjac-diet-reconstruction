# =============================================================================
# ms-0602 (定稿-v2)  Muntjac faecal-metagenomics diet — figure package
# Target journal : ECOLOGY AND EVOLUTION (Wiley)
# Backend        : R.  Vector PDF + SVG (editable) + 600 dpi TIFF (RGB).
#
#    valid diet: ASM95 -> 40 Reeves + 6 black; ASM90 -> 43 Reeves + 7 black.
#   (collected : 52 Reeves + 13 black = 65.)
#
# ONE FIGURE PER RESULTS SUBSECTION (renumbered; see old->new map in README):
#   Fig 1  Study area / sampling workflow            (Methods 2.1)
#   Fig 2  Sequencing QC & host removal              (§3.1, Table S2)     [thr-independent]
#   Fig 3  Cross-method Jaccard robustness (A–F)     (§3.2, per species)  [shows all thr]
#   Fig 4  Dietary composition (stacked mRRA + chord)(§3.3)               [ASM90 & ASM95]
#   Fig 5  Dietary diversity & niche breadth         (§3.4)               [ASM90 & ASM95]
#   Fig 6  Dietary differentiation (Venn/PCoA/SIMPER/IndVal/overlap) (§3.5)[ASM90 & ASM95]
#   App S1 25-family detection heatmap (per species) (§3.2 supporting)
#
# E&E figure spec honoured: vector PDF/SVG (resolution-independent) + 600 dpi TIFF;
#   RGB; Okabe–Ito colour-blind-safe palette; single(80mm)/double(166mm) column
#   widths; lowercase (a)(b) panel tags; figures numbered by order of appearance.
#
# Also exports: ASM90 & ASM95 () sample×family matrices + per-metric result
#   tables (tables/), for reproducibility.
#
# Run:  Rscript figures_main.R
# =============================================================================
suppressPackageStartupMessages({
  library(ggplot2); library(patchwork); library(dplyr); library(tidyr)
  library(scales);  library(forcats);  library(ggrepel); library(vegan)
  library(readxl);  library(circlize); library(ggplotify)
})
set.seed(1)

# Paths are relative to the archived data package (Zenodo).
# Set DATA_DIR to the directory containing the downloaded data package.
DATA_DIR <- Sys.getenv("MUNTJAC_DATA", unset = ".")  # 指向 Zenodo 数据包根目录
RES_DIR  <- file.path(DATA_DIR, "03_diet_composition")
SUP_DIR  <- file.path(DATA_DIR, "04_statistics")
ABU_DIR  <- file.path(DATA_DIR, "02_abundance_matrices")
OUT_DIR <- file.path(DATA_DIR, "figures")
TAB_DIR <- file.path(OUT_DIR, "tables")
SUPP_DIR<- file.path(OUT_DIR, "Supplementary")
dir.create(SUPP_DIR, showWarnings=FALSE, recursive=TRUE)
QC_XLSX <- "C:/Users/Administrator/Downloads/Table_S2_QC_statistics_2.xlsx"
dir.create(TAB_DIR, showWarnings = FALSE, recursive = TRUE)
rd <- function(p) tryCatch(read.csv(p, fileEncoding="UTF-8-BOM", check.names=FALSE, stringsAsFactors=FALSE),
                           error=function(e) read.csv(p, check.names=FALSE, stringsAsFactors=FALSE))

# ---- E&E style --------------------------------------------------------------
BASE_FONT <- "Arial"
theme_ee <- function(base_size = 7) {
  theme_classic(base_size=base_size, base_family=BASE_FONT) +
    theme(axis.line=element_line(linewidth=0.35, colour="black"), axis.ticks=element_line(linewidth=0.35, colour="black"),
          axis.text=element_text(colour="black", size=base_size), axis.title=element_text(size=base_size+0.5),
          legend.title=element_text(size=base_size), legend.text=element_text(size=base_size-0.5),
          legend.key.size=unit(3.2,"mm"), strip.text=element_text(size=base_size+0.5, face="bold"),
          strip.background=element_blank(), plot.title=element_text(size=base_size+1, face="bold"),
          plot.subtitle=element_text(size=base_size-0.5, colour="grey30"),
          plot.tag=element_text(size=base_size+2.5, face="bold"), panel.grid=element_blank())
}
theme_set(theme_ee())
SP_COL <- c("M. reevesi"="#0072B2", "M. crinifrons"="#D55E00")           # Okabe–Ito, RGB
SP_COMMON <- c("M. reevesi"="Reeves's muntjac", "M. crinifrons"="black muntjac")
SP_SHORT  <- c("M. reevesi"="Reeves's", "M. crinifrons"="black")
as_sp   <- function(x) factor(x, levels=names(SP_COL))
sp_norm <- function(x) recode(as.character(x), "Reeves"="M. reevesi","black"="M. crinifrons",
                              "M. reevesi"="M. reevesi","M. crinifrons"="M. crinifrons")

# Export helper — vector PDF + SVG + 600 dpi TIFF (E&E combo line-art spec)
save_ee <- function(plot, name, width_mm, height_mm, dpi=600, dir=OUT_DIR) {
  w<-width_mm/25.4; h<-height_mm/25.4; f<-function(e) file.path(dir, paste0(name,e))
  ok<-function(expr) tryCatch({expr;TRUE}, error=function(e){message("   ! ",conditionMessage(e));FALSE})
  ok({svglite::svglite(f(".svg"), width=w, height=h); print(plot); dev.off()})
  ok({grDevices::cairo_pdf(f(".pdf"), width=w, height=h, family=BASE_FONT); print(plot); dev.off()})
  ok({ragg::agg_tiff(f(".tiff"), width=w, height=h, units="in", res=dpi, compression="lzw"); print(plot); dev.off()})
  message("  saved: ", name, "  (", width_mm, " x ", height_mm, " mm, ", dpi, " dpi TIFF + vector)")
}

thr_lab <- function(thr) paste0("BLAST ", sub("ASM","",thr), "% identity (", thr, ")")
d73_lab <- "final dataset"
scen_key<- function(thr) paste0(thr, "")
tok     <- ""
tot_N   <- c("M. reevesi"=52L, "M. crinifrons"=13L)   # collected, 

# =============================================================================
# Shared data loaders
# =============================================================================
load_diet <- function(thr) {
  fr <- rd(file.path(RES_DIR, sprintf("Diet_%s_Reeves.csv", thr)))
  fb <- rd(file.path(RES_DIR, sprintf("Diet_%s_black.csv",  thr)))
  list(fr=fr, fb=fb, nR=fr$S_valid[1], nB=fb$S_valid[1])
}
load_otu <- function(thr) {                      #
  t2 <- sub("ASM","",thr)
  fr <- rd(file.path(RES_DIR, sprintf("Matrix_ASM%s_sample_x_family.csv", t2)))
  fb <- rd(file.path(RES_DIR, sprintf("Matrix_ASM%s_sample_x_family.csv", t2)))
  rownames(fr)<-fr$SampleID; fr$SampleID<-NULL; rownames(fb)<-fb$SampleID; fb$SampleID<-NULL
  fams<-union(colnames(fr), colnames(fb)); fr[setdiff(fams,colnames(fr))]<-0; fb[setdiff(fams,colnames(fb))]<-0
  m<-rbind(fr[,fams,drop=FALSE], fb[,fams,drop=FALSE]); m[is.na(m)]<-0
  host<-c(rep("M. reevesi",nrow(fr)), rep("M. crinifrons",nrow(fb)))
  ok<-rowSums(m)>0; list(m=as.matrix(m[ok,,drop=FALSE]), host=as_sp(host[ok]))
}

# =============================================================================
# FIGURE 1 — sampling & analytical workflow ()
# =============================================================================
build_fig1 <- function() {
  box<-function(x,y,w,h,label,fill,tcol="black",tsize=2.2) list(
    r=data.frame(xmin=x-w/2,xmax=x+w/2,ymin=y-h/2,ymax=y+h/2,fill=fill),
    t=data.frame(x=x,y=y,label=label,tsize=tsize,tcol=tcol))
  B<-list(box(3,9,4.8,1.2,"65 fresh faecal samples\n30 plots · 4 forest types · Dec 2024","grey92"),
          box(3,7.2,4.8,1.0,"Shotgun metagenomics + competitive host mapping","grey92"),
          box(1.4,5.3,2.4,1.2,"M. reevesi\n52 samples","#CDE3F1",tcol="#0072B2"),
          box(4.6,5.3,2.4,1.2,"M. crinifrons\n13 samples","#F6D9C7",tcol="#D55E00"),
          box(1.4,3.2,2.4,1.3,"valid diet\n40 (ASM95) / 43 (ASM90)","#CDE3F1",tcol="#0072B2"),
          box(4.6,3.2,2.4,1.3,"valid diet\n6 (ASM95) / 7 (ASM90)","#F6D9C7",tcol="#D55E00"),
          box(3,1.2,4.8,1.0,"Family-level diet matrices → ecological analyses","grey92"))
  rects<-do.call(rbind, lapply(B,`[[`,"r")); texts<-do.call(rbind, lapply(B,`[[`,"t"))
  seg<-data.frame(x=c(3,3,3,3,1.4,4.6,1.4,4.6),xend=c(3,3,1.4,4.6,1.4,4.6,1.4,4.6),
                  y=c(8.4,6.7,6.7,6.7,4.7,4.7,2.55,2.55),yend=c(7.7,5.9,5.9,5.9,3.85,3.85,1.7,1.7))
  ggplot() +
    geom_segment(data=seg, aes(x,y,xend=xend,yend=yend), linewidth=0.3, colour="grey45") +
    geom_rect(data=rects, aes(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax), fill=rects$fill, colour="grey35", linewidth=0.3) +
    geom_text(data=texts, aes(x,y,label=label), size=texts$tsize, colour=texts$tcol, lineheight=0.95) +
    coord_cartesian(xlim=c(-0.2,6.2), ylim=c(0.4,9.9)) +
    labs(title="Sampling design and analytical workflow",
         subtitle="Data-driven sample-flow schematic (a geographic study-area map must be produced in GIS)") +
    theme_void(base_family=BASE_FONT) + theme(plot.title=element_text(size=8, face="bold"),
      plot.subtitle=element_text(size=6.2, colour="grey30"), plot.margin=margin(6,6,6,6))
}

# =============================================================================
# FIGURE 2 (§3.1) — sequencing QC & host removal (Table S2, )
# =============================================================================
build_fig2 <- function() {
  raw <- read_excel(QC_XLSX, sheet=1, skip=1)
  names(raw) <- c("No","SampleID","Species","Site","RawReads","RawGb","Q20","Q30","HostReads","HostRate")
  qc <- raw |> filter(!is.na(suppressWarnings(as.numeric(No)))) |>
    mutate(across(c(No,RawReads,RawGb,Q20,Q30,HostReads,HostRate), ~suppressWarnings(as.numeric(.))),
           Host=as_sp(recode(Species, "Muntiacus reevesi"="M. reevesi","Muntiacus crinifrons"="M. crinifrons")))
  retained <- rd(file.path(RES_DIR, sprintf("Alpha_perSample_ASM95.csv")))$SampleID
  qc$Retained <- ifelse(qc$SampleID %in% retained, "analysed", "below detection")
  qc_ord <- qc |> arrange(desc(RawGb)) |> mutate(rank=row_number())
  pa <- ggplot(qc_ord, aes(rank, RawGb, fill=Host)) + geom_col(width=0.8) +
    geom_hline(yintercept=median(qc$RawGb), linetype="dashed", linewidth=0.3, colour="grey40") +
    annotate("text", x=45, y=median(qc$RawGb), vjust=-0.5, size=1.9, colour="grey35", label=sprintf("median %.1f Gb", median(qc$RawGb))) +
    scale_fill_manual(values=SP_COL, name=NULL, labels=function(x) SP_COMMON[x]) +
    scale_x_continuous(expand=expansion(mult=c(0.01,0.02))) +
    labs(x="Samples ranked by sequencing output", y="Raw data (Gb)", tag="A") +
    theme(legend.position="top", legend.text=element_text(face="italic"), axis.text.x=element_blank(), axis.ticks.x=element_blank())
  pb <- ggplot(qc, aes(Host, Q30, fill=Host)) +
    geom_boxplot(width=0.55, outlier.shape=NA, linewidth=0.3, alpha=0.55) + geom_jitter(aes(colour=Host), width=0.14, size=0.7, alpha=0.85) +
    scale_fill_manual(values=SP_COL, guide="none") + scale_colour_manual(values=SP_COL, guide="none") +
    scale_x_discrete(labels=function(x) SP_SHORT[x]) + labs(x=NULL, y="Q30 (%)", tag="B") +
    theme(axis.text.x=element_text(face="italic"))
  pc <- ggplot(qc, aes(Host, HostRate, fill=Host)) +
    geom_boxplot(width=0.55, outlier.shape=NA, linewidth=0.3, alpha=0.55) + geom_jitter(aes(colour=Host, shape=Retained), width=0.14, size=0.9, alpha=0.85) +
    geom_hline(yintercept=mean(qc$HostRate), linetype="dashed", linewidth=0.3, colour="grey40") +
    scale_fill_manual(values=SP_COL, guide="none") + scale_colour_manual(values=SP_COL, guide="none") +
    scale_shape_manual(values=c("analysed"=19,"below detection"=1), name=NULL) +
    scale_x_discrete(labels=function(x) SP_SHORT[x]) +
    labs(x=NULL, y="Host alignment rate (%)", tag="C", subtitle="dashed = overall mean host rate") +
    theme(axis.text.x=element_text(face="italic"), legend.position="bottom")
  write.csv(qc, file.path(TAB_DIR, "Table_S2_QC.csv"), row.names=FALSE)
  pa / (pb | pc) + plot_layout(heights=c(1,1.05)) +
    plot_annotation(title="Sequencing quality control and host-read removal",
      subtitle="§3.1 · 65 samples · 52 Reeves's + 13 black muntjac · Q30 > 95.7% · host reads 46.6%",
      theme=theme_ee())
}

# =============================================================================
# FIGURE 3 (§3.2) — per-species cross-method Jaccard (A–F) + App S1 heatmap
# =============================================================================
meth_levels <- c("Traditional","DIR_Loose","DIR_Medium","DIR_Strict","ASM_90","ASM_95","ASM_99")
meth_pretty <- c(Traditional="Trad.", DIR_Loose="DIR-L", DIR_Medium="DIR-M", DIR_Strict="DIR-S", ASM_90="ASM-90", ASM_95="ASM-95", ASM_99="ASM-99")
jacc_mat <- function(csv, subset) {
  m<-rd(csv); rownames(m)<-m$Family_Latin; M<-as.matrix(m[,meth_levels])
  if (subset!="all"){n<-as.integer(sub("top","",subset)); M<-M[order(rowSums(M), decreasing=TRUE)[seq_len(min(n,nrow(M)))],,drop=FALSE]}
  PA<-(M>0)*1; J<-matrix(NA_real_,7,7,dimnames=list(meth_levels,meth_levels))
  for(i in meth_levels)for(j in meth_levels){a<-PA[,i]>0;b<-PA[,j]>0;u<-sum(a|b);J[i,j]<-ifelse(u==0,NA,sum(a&b)/u)}
  J
}
jaccard_tile <- function(csv, subset, sp_name, tag) {
  J<-jacc_mat(csv, subset); df<-as.data.frame(as.table(J)); names(df)<-c("m1","m2","J")
  df$m1<-factor(meth_pretty[as.character(df$m1)], levels=meth_pretty[meth_levels])
  df$m2<-factor(meth_pretty[as.character(df$m2)], levels=rev(meth_pretty[meth_levels]))
  sub_lab<-c(all="all families", top20="top 20 families", top10="top 10 families")[subset]
  ggplot(df, aes(m1, m2, fill=J)) + geom_tile(colour="white", linewidth=0.4) +
    geom_text(aes(label=ifelse(is.na(J),"",sprintf("%.2f",J))), size=1.5, colour=ifelse(!is.na(df$J)&df$J>0.6,"white","grey20")) +
    scale_fill_gradient(low="#F7FBFF", high="#08519C", name="Jaccard", limits=c(0,1), labels=percent, na.value="grey92") +
    coord_equal() + labs(x=NULL, y=NULL, tag=tag, title=sprintf("%s · %s", SP_SHORT[[sp_name]], sub_lab)) +
    theme(axis.text.x=element_text(angle=45, hjust=1, size=5.6), axis.text.y=element_text(size=5.6),
          axis.line=element_blank(), axis.ticks=element_blank(), plot.title=element_text(size=6.8, face="bold"))
}
build_fig3 <- function() {
  fR<-file.path(ABU_DIR,"Matrix_detections_Reeves.csv"); fB<-file.path(ABU_DIR,"Matrix_detections_crinifrons.csv")
  pa<-jaccard_tile(fR,"all","M. reevesi","A");   pb<-jaccard_tile(fB,"all","M. crinifrons","B")
  pc<-jaccard_tile(fR,"top20","M. reevesi","C"); pd<-jaccard_tile(fB,"top20","M. crinifrons","D")
  pe<-jaccard_tile(fR,"top10","M. reevesi","E"); pf<-jaccard_tile(fB,"top10","M. crinifrons","F")
  # export the Jaccard matrices
  for(sp in list(c(fR,"reevesi"),c(fB,"black"))) for(su in c("all","top20","top10"))
    write.csv(round(jacc_mat(sp[1],su),3), file.path(TAB_DIR, sprintf("Jaccard_%s_%s.csv", sp[2], su)))
  (pa|pb)/(pc|pd)/(pe|pf) + plot_layout(guides="collect") +
    plot_annotation(title="Cross-method and cross-threshold agreement of the dietary catalogue",
      subtitle="§3.2 · family-level Jaccard similarity (presence/absence) among the traditional survey, DIR and ASM methods · left M. reevesi, right M. crinifrons",
      theme=theme_ee()) & theme(legend.position="right")
}
heat_species <- function(csv, sp_name, tag, topn=25) {
  m<-rd(csv); rownames(m)<-m$Family_Latin; M<-as.matrix(m[,meth_levels])
  ord<-names(sort(rowSums(M), decreasing=TRUE))[seq_len(min(topn,nrow(M)))]; sub<-M[ord,,drop=FALSE]
  df<-as.data.frame(as.table(sub)); names(df)<-c("Family","Method","count")
  df$Method<-factor(meth_pretty[as.character(df$Method)], levels=meth_pretty[meth_levels]); df$Family<-factor(df$Family, levels=rev(ord)); df$present<-df$count>0
  ggplot(df, aes(Method, Family)) + geom_tile(aes(fill=ifelse(present,count,NA)), colour="grey88", linewidth=0.4) +
    scale_fill_gradient(low="#DEEBF7", high="#08519C", na.value="grey96", name="Detections", trans="log1p", breaks=c(1,2,4,8,16)) +
    scale_x_discrete(expand=c(0,0)) + scale_y_discrete(expand=c(0,0)) +
    labs(x=NULL, y=NULL, tag=tag, title=sprintf("%s (%s)", sp_name, SP_COMMON[[sp_name]])) +
    theme(axis.text.y=element_text(face="italic", size=5.6), axis.text.x=element_text(angle=45, hjust=1, size=6),
          axis.line=element_blank(), axis.ticks=element_blank(), plot.title=element_text(size=7.5, face="bold.italic"))
}
build_appS1 <- function() {
  (heat_species(file.path(ABU_DIR,"Matrix_detections_Reeves.csv"),"M. reevesi","a") |
   heat_species(file.path(ABU_DIR,"Matrix_detections_crinifrons.csv"),"M. crinifrons","b")) + plot_layout(guides="collect") +
    plot_annotation(title="Detection of dominant plant families across methods and thresholds",
      subtitle="§3.2 supporting (Appendix S1) · 25 most-recorded families per species · grey = not detected",
      theme=theme_ee()) & theme(legend.position="right")
}

# =============================================================================
# FIGURE 4 (§3.3) — composition: stacked mRRA (top10 + Others) + chord
# =============================================================================
qual_pal <- function(n) {
  pool <- unique(c(ggsci::pal_d3("category20")(20), ggsci::pal_d3("category20b")(20), ggsci::pal_d3("category20c")(20)))
  is_grey <- vapply(pool, function(h){r<-grDevices::col2rgb(h); (max(r)-min(r))<14}, logical(1))
  pool[!is_grey][seq_len(n)]
}
panel_bars <- function(D, tag="B", topn=10, show_foo=TRUE) {
  fams <- union(D$fr$Family, D$fb$Family); rr <- setNames(rep(0,length(fams)),fams); bb<-rr
  rr[D$fr$Family] <- D$fr$mRRA_pct; bb[D$fb$Family] <- D$fb$mRRA_pct; comb <- rr+bb
  foo_r <- setNames(rep(NA_real_,length(fams)),fams); foo_b<-foo_r
  foo_r[D$fr$Family] <- D$fr$FOO_pct; foo_b[D$fb$Family] <- D$fb$FOO_pct
  top <- names(sort(comb, decreasing=TRUE))[seq_len(min(topn,length(fams)))]
  fam_pal <- setNames(qual_pal(length(top)), top); fam_pal["Others"] <- "grey80"
  seg_one <- function(v, foo, host, xc) {
    top_v<-v[top]; others<-sum(v[setdiff(names(v),top)]); ord<-names(sort(top_v, decreasing=TRUE))
    vals<-c(top_v[ord], Others=others); vals<-vals[vals>0]; cum<-cumsum(vals)
    foov <- c(foo[ord], Others=NA_real_)[names(vals)]
    data.frame(Host=host, Family=names(vals), x=xc, ymax=100-c(0,head(cum,-1)), ymin=100-cum, foo=unname(foov))
  }
  df <- rbind(seg_one(bb,foo_b,"M. crinifrons",1), seg_one(rr,foo_r,"M. reevesi",2))
  df$Family<-factor(df$Family, levels=c(top,"Others"))
  df$ymid <- (df$ymin+df$ymax)/2; df$thick <- df$ymax-df$ymin
  df$foolab <- ifelse(show_foo & !is.na(df$foo) & df$thick>=6, sprintf("%.0f%%", df$foo), "")
  xlab <- c(sprintf("M. crinifrons\n(n=%d)", D$nB), sprintf("M. reevesi\n(n=%d)", D$nR))
  ylab <- if (show_foo) "Relative abundance (mRRA) · in-bar labels = FOO" else "Relative abundance (mRRA)"
  ggplot(df) + geom_rect(aes(xmin=x-0.31, xmax=x+0.31, ymin=ymin, ymax=ymax, fill=Family), colour="white", linewidth=0.12) +
    geom_text(aes(x=x, y=ymid, label=foolab), size=1.45, colour="grey12") +
    scale_fill_manual(values=fam_pal, name="Taxa", breaks=c(top,"Others")) +
    scale_x_continuous(breaks=c(1,2), labels=xlab, limits=c(0.4,2.6)) +
    scale_y_continuous(labels=function(x) paste0(x,"%"), expand=expansion(mult=c(0,0.01))) +
    guides(fill=guide_legend(ncol=1, keyheight=unit(2.9,"mm"))) +
    labs(x=NULL, y=ylab, tag=tag) +
    theme(axis.text.x=element_text(face="italic", size=6), legend.text=element_text(face="italic", size=5.6),
          legend.title=element_text(size=6.2), legend.key.size=unit(2.9,"mm"),
          plot.tag=element_text(size=9, face="bold"))
}
CAT_COL <- c(Shared="grey60", crinifrons=unname(SP_COL["M. crinifrons"]), reevesi=unname(SP_COL["M. reevesi"]))
chord_all <- function(D) {
  fams<-union(D$fr$Family, D$fb$Family); rr<-setNames(rep(0,length(fams)),fams); bb<-rr
  rr[D$fr$Family]<-D$fr$mRRA_pct; bb[D$fb$Family]<-D$fb$mRRA_pct; comb<-rr+bb
  sh<-intersect(D$fr$Family, D$fb$Family); ro<-setdiff(D$fr$Family, D$fb$Family); bo<-setdiff(D$fb$Family, D$fr$Family)
  ord<-function(f) f[order(comb[f], decreasing=TRUE)]; fam_order<-c(ord(sh), ord(ro), ord(bo))
  cat<-setNames(rep("Shared", length(fams)), fams); cat[ro]<-"reevesi"; cat[bo]<-"crinifrons"
  M<-matrix(0,2,length(fam_order), dimnames=list(c("M. reevesi","M. crinifrons"), fam_order))
  M["M. reevesi",]<-rr[fam_order]; M["M. crinifrons",]<-bb[fam_order]
  list(M=M, cat=cat, sh=sh, ro=ro, bo=bo, fam_order=fam_order)
}
draw_chord <- function(X) {
  M<-X$M; cat<-X$cat; fo<-X$fam_order
  colmat<-matrix(CAT_COL[cat[fo]], nrow=2, ncol=length(fo), byrow=TRUE, dimnames=dimnames(M))
  grid.col<-c("M. reevesi"=unname(SP_COL["M. reevesi"]), "M. crinifrons"=unname(SP_COL["M. crinifrons"]), setNames(CAT_COL[cat[fo]], fo))
  sect<-c("M. reevesi", fo, "M. crinifrons")
  # wide inter-family gaps (2.6°) are angular whitespace INDEPENDENT of ribbon width,
  # so even near-zero-mRRA family sectors are pushed ≥2.6° apart -> their radial
  # labels no longer collide. Species gaps kept modest so families get more arc.
  gaps<-setNames(rep(2.6, length(sect)), sect); gaps[c("M. reevesi","M. crinifrons")]<-4
  if(length(X$sh)) gaps[X$fam_order[length(X$sh)]]<-5
  if(length(X$ro)) gaps[X$fam_order[length(X$sh)+length(X$ro)]]<-5
  circos.clear()
  circos.par(start.degree=90, gap.after=gaps, points.overflow.warning=FALSE, canvas.xlim=c(-1.75,1.75), canvas.ylim=c(-1.75,1.75))
  chordDiagram(M, order=sect, grid.col=grid.col, col=colmat, transparency=0.32, annotationTrack="grid",
               preAllocateTracks=list(track.height=mm_h(2.5)), directional=1, direction.type="diffHeight", diffHeight=mm_h(0.8))
  circos.track(track.index=1, panel.fun=function(x,y){ s<-CELL_META$sector.index; isp<-s %in% names(SP_COL)
    circos.text(CELL_META$xcenter, CELL_META$ylim[1], s, facing="clockwise", niceFacing=TRUE, adj=c(0,0.5),
                cex=ifelse(isp,0.68,0.42), font=ifelse(isp,3,1), col=ifelse(isp,"black","grey15")) }, bg.border=NA)
  legend(x=-1.78, y=1.8, legend=c("Shared","M. crinifrons excl.","M. reevesi excl."),
         fill=c(CAT_COL[["Shared"]], CAT_COL[["crinifrons"]], CAT_COL[["reevesi"]]), border=NA, bty="n", cex=0.6, text.font=c(1,3,3))
  circos.clear()
}
panel_chord <- function(D, tag="A") as.ggplot(function() draw_chord(chord_all(D))) +
  labs(tag=tag) + theme(plot.tag=element_text(size=9, face="bold"), plot.margin=margin(2,2,2,2))

# mRRA + FOO: per-species top families, mRRA bars with FOO% labels at bar ends
panel_mrra_foo <- function(thr, tag="B", topn=8) {
  D <- load_diet(thr)
  mk <- function(dd, sp){ dd<-dd[order(-dd$mRRA_pct),][seq_len(min(topn,nrow(dd))),]
    data.frame(Host=sp, Family=dd$Family, mRRA=dd$mRRA_pct, FOO=dd$FOO_pct) }
  df <- rbind(mk(D$fr,"M. reevesi"), mk(D$fb,"M. crinifrons"))
  df$key <- factor(paste(df$Family, df$Host, sep="___"),
                   levels=paste(df$Family, df$Host, sep="___")[order(df$Host, df$mRRA)])
  hl <- setNames(c(sprintf("M. reevesi (n=%d)", D$nR), sprintf("M. crinifrons (n=%d)", D$nB)), names(SP_COL))
  df$Host <- factor(df$Host, levels=names(SP_COL))
  ggplot(df, aes(mRRA, key, fill=Host)) +
    geom_col(width=0.72) +
    geom_text(aes(label=sprintf("%.0f%%", FOO)), hjust=-0.15, size=1.75, colour="grey25") +
    facet_wrap(~Host, scales="free_y", ncol=1, labeller=labeller(Host=hl)) +
    scale_fill_manual(values=SP_COL, guide="none") +
    scale_y_discrete(labels=function(x) sub("___.*","",x)) +
    scale_x_continuous(limits=c(0, max(df$mRRA)*1.22), expand=expansion(mult=c(0,0.02))) +
    labs(x="mRRA (%)   ·   labels = FOO", y=NULL, tag=tag) +
    theme(axis.text.y=element_text(face="italic", size=5.6),
          strip.text=element_text(face="bold.italic", size=6), plot.tag=element_text(size=9, face="bold"))
}
# Two panel-B options for Fig 4 (user chooses which becomes canonical):
#   _stacked : B = stacked mRRA composition (top10 + Others), FOO% shown in-bar
#   _barfoo  : B = per-species mRRA horizontal bars with FOO% end-labels (current)
build_fig4_stacked <- function(thr) {
  D<-load_diet(thr)
  wrap_plots(panel_chord(D,"A"), panel_bars(D,"B",topn=10), panel_venn(thr,"C"),
             design="AB\nAB\nAC", widths=c(1.4, 1))
}
build_fig4_barfoo <- function(thr) {
  D<-load_diet(thr)
  wrap_plots(panel_chord(D,"A"), panel_mrra_foo(thr,"B"), panel_venn(thr,"C"),
             design="AB\nAB\nAC", widths=c(1.4, 1))
}

# =============================================================================
# FIGURE 5 (§3.4) — dietary diversity (alpha) + Levins niche breadth
# =============================================================================
panel_alpha <- function(thr, tag="a") {
  a <- rd(file.path(RES_DIR, sprintf("Alpha_perSample_%s.csv", thr))); a$Host<-as_sp(sp_norm(a$Host))
  nval<-a |> count(Host); nmapv<-setNames(nval$n, as.character(nval$Host))
  long<-a |> pivot_longer(c(Richness,Shannon,Evenness), names_to="metric", values_to="value") |> filter(!is.na(value)) |>
    mutate(metric=factor(metric, levels=c("Richness","Shannon","Evenness"), labels=c("Family richness","Shannon index","Pielou evenness")))
  ann<-long |> group_by(metric) |> summarise(p=tryCatch(wilcox.test(value~Host)$p.value, error=function(e) NA_real_), y=max(value)*1.06, .groups="drop") |>
    mutate(lab=ifelse(is.na(p),"",sprintf("Wilcoxon P = %.3f", p)))
  ggplot(long, aes(Host, value, fill=Host)) +
    geom_boxplot(width=0.55, outlier.shape=NA, linewidth=0.3, alpha=0.55) + geom_jitter(aes(colour=Host), width=0.14, size=0.6, alpha=0.8) +
    geom_text(data=ann, aes(x=1.5, y=y, label=lab), inherit.aes=FALSE, size=1.9, colour="grey20") +
    facet_wrap(~metric, scales="free_y", nrow=1) + scale_fill_manual(values=SP_COL, guide="none") + scale_colour_manual(values=SP_COL, guide="none") +
    scale_x_discrete(labels=function(x) paste0(x,"\n(n=",nmapv[x]," of ",tot_N[x],")")) + scale_y_continuous(expand=expansion(mult=c(0.05,0.13))) +
    labs(x=NULL, y=NULL, tag=tag, title="Within-sample (alpha) diet diversity") + theme(axis.text.x=element_text(face="italic", size=6))
}
panel_levins <- function(thr, tag="b") {
  nb <- rd(file.path(SUP_DIR,sprintf("NicheBreadth_Levins_%s.csv", thr))) |> filter(scenario==scen_key(thr)) |> mutate(Host=as_sp(sp_norm(Host)))
  ggplot(nb, aes(Host, Levins_Ba_commonpool, fill=Host)) +
    geom_col(width=0.6) +
    geom_text(aes(label=sprintf("Ba = %.2f\n(%.1f eff. fam.)", Levins_Ba_commonpool, effective_families)), vjust=-0.3, size=1.9, colour="grey20", lineheight=0.9) +
    scale_fill_manual(values=SP_COL, guide="none") + scale_x_discrete(labels=function(x) SP_SHORT[x]) +
    scale_y_continuous(limits=c(0, max(nb$Levins_Ba_commonpool)*1.30), expand=expansion(mult=c(0,0))) +
    labs(x=NULL, y="Levins' standardised niche breadth (Ba)", tag=tag, title="Population-level niche breadth") +
    theme(axis.text.x=element_text(face="italic", size=6.5))
}
build_figS_div <- function(thr) {              # SUPPLEMENTARY (diversity + niche breadth)
  wrap_plots(panel_alpha(thr,"A"), panel_levins(thr,"B"), widths=c(2.3,1)) +
    plot_annotation(subtitle=paste0(thr_lab(thr), " · ", d73_lab), theme=theme_ee())
}

# =============================================================================
# FIGURE 6 (§3.5) — differentiation: Venn · PCoA · SIMPER · IndVal · overlap
# =============================================================================
panel_venn <- function(thr, tag="a") {
  fr<-unique(rd(file.path(RES_DIR, sprintf("Diet_%s_Reeves.csv", thr)))$Family)
  fb<-unique(rd(file.path(RES_DIR, sprintf("Diet_%s_black.csv",  thr)))$Family)
  ns<-length(intersect(fr,fb)); nR<-length(setdiff(fr,fb)); nB<-length(setdiff(fb,fr)); nt<-length(union(fr,fb))
  cir<-function(cx,r=1.15,n=200){t<-seq(0,2*pi,length.out=n); data.frame(x=cx+r*cos(t), y=r*sin(t))}
  cR<-transform(cir(-0.55),grp="M. reevesi"); cB<-transform(cir(0.55),grp="M. crinifrons")
  lab<-data.frame(x=c(-1.02,0,1.02), y=0, n=c(nR,ns,nB), sub=c(SP_SHORT[["M. reevesi"]],"shared",SP_SHORT[["M. crinifrons"]]))
  ggplot() + geom_polygon(data=cR, aes(x,y,fill=grp), alpha=0.40, colour=NA) + geom_polygon(data=cB, aes(x,y,fill=grp), alpha=0.40, colour=NA) +
    geom_text(data=lab, aes(x,y,label=n), size=3.4, fontface="bold") + geom_text(data=lab, aes(x,y=-0.32,label=sub), size=2.0, colour="grey25") +
    scale_fill_manual(values=SP_COL, guide="none") + coord_equal(clip="off") +
    labs(tag=tag) +
    theme_void(base_family=BASE_FONT) + theme(plot.tag=element_text(size=9, face="bold"), plot.margin=margin(4,4,4,4))
}
panel_pcoa <- function(thr, tag="b") {
  D<-load_otu(thr); m<-D$m; host<-D$host; bc<-vegdist(m,"bray"); per<-adonis2(bc~host, permutations=999)
  pd<-rd(file.path(SUP_DIR,sprintf("PERMDISP_%s.csv", thr))) |> filter(scenario==scen_key(thr))
  pc<-cmdscale(bc,k=2,eig=TRUE); eig<-pmax(pc$eig,0); ax1<-100*eig[1]/sum(eig); ax2<-100*eig[2]/sum(eig)
  sc<-as.data.frame(pc$points); names(sc)<-c("PCoA1","PCoA2"); sc$Host<-host; hulls<-sc |> group_by(Host) |> slice(chull(PCoA1,PCoA2)) |> ungroup()
  lab1<-sprintf("PERMANOVA: F = %.2f, R² = %.1f%%, P = %.3f", per$F[1], 100*per$R2[1], per$`Pr(>F)`[1])
  lab2<-if(nrow(pd)) sprintf("PERMDISP: F = %.3f, P = %.3f", pd$PERMDISP_F[1], pd$PERMDISP_p[1]) else ""
  ggplot(sc, aes(PCoA1, PCoA2, colour=Host, fill=Host)) +
    geom_polygon(data=hulls, aes(group=Host), alpha=0.12, colour=NA) +
    geom_vline(xintercept=0, linewidth=0.25, colour="grey85") + geom_hline(yintercept=0, linewidth=0.25, colour="grey85") + geom_point(size=1.4, alpha=0.9) +
    annotate("text", x=-Inf, y=Inf, hjust=-0.04, vjust=1.4, label=lab1, size=1.9, colour="grey20") +
    annotate("text", x=-Inf, y=Inf, hjust=-0.04, vjust=3.1, label=lab2, size=1.9, colour="grey35") +
    scale_colour_manual(values=SP_COL, name=NULL, labels=function(x) SP_COMMON[x]) + scale_fill_manual(values=SP_COL, guide="none") +
    labs(x=sprintf("PCoA1 (%.1f%%)",ax1), y=sprintf("PCoA2 (%.1f%%)",ax2), tag=tag) +
    theme(legend.position="top", legend.text=element_text(face="italic"), plot.tag=element_text(size=9, face="bold"))
}
panel_simper <- function(thr, tag="c") {
  s<-rd(file.path(SUP_DIR, sprintf("SIMPER_%s.csv", thr))) |> slice_max(contrib_pct, n=10, with_ties=FALSE) |>
    mutate(higher_in=as_sp(sp_norm(higher_in)), signed=ifelse(higher_in=="M. crinifrons",contrib_pct,-contrib_pct), Family=fct_reorder(Family, contrib_pct))
  ggplot(s, aes(signed, Family, fill=higher_in)) + geom_col(width=0.72) + geom_vline(xintercept=0, linewidth=0.3, colour="grey40") +
    scale_fill_manual(values=SP_COL, name="Higher in", labels=function(x) SP_SHORT[x]) + scale_x_continuous(labels=function(x) abs(x)) +
    labs(x="SIMPER contribution (%)", y=NULL, tag=tag) +
    theme(axis.text.y=element_text(face="italic", size=6), legend.position="bottom", plot.tag=element_text(size=9, face="bold"))
}
panel_indval <- function(thr, tag="d") {
  iv<-rd(file.path(SUP_DIR, sprintf("IndVal_%s.csv", thr))) |> mutate(Indicator_of=as_sp(sp_norm(Indicator_of)), sig=p_value<0.05) |>
    slice_max(IndVal, n=10, with_ties=FALSE) |> mutate(Family=fct_reorder(Family, IndVal), plab=ifelse(sig, sprintf("P=%.3f", p_value), ""))
  ggplot(iv, aes(IndVal, Family)) + geom_segment(aes(x=0, xend=IndVal, yend=Family, colour=Indicator_of), linewidth=0.4) +
    geom_point(aes(colour=Indicator_of, shape=sig, size=sig)) + geom_text(aes(label=plab), hjust=-0.15, size=1.8, colour="grey30") +
    scale_colour_manual(values=SP_COL, name="Indicator of", labels=function(x) SP_SHORT[x]) +
    scale_shape_manual(values=c(`FALSE`=21,`TRUE`=19), labels=c("n.s.","P < 0.05"), name=NULL, drop=FALSE) +
    scale_size_manual(values=c(`FALSE`=1.3,`TRUE`=2.1), guide="none", drop=FALSE) +
    scale_x_continuous(limits=c(0, max(iv$IndVal)*1.18), expand=expansion(mult=c(0,0.02))) +
    guides(colour=guide_legend(order=1, override.aes=list(shape=19)), shape=guide_legend(order=2)) +
    labs(x="Indicator value (IndVal)", y=NULL, tag=tag) +
    theme(axis.text.y=element_text(face="italic", size=6), legend.position="bottom", legend.box="vertical", legend.spacing.y=unit(0.4,"mm"), plot.tag=element_text(size=9, face="bold"))
}
pianka_stats <- function(thr, nperm=1999) {
  D<-load_diet(thr); fams<-union(D$fr$Family, D$fb$Family); px<-setNames(rep(0,length(fams)),fams); py<-px
  px[D$fr$Family]<-D$fr$mRRA_pct/100; py[D$fb$Family]<-D$fb$mRRA_pct/100
  pk<-function(a,b) sum(a*b)/sqrt(sum(a^2)*sum(b^2)); O<-pk(px,py); set.seed(1); nul<-replicate(nperm, pk(sample(px), sample(py)))
  list(O=O, E=mean(nul), SES=(O-mean(nul))/sd(nul), p=mean(nul>=O))
}
panel_overlap <- function(thr, tag="e") {
  bp<-rd(file.path(SUP_DIR,sprintf("Bipartite_specialization_%s.csv", thr))) |> filter(scenario==scen_key(thr)); ps<-pianka_stats(thr)
  df<-data.frame(Host=as_sp(c("M. reevesi","M. crinifrons")), dprime=c(bp$dprime_raw_Reeves[1], bp$dprime_raw_black[1]))
  lab <- sprintf("Pianka O = %.3f (RA3 P = %.3f%s) · H2' = %.3f", ps$O, ps$p, ifelse(ps$p<0.05," *"," n.s."), bp$H2prime[1])
  ggplot(df, aes(Host, dprime, fill=Host)) + geom_col(width=0.6) + geom_text(aes(label=sprintf("d' = %.2f", dprime)), vjust=-0.4, size=1.9, colour="grey20") +
    annotate("text", x=1.5, y=0.985, label=lab, size=1.75, colour="grey25") +
    scale_fill_manual(values=SP_COL, guide="none") + scale_x_discrete(labels=function(x) SP_SHORT[x]) + scale_y_continuous(limits=c(0,1), expand=expansion(mult=c(0,0.02))) +
    labs(x=NULL, y="Specialisation d'", tag=tag) +
    theme(axis.text.x=element_text(face="italic", size=6.5), plot.tag=element_text(size=9, face="bold"))
}
build_fig5 <- function(thr) {                  # differentiation (Venn moved to Fig4): A PCoA · B SIMPER · C IndVal · D overlap
  wrap_plots(panel_pcoa(thr,"A"), panel_simper(thr,"B"), panel_indval(thr,"C"), panel_overlap(thr,"D"),
             design="AB\nCD", heights=c(1.15, 1))
}

# =============================================================================
# EXPORT matrices + per-metric result tables (, ASM90 & ASM95)
# =============================================================================
export_tables <- function(thr) {
  D <- load_otu(thr)
  mat <- data.frame(SampleID=rownames(D$m), Host=as.character(D$host), D$m, check.names=FALSE)
  write.csv(mat, file.path(TAB_DIR, sprintf("Matrix_%s_sample_x_family.csv", thr)), row.names=FALSE)
  # composition (mRRA/FOO/POO/wPOO) already per-species
  file.copy(file.path(RES_DIR, sprintf("Diet_%s_Reeves.csv", thr)), file.path(TAB_DIR, sprintf("Composition_%s_Reeves.csv", thr)), overwrite=TRUE)
  file.copy(file.path(RES_DIR, sprintf("Diet_%s_black.csv",  thr)), file.path(TAB_DIR, sprintf("Composition_%s_black.csv",  thr)), overwrite=TRUE)
  # consolidated per-metric summary
  bc<-vegdist(D$m,"bray"); per<-adonis2(bc~D$host, permutations=999)
  pd<-rd(file.path(SUP_DIR,sprintf("PERMDISP_%s.csv", thr))) |> filter(scenario==scen_key(thr))
  nb<-rd(file.path(SUP_DIR,sprintf("NicheBreadth_Levins_%s.csv", thr))) |> filter(scenario==scen_key(thr))
  bp<-rd(file.path(SUP_DIR,sprintf("Bipartite_specialization_%s.csv", thr))) |> filter(scenario==scen_key(thr)); ps<-pianka_stats(thr)
  a<-rd(file.path(RES_DIR, sprintf("Alpha_perSample_%s.csv", thr))); a$Host<-sp_norm(a$Host)
  wilx<-function(mn){v<-a[[mn]];h<-a$Host;k<-!is.na(v);suppressWarnings(wilcox.test(v[k]~h[k])$p.value)}
  summ <- data.frame(
    metric=c("n_valid_Reeves","n_valid_black","PERMANOVA_F","PERMANOVA_R2","PERMANOVA_P",
             "PERMDISP_F","PERMDISP_P","alpha_Richness_WilcoxP","alpha_Shannon_WilcoxP","alpha_Evenness_WilcoxP",
             "Levins_Ba_Reeves","Levins_Ba_black","Pianka_O_mRRA","Pianka_RA3_SES","Pianka_RA3_P","H2prime","dprime_Reeves","dprime_black"),
    value=c(sum(D$host=="M. reevesi"), sum(D$host=="M. crinifrons"), round(per$F[1],3), round(per$R2[1],4), per$`Pr(>F)`[1],
            round(pd$PERMDISP_F[1],3), pd$PERMDISP_p[1], round(wilx("Richness"),3), round(wilx("Shannon"),3), round(wilx("Evenness"),3),
            round(nb$Levins_Ba_commonpool[nb$Host=="M. reevesi"],3), round(nb$Levins_Ba_commonpool[nb$Host=="M. crinifrons"],3),
            round(ps$O,3), round(ps$SES,3), round(ps$p,3), round(bp$H2prime[1],3), round(bp$dprime_raw_Reeves[1],3), round(bp$dprime_raw_black[1],3)))
  write.csv(summ, file.path(TAB_DIR, sprintf("Results_summary_%s.csv", thr)), row.names=FALSE)
  for(m in c("SIMPER","IndVal")) file.copy(file.path(SUP_DIR, sprintf("%s_%s.csv", m, thr)),
                                           file.path(TAB_DIR, sprintf("%s_%s.csv", m, thr)), overwrite=TRUE)
}

# =============================================================================
# RENDER
# =============================================================================
message("== Rendering E&E (ms-0602) figures to ", OUT_DIR, " ==")
# Fig 1 = study-area map (built by the separate script figure1_map_EE.R)
save_ee(build_fig2(), "Fig2_QC_sequencing", 166, 140)                    # §3.1
save_ee(build_fig3(), "Fig3_MethodJaccard_perSpecies", 166, 210)         # §3.2
for (thr in c("ASM95","ASM90")) {
  message("-- ", thr, "  --")
  save_ee(build_fig4_stacked(thr), sprintf("Fig4_Composition_%s_stacked", thr), 183, 150) # A chord · B stacked mRRA(top10+Others)+FOO · C Venn
  save_ee(build_fig4_barfoo(thr),  sprintf("Fig4_Composition_%s_barfoo",  thr), 183, 150) # A chord · B mRRA bars + FOO · C Venn
  save_ee(build_fig5(thr), sprintf("Fig5_Differentiation_%s", thr), 166, 150)  # A PCoA · B SIMPER · C IndVal · D overlap
  save_ee(build_figS_div(thr), sprintf("FigS1_Diversity_%s", thr),  166, 90, dir=SUPP_DIR)
  export_tables(thr)
}
save_ee(build_appS1(), "FigS2_FamilyMethodHeatmap", 166, 150, dir=SUPP_DIR)
# tidy: remove superseded main-folder outputs from the previous layout
for (p in c("Fig5_Diversity","Fig6_Differentiation","AppS1_FamilyMethodHeatmap","Fig1_StudyWorkflow",
            "Fig4_Composition_ASM95","Fig4_Composition_ASM90"))              # superseded unsuffixed Fig4
  for (e in c(".pdf",".svg",".tiff",".png")) unlink(file.path(OUT_DIR, paste0(p, e)))
message("== DONE. Main figs in ", OUT_DIR, " ; supplementary in ", SUPP_DIR, " ; tables in ", TAB_DIR, " ==")
