# =============================================================================
# ms-0602  FIGURE 1 — Study-area / sampling map  (Ecology and Evolution)
#   Jiulongshan National Nature Reserve, Zhejiang, China.
#   Red circles = sampling plots, size ∝ number of faecal samples (D73 EXCLUDED).
#   Grayscale = elevation (SRTM via elevatr). Inset = study-area location in China.
#   Photos: (A) black muntjac (M. crinifrons), (B) Reeves's muntjac (M. reevesi).
#   Backend R. Vector PDF + SVG + 600 dpi TIFF (RGB).
#   Run: "E:/R/R-4.5.3/bin/Rscript.exe" figure1_map_EE.R
# =============================================================================
suppressPackageStartupMessages({
  library(readxl); library(dplyr); library(sf); library(terra); library(tidyterra)
  library(elevatr); library(rnaturalearth); library(ggplot2); library(ggspatial)
  library(ggrepel); library(patchwork); library(scales); library(ggnewscale); library(magick)
})
# Rscript's default device is a null PDF that lacks font metrics -> cowplot's eager
# grob conversion throws "failed to find or load PDF CID font". Open a real ragg
# device so text is measured with proper Arial metrics during plot building.
ragg::agg_png(file.path(tempdir(), "fontdev.png"), width=1200, height=1200, res=150)

DATA_DIR <- Sys.getenv("MUNTJAC_DATA", unset=".")
BASE    <- DATA_DIR
OUT_DIR <- file.path(BASE, "SCI_figures_EE")
PHOTO_A <- file.path(DATA_DIR, "01_plant_reference_database", "photos", "black_muntjac.png")   # black muntjac (camera trap, cropped v2)
PHOTO_B <- file.path(DATA_DIR, "01_plant_reference_database", "photos", "reeves_muntjac.png")   # Reeves's muntjac (camera trap, cropped v2)
SHP     <- file.path(DATA_DIR, "gis", "reserve_boundary.shp")   # reserve boundary (full shapefile w/ sidecars)
BASE_FONT <- "Arial"

# ---- sites (remove D73: P29/JLSK28 5 -> 4, total 66 -> 65) -------------------
f <- file.path(DATA_DIR, "gis", "Sampling_sites_inventory.xlsx")
d <- suppressMessages(read_excel(f, skip=3, col_names=FALSE))
names(d) <- c("plot","site","zone","lat","lon","elev","veg","dom","nE","Esamp","nD","Dsamp","ntot")
d <- d |> filter(grepl("^JLSK", site)) |>
  mutate(across(c(lat,lon,elev,ntot), as.numeric),
         has73 = grepl("D73", ifelse(is.na(Dsamp),"",Dsamp)),
         n = ntot - as.integer(has73),
         zone2 = case_when(grepl("Original", zone) ~ "Original",
                           grepl("Expansion", zone) ~ "Expansion", TRUE ~ "Shared"))
message(sprintf("plots=%d  total samples (noD73)=%d  (D73 was in %s)", nrow(d), sum(d$n), d$site[d$has73]))
pts <- st_as_sf(d, coords=c("lon","lat"), crs=4326, remove=FALSE)

# ---- reserve boundary (WGS84; 2 polygons = reserve areas) --------------------
reserve <- st_read(SHP, quiet=TRUE) |>
  st_transform(4326) |> st_make_valid()
res_u <- st_union(reserve)
rb <- st_bbox(reserve); pad <- 0.008
xr <- c(rb[["xmin"]]-pad, rb[["xmax"]]+pad); yr <- c(rb[["ymin"]]-pad, rb[["ymax"]]+pad)

# ---- DEM (SRTM) for reserve extent, masked to the reserve -------------------
DEM_CACHE <- file.path(DATA_DIR, "gis", "jls_dem_z12.tif")
if (file.exists(DEM_CACHE)) {
  message("loading cached DEM…"); dem <- terra::rast(DEM_CACHE)
} else {
  message("downloading DEM (elevatr)…")
  dem_r <- NULL
  for (att in 1:4) {
    dem_r <- tryCatch(get_elev_raster(locations=reserve, z=12, clip="bbox", expand=0.01, verbose=FALSE),
                      error=function(e){message("   retry ", att, ": ", conditionMessage(e)); NULL})
    if (!is.null(dem_r)) break; Sys.sleep(3)
  }
  if (is.null(dem_r)) stop("DEM download failed after retries")
  dem <- terra::rast(dem_r); names(dem) <- "elev"
  terra::writeRaster(dem, DEM_CACHE, overwrite=TRUE)
}
names(dem) <- "elev"
dem <- terra::crop(dem, terra::ext(xr[1], xr[2], yr[1], yr[2]))
dem <- terra::crop(dem, terra::ext(xr[1], xr[2], yr[1], yr[2]))   # full extent (no mask): show surrounding terrain
hs <- terra::shade(terra::terrain(dem, "slope", unit="radians"),
                   terra::terrain(dem, "aspect", unit="radians"), angle=45, direction=315)
names(hs) <- "hillshade"

# ---- MAIN MAP ---------------------------------------------------------------
p_map <- ggplot() +
  geom_spatraster(data=dem, maxcell=6e5) +
  scale_fill_gradient(low="#efedeb", high="#5a5a5a", na.value=NA, name="Elevation\n(m a.s.l.)",
                      guide=guide_colorbar(barwidth=unit(2.6,"mm"), barheight=unit(20,"mm"))) +
  ggnewscale::new_scale_fill() +
  geom_spatraster(data=hs, maxcell=6e5, alpha=0.35, show.legend=FALSE) +
  scale_fill_gradient(low="black", high="white", na.value=NA, guide="none") +
  geom_sf(data=reserve, fill=NA, colour="grey15", linewidth=0.45) +
  geom_point(data=d, aes(lon, lat, size=n), shape=21, fill="red", colour="black",
             stroke=0.35, alpha=0.85) +
  geom_text_repel(data=d, aes(lon, lat, label=site), size=1.9, colour="grey15",
                  segment.size=0.2, min.segment.length=0, box.padding=0.35, max.overlaps=Inf, seed=1) +
  scale_size_area(max_size=7, breaks=c(1,2,5,12,20), name="Faecal\nsamples (n)") +
  annotation_scale(location="bl", width_hint=0.28, height=unit(1.2,"mm"),
                   text_cex=0.55, line_width=0.4, pad_y=unit(1.5,"mm")) +
  annotation_north_arrow(location="tr", height=unit(8,"mm"), width=unit(8,"mm"),
                         style=north_arrow_orienteering(text_size=6, line_width=0.5)) +
  coord_sf(xlim=xr, ylim=yr, expand=FALSE, crs=4326) +
  labs(x="Longitude (°E)", y="Latitude (°N)") +
  theme_bw(base_size=7, base_family=BASE_FONT) +
  theme(panel.grid=element_line(linewidth=0.2, colour="grey88"),
        axis.text=element_text(colour="black", size=6), legend.title=element_text(size=6),
        legend.text=element_text(size=5.5), legend.key.size=unit(3,"mm"),
        legend.background=element_rect(fill=alpha("white",0.7), colour=NA),
        legend.box.just="left", legend.spacing.y=unit(0.5,"mm"))

# (China inset removed per request — figure = reserve sampling map + camera photos only)

# ---- PHOTOS (ggplot + rasterGrob; deferred render avoids the null-device font bug)
photo_panel <- function(path, tag, common, latin) {
  img  <- magick::image_read(path); info <- magick::image_info(img); ar <- info$height/info$width
  grob <- grid::rasterGrob(grDevices::as.raster(img), width=grid::unit(1,"npc"), height=grid::unit(1,"npc"))
  ggplot() +
    annotation_custom(grob, xmin=0, xmax=1, ymin=0, ymax=1) +
    annotate("rect", xmin=0, xmax=0.52, ymin=0.82, ymax=1, fill="black", alpha=0.32) +
    annotate("text", x=0.02, y=0.965, label=sprintf("(%s) %s", tag, common), hjust=0, vjust=1,
             fontface="bold", size=2.7, colour="white") +
    annotate("text", x=0.02, y=0.885, label=latin, hjust=0, vjust=1, fontface="italic", size=2.4, colour="white") +
    scale_x_continuous(limits=c(0,1), expand=c(0,0)) + scale_y_continuous(limits=c(0,1), expand=c(0,0)) +
    coord_fixed(ratio=ar, clip="off") + theme_void() +
    theme(plot.margin=margin(1,1,1,1), plot.background=element_blank())
}
pA <- photo_panel(PHOTO_A, "B", "Black muntjac",   "Muntiacus crinifrons")
pB <- photo_panel(PHOTO_B, "C", "Reeves's muntjac","Muntiacus reevesi")

# equal-height photo row: panel width ∝ 1/aspect (aspect = h/w), so both photos
# render at the same height with widths matched to their (different) crops.
arA <- {i<-magick::image_info(magick::image_read(PHOTO_A)); i$height/i$width}
arB <- {i<-magick::image_info(magick::image_read(PHOTO_B)); i$height/i$width}
photo_row <- (pA | pB) + plot_layout(widths=c(1/arA, 1/arB))

# ---- COMPOSE (patchwork: reserve map over the two camera photos) ------------
# Landscape layout: map (A) on the left, two photos (B, C) stacked on the right
photo_col <- (pB / pA) + plot_layout(heights=c(1, 1))   # C 在上（与 A 顶部对齐），B 在下
fig1 <- (p_map | photo_col) + plot_layout(widths=c(2.2, 1)) +
  plot_annotation(tag_levels="A") &
  theme(plot.tag=element_text(size=12, face="bold"))

# ---- EXPORT (E&E: vector PDF + SVG + 600 dpi TIFF) --------------------------
save_map <- function(plot, name, width_mm, height_mm, dpi=600) {
  w<-width_mm/25.4; h<-height_mm/25.4; f<-function(e) file.path(OUT_DIR, paste0(name,e))
  ok<-function(expr) tryCatch({expr;TRUE}, error=function(e){message("   ! ",conditionMessage(e));FALSE})
  ok({svglite::svglite(f(".svg"), width=w, height=h); print(plot); dev.off()})
  ok({grDevices::cairo_pdf(f(".pdf"), width=w, height=h, family=BASE_FONT); print(plot); dev.off()})
  ok({ragg::agg_tiff(f(".tiff"), width=w, height=h, units="in", res=dpi, compression="lzw"); print(plot); dev.off()})
  ok({ragg::agg_png(f(".png"), width=w, height=h, units="in", res=150); print(plot); dev.off()})  # QA preview
  message("  saved: ", name)
}
save_map(fig1, "Fig1_landscape", 171, 120)
message("== DONE Figure 1 ==")
