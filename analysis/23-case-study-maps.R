# Figure 1: maps of the four case-study surveys, coloured by the regions that
# alternate between years (the case studies in 09-case-study-summary.R):
#   a GOA + BC trawl: GOA and QCS + HS (odd) vs. WCHG + WCVI (even)
#   b Norwegian gillnet-fyke net: North vs. South
#   c BC trawl: QCS + HS (odd years) vs. WCHG + WCVI (even years)
#   d BC longline: North vs. South
#
# The 2 x 2 km (BC) and 2 x 2 nmi (GOA) survey grid cells are drawn as squares
# in their native projections and dissolved into one polygon per survey, so
# the figure stays a small vector PDF. Norway's prediction grid is 50 sites,
# so the set locations are shown as points.

library(dplyr)
library(ggplot2)
library(sf)
library(patchwork)

sf_use_s2(FALSE)

# Okabe-Ito blue and vermillion (colour-blind safe) for the two alternating
# regions, plus a fainter blue for GOA (same odd-year group as QCS + HS)
region_cols <- c(goa = "#A6CBE3", a = "#0072B2", b = "#D55E00")
land_col <- "grey82"
crs_bc <- 3005 # BC Albers
crs_ak <- 3338 # Alaska Albers
# Norwegian coast: an Albers equal-area conic like the other two, so all
# panels share one projection family
crs_no <- "+proj=aea +lat_0=60 +lat_1=62.5 +lat_2=67 +lon_0=10 +datum=WGS84 +units=m +no_defs"

# Grids --------------------------------------------------------------------

# Square cells of side `size` (m) around the grid points in `crs`, dissolved
# by `group`. The small overlap closes floating-point slivers between cells.
dissolve_cells <- function(d, group, crs, size, lon = "lon", lat = "lat") {
  st_as_sf(d, coords = c(lon, lat), crs = 4326) |>
    st_transform(crs) |>
    st_buffer(size / 2 + 5, endCapStyle = "SQUARE") |>
    group_by(.data[[group]]) |>
    summarise(.groups = "drop") |>
    rename(survey = 1)
}

syn <- surveyjoin::dfo_synoptic_grid |>
  dissolve_cells("survey", crs = 32609, size = 2000) |>
  mutate(survey = sub("SYN ", "", survey))

hbll <- bind_rows(
  gfplot::hbll_n_grid$grid |> mutate(survey = "North"),
  gfplot::hbll_s_grid$grid |> mutate(survey = "South")
) |>
  dissolve_cells("survey", crs = 32609, size = 2000, lon = "X", lat = "Y")

goa <- surveyjoin::afsc_grid |>
  filter(survey == "Gulf of Alaska Bottom Trawl Survey") |>
  mutate(survey = "GOA") |>
  dissolve_cells("survey", crs = crs_ak, size = 3704)

norway_sites <- readRDS(here::here("data-raw/garn_ruse_predgrid.rds")) |>
  ungroup() |>
  distinct(site, subarea)
norway_sets <- readRDS(here::here("data-raw/garn_ruse_catch_21sp-2026.rds")) |>
  inner_join(norway_sites, by = "site") |>
  distinct(lon, lat, subarea) |>
  st_as_sf(coords = c("lon", "lat"), crs = 4326) |>
  st_transform(crs_no)

# Coastlines ---------------------------------------------------------------

# Natural Earth high-resolution land, cropped to a lon/lat box before
# projecting
land <- rnaturalearth::ne_countries(scale = 10, returnclass = "sf") |>
  st_geometry() |>
  st_make_valid()

coast <- function(xlim, ylim, crs) {
  box <- st_bbox(c(xmin = xlim[1], xmax = xlim[2], ymin = ylim[1], ymax = ylim[2]), crs = 4326)
  suppressWarnings(st_crop(land, box)) |>
    st_transform(crs)
}

# Boxes well beyond the map limits so the crop edges are never in view
coast_bc <- coast(c(-145, -110), c(43, 62), crs_bc)
coast_ak <- coast(c(-180, -100), c(40, 75), crs_ak)
coast_no <- coast(c(-5, 30), c(55, 75), crs_no)

# Plotting -----------------------------------------------------------------

# Projected limits from a lon/lat box (the box's projected bounding box);
# with no box, the bounding box of `x` with margins `pad` (fractions of the
# width and height: left, right, bottom, top)
lims <- function(xlim, ylim, crs, x = NULL, pad = c(0.03, 0.03, 0.4, 0.12)) {
  if (is.null(xlim)) {
    b <- st_bbox(x)
    w <- b[["xmax"]] - b[["xmin"]]
    h <- b[["ymax"]] - b[["ymin"]]
    return(list(
      x = c(b[["xmin"]] - pad[1] * w, b[["xmax"]] + pad[2] * w),
      y = c(b[["ymin"]] - pad[3] * h, b[["ymax"]] + pad[4] * h)
    ))
  }
  b <- st_bbox(c(xmin = xlim[1], xmax = xlim[2], ymin = ylim[1], ymax = ylim[2]), crs = 4326) |>
    st_as_sfc() |>
    st_segmentize(0.1) |>
    st_transform(crs) |>
    st_bbox()
  list(x = b[c("xmin", "xmax")], y = b[c("ymin", "ymax")])
}

# One polygon per phase; the closing (grow, then shrink) seals the seam where
# the QCS and HS grids meet, which are offset by 1 km
by_phase <- function(x, close = 600) {
  x |>
    group_by(phase) |>
    summarise(.groups = "drop") |>
    st_buffer(close) |>
    st_buffer(-close)
}

# Text labels placed at lon/lat
lab_df <- function(..., crs) {
  d <- tibble::tribble(~label, ~lon, ~lat, ~kind, ...)
  xy <- st_coordinates(st_transform(st_as_sf(d, coords = c("lon", "lat"), crs = 4326), crs))
  mutate(d, X = xy[, 1], Y = xy[, 2])
}

map_panel <- function(regions, land, crs, xlim, ylim, labels, fill_labels, title,
                      legend_pos = c(0.02, 0.08), scale_loc = "bl",
                      arrow_loc = "tr", arrow_pad_x = unit(0.2, "cm"), arrow_pad_y = unit(0.2, "cm"),
                      legend_just = round(legend_pos), breaks_x = waiver(), breaks_y = waiver()) {
  l <- lims(xlim, ylim, crs, x = regions)
  is_points <- all(st_geometry_type(regions) == "POINT")
  # Survey grid cells go under the land (hiding cells that straddle the
  # coast); sampling locations go on top
  g <- if (is_points) {
    ggplot() +
      geom_sf(data = land, fill = land_col, colour = NA) +
      geom_sf(aes(colour = phase), data = regions) +
      scale_colour_manual(values = region_cols, breaks = names(fill_labels), labels = fill_labels, name = NULL)
  } else {
    ggplot() +
      geom_sf(aes(fill = phase), data = regions, colour = NA) +
      geom_sf(data = land, fill = land_col, colour = NA) +
      scale_fill_manual(values = region_cols, breaks = names(fill_labels), labels = fill_labels, name = NULL)
  }
  g +
    ggrepel::geom_text_repel(
      aes(X, Y, label = label),
      data = filter(labels, kind == "survey"), size = 3, fontface = "bold", colour = "grey15",
      bg.colour = "white", bg.r = 0.15, force = 0, force_pull = 0, min.segment.length = Inf
    ) +
    geom_text(
      aes(X, Y, label = label),
      data = filter(labels, kind == "place"), size = 2.6, fontface = "italic", colour = "grey35"
    ) +
    scale_x_continuous(breaks = breaks_x) +
    scale_y_continuous(breaks = breaks_y) +
    ggspatial::annotation_scale(
      location = scale_loc, width_hint = 0.25, text_cex = 0.65,
      bar_cols = c("grey70", "white"), line_col = "grey40", height = unit(0.18, "cm"), pad_x = unit(0.3, "cm"), pad_y = unit(0.25, "cm")
    ) +
    # Grid north; each arrow sits near its projection's central meridian
    # (126 W BC, 154 W Alaska, 10 E Norway), where grid north
    # is true north
    ggspatial::annotation_north_arrow(
      location = arrow_loc, which_north = "grid", style = ggspatial::north_arrow_orienteering(
        fill = c("white", "grey35"), line_col = "grey35", text_col = "grey20", text_size = 6
      ),
      height = unit(0.7, "cm"), width = unit(0.7, "cm"), pad_x = arrow_pad_x, pad_y = arrow_pad_y
    ) +
    coord_sf(crs = crs, xlim = l$x, ylim = l$y, expand = FALSE) +
    ggtitle(title) +
    theme_light() +
    theme(
      axis.title = element_blank(),
      axis.text = element_text(size = rel(0.75)),
      panel.grid = element_line(colour = "grey92"),
      panel.border = element_rect(colour = "grey60", fill = NA),
      plot.title = element_text(size = rel(0.95)),
      legend.position = "inside",
      legend.position.inside = legend_pos,
      legend.justification = legend_just,
      legend.background = element_rect(fill = scales::alpha("white", 0.85), colour = NA),
      legend.key.size = unit(10, "pt"),
      legend.text = element_text(size = rel(0.75)),
      legend.margin = margin(3, 4, 3, 4)
    )
}

bc_xlim <- c(-134.4, -124.6)
bc_ylim <- c(48.2, 54.8)

p_syn <- map_panel(
  syn |>
    mutate(phase = ifelse(survey %in% c("QCS", "HS"), "a", "b")) |>
    by_phase() |>
    st_transform(crs_bc),
  coast_bc, crs_bc, bc_xlim, bc_ylim,
  labels = lab_df(
    "QCS", -129.75, 51.7, "survey",
    "HS", -130.75, 53.4, "survey",
    "WCHG", -133.85, 52.75, "survey",
    "WCVI", -127.2, 48.85, "survey",
    "Haida\nGwaii", -132.25, 53.3, "place",
    "Vancouver\nIsland", -125.75, 49.75, "place",
    "British\nColumbia", -126.3, 53.1, "place",
    crs = crs_bc
  ),
  fill_labels = c(a = "QCS + HS (odd years)", b = "WCHG + WCVI (even years)"),
  title = "BC trawl",
  breaks_x = seq(-134, -124, 2)
)

p_hbll <- map_panel(
  hbll |>
    mutate(phase = ifelse(survey == "North", "a", "b")) |>
    st_transform(crs_bc),
  coast_bc, crs_bc, bc_xlim, bc_ylim,
  labels = lab_df(
    "North", -129.6, 52.8, "survey",
    "South", -127.0, 49.3, "survey",
    "Haida\nGwaii", -132.25, 53.3, "place",
    "Vancouver\nIsland", -125.75, 49.75, "place",
    "British\nColumbia", -126.3, 53.1, "place",
    crs = crs_bc
  ),
  fill_labels = c(a = "North", b = "South"),
  title = "BC longline",
  breaks_x = seq(-134, -124, 2)
)

p_goa <- map_panel(
  bind_rows(st_transform(goa, crs_ak), st_transform(syn, crs_ak)) |>
    mutate(phase = case_when(survey == "GOA" ~ "goa", survey %in% c("QCS", "HS") ~ "a", .default = "b")) |>
    by_phase(),
  coast_ak, crs_ak, xlim = NULL, ylim = NULL,
  labels = lab_df(
    "GOA", -152, 55.6, "survey",
    "Alaska", -155.5, 60.6, "place",
    "British\nColumbia", -127, 55.8, "place",
    crs = crs_ak
  ),
  fill_labels = c(goa = "GOA (odd years)", a = "QCS + HS (odd years)", b = "WCHG + WCVI (even years)"),
  legend_pos = c(0.6, 0.04), legend_just = c(0.5, 0),
  arrow_loc = "bl", arrow_pad_x = unit(3.5, "cm"), arrow_pad_y = unit(0.9, "cm"),
  title = "Gulf of Alaska (GOA) + BC trawl",
  breaks_x = seq(-170, -120, 10)
)

p_no <- map_panel(
  norway_sets |> mutate(phase = ifelse(subarea == "North", "a", "b")),
  coast_no, crs_no, c(4.6, 16.3), c(61.9, 68.3),
  labels = lab_df(
    "North", 10.9, 65.9, "survey",
    "South", 5.7, 63.6, "survey",
    "Norway", 13.0, 64.4, "place",
    crs = crs_no
  ),
  fill_labels = c(a = "North", b = "South"),
  title = "Norwegian gillnet–fyke net",
  legend_pos = c(0.98, 0.1),
  scale_loc = "br",
  arrow_loc = "tl", arrow_pad_x = unit(1.6, "cm"),
  breaks_x = seq(6, 16, 4)
)

# Layout: each map has a fixed aspect ratio, so patchwork would centre it in
# its cell and leave ragged outer edges. Instead, size every panel in inches
# so both rows have the same total panel width (row 2 with c and d at equal
# heights), then save at the layout's natural size; all slack then sits in the
# gutters between panels.
aspect <- function(p) {
  pp <- ggplot_build(p)$layout$panel_params[[1]]
  diff(pp$y_range) / diff(pp$x_range)
}
panel_w <- 6.21 # total panel width per row (in)
r <- vapply(list(p_goa, p_no, p_syn, p_hbll), aspect, numeric(1))
w_top <- panel_w * (1 / r[1:2]) / sum(1 / r[1:2])
w_bot <- panel_w * (1 / r[3:4]) / sum(1 / r[3:4])
size <- function(p, w, r, tag) {
  p + ggh4x::force_panelsizes(cols = unit(w, "in"), rows = unit(w * r, "in")) + labs(tag = tag)
}
row_top <- size(p_goa, w_top[1], r[1], "a") + size(p_no, w_top[2], r[2], "b")
row_bot <- size(p_syn, w_bot[1], r[3], "c") + size(p_hbll, w_bot[2], r[4], "d")
# The outer heights set each row's panel area (the rows' titles and axes are
# added around it), so with them in inches too, the whole gtable is absolute
# and sums to the figure's natural size
g <- wrap_plots(row_top, row_bot, ncol = 1, heights = unit(c(w_top[1] * r[1], w_bot[1] * r[3]), "in"))
gt <- patchwork::patchworkGrob(g)
fig_w <- grid::convertWidth(sum(gt$widths), "in", valueOnly = TRUE)
fig_h <- grid::convertHeight(sum(gt$heights), "in", valueOnly = TRUE)

ggsave(here::here("figs/case-study-maps.pdf"), g, width = fig_w, height = fig_h, device = cairo_pdf)
ggsave(here::here("figs/case-study-maps.png"), g, width = fig_w, height = fig_h, dpi = 300, device = ragg::agg_png)
