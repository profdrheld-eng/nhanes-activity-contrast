#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE)

args <- commandArgs(trailingOnly=TRUE)
if(length(args)!=1) stop('Supply exactly one local work directory')
WORK <- normalizePath(args[1],mustWork=TRUE)
OUT <- file.path(WORK,'displays');dir.create(OUT,showWarnings=FALSE)
targets <- file.path(OUT,c('figure_1_two_regime_binned_source.csv',
  'figure_1_two_regime_summary.csv','Figure_1_two_regime_density.png'))
if(any(file.exists(targets))) stop('Refusing to overwrite completed or partial stage outputs')
HIP_2003 <- file.path(WORK,'2003-2004/cohort/analytic_cohort.csv.gz')
HIP_2005 <- file.path(WORK,'2005-2006/cohort/analytic_cohort.csv.gz')
MALE <- "#5AA9E6"
FEMALE <- "#E888A5"
LOWER_FILL <- "#E8F3EF"
HIGHER_FILL <- "#F2ECF8"
BLACK <- "#000000"

as_flag <- function(x) {
  y <- tolower(as.character(x))
  !is.na(y) & y %in% c("true", "1", "yes")
}

read_cycle <- function(path, cycle, regime) {
  d <- read.csv(gzfile(path), check.names = FALSE)
  d <- d[
    as_flag(d$reference_domain) &
      d$sex %in% c("male", "female") &
      is.finite(d$S) & is.finite(d$D) & is.finite(d$WTMEC2YR),
  ]
  data.frame(
    cycle = cycle,
    regime = regime,
    sex = ifelse(d$sex == "male", "Male", "Female"),
    S = d$S,
    D = d$D,
    pooled_weight = d$WTMEC2YR / 2,
    stringsAsFactors = FALSE
  )
}

hip <- rbind(
  read_cycle(HIP_2003, "2003-2004", "Hip, uniaxial counts"),
  read_cycle(HIP_2005, "2005-2006", "Hip, uniaxial counts")
)
wrist_all <- do.call(rbind,lapply(c('2011-2012','2013-2014'),function(cycle)read.csv(gzfile(file.path(WORK,cycle,'cohort/analytic_cohort.csv.gz')),check.names=FALSE)))
wrist <- wrist_all[
  as_flag(wrist_all$reference_domain) &
    wrist_all$sex %in% c("male", "female") &
    is.finite(wrist_all$S) & is.finite(wrist_all$D) &
    is.finite(wrist_all$WTMEC2YR),
]
wrist <- data.frame(
  cycle = wrist$cycle,
  regime = "Wrist, triaxial MIMS",
  sex = ifelse(wrist$sex == "male", "Male", "Female"),
  S = wrist$S,
  D = wrist$D,
  pooled_weight = wrist$WTMEC2YR / 2,
  stringsAsFactors = FALSE
)

stopifnot(nrow(hip) == 5715, nrow(wrist) == 8514)

all_scores <- c(hip$S, hip$D, wrist$S, wrist$D)
bounds <- quantile(all_scores, c(0.003, 0.997), na.rm = TRUE, names = FALSE)
display_limit <- max(abs(bounds))
edges <- seq(-display_limit, display_limit, length.out = 43)

weighted_bins <- function(x, y, w, edges) {
  ix <- findInterval(x, edges, all.inside = TRUE)
  iy <- findInterval(y, edges, all.inside = TRUE)
  keep <- is.finite(x) & is.finite(y) & is.finite(w)
  out <- matrix(0, nrow = length(edges) - 1, ncol = length(edges) - 1)
  # Explicit levels preserve empty coordinate bins; tapply's default factors
  # would compress gaps and shift density mass to different plotted coordinates.
  bins <- seq_len(length(edges) - 1)
  totals <- tapply(w[keep], list(factor(ix[keep], levels = bins),
                               factor(iy[keep], levels = bins)), sum)
  index <- which(!is.na(totals), arr.ind = TRUE)
  out[index] <- totals[index]
  out
}

panel_specs <- list(
  list(data = hip[hip$sex == "Male", ], sex = "Male",
       regime = "Hip, uniaxial counts", years = "2003-2006"),
  list(data = hip[hip$sex == "Female", ], sex = "Female",
       regime = "Hip, uniaxial counts", years = "2003-2006"),
  list(data = wrist[wrist$sex == "Male", ], sex = "Male",
       regime = "Wrist, triaxial MIMS", years = "2011-2014"),
  list(data = wrist[wrist$sex == "Female", ], sex = "Female",
       regime = "Wrist, triaxial MIMS", years = "2011-2014")
)

panels <- lapply(panel_specs, function(spec) {
  totals <- weighted_bins(
    spec$data$S, spec$data$D, spec$data$pooled_weight, edges
  )
  list(
    totals = totals,
    shares = totals / sum(totals),
    spec = spec
  )
})
common_max_share <- max(unlist(lapply(panels, function(x) x$shares)))

bin_frame <- function(panel) {
  z <- panel$totals
  shares <- panel$shares
  index <- expand.grid(
    x_bin = seq_len(nrow(z)),
    y_bin = seq_len(ncol(z))
  )
  index$regime <- panel$spec$regime
  index$years <- panel$spec$years
  index$sex <- panel$spec$sex
  index$x_low <- edges[index$x_bin]
  index$x_high <- edges[index$x_bin + 1]
  index$y_low <- edges[index$y_bin]
  index$y_high <- edges[index$y_bin + 1]
  index$weighted_total <- as.vector(z)
  index$within_panel_weighted_share <- as.vector(shares)
  index$common_scale_intensity <- sqrt(
    index$within_panel_weighted_share / common_max_share
  )
  index
}
bin_source <- do.call(rbind, lapply(panels, bin_frame))
write.csv(
  bin_source,
  file.path(OUT, "figure_1_two_regime_binned_source.csv"),
  row.names = FALSE
)

summary_frame <- function(data, years, regime) {
  data.frame(
    years = years,
    regime = regime,
    n = nrow(data),
    n_male = sum(data$sex == "Male"),
    n_female = sum(data$sex == "Female"),
    pooled_weight_total = sum(data$pooled_weight),
    display_lower_bound = -display_limit,
    display_upper_bound = display_limit,
    stringsAsFactors = FALSE
  )
}
summary_source <- rbind(
  summary_frame(hip, "2003-2006", "Hip, uniaxial counts"),
  summary_frame(wrist, "2011-2014", "Wrist, triaxial MIMS")
)
write.csv(
  summary_source,
  file.path(OUT, "figure_1_two_regime_summary.csv"),
  row.names = FALSE
)

draw_density_panel <- function(panel, color, label) {
  z <- panel$shares
  par(
    mar = c(4.8, 5.0, 2.2, 0.7),
    las = 1,
    col.axis = BLACK,
    col.lab = BLACK,
    font.axis = 2,
    font.lab = 2,
    fg = BLACK
  )
  plot(
    NA,
    xlim = range(edges), ylim = range(edges), asp = 1, bty = "l",
    xaxs = "i", yaxs = "i",
    xlab = "Transformed self-report activity score",
    ylab = "Transformed device activity score",
    cex.lab = 1.08, font.lab = 2, axes = FALSE
  )
  limits <- range(edges)
  plot_limits <- par("usr")
  x_low <- plot_limits[1]
  x_high <- plot_limits[2]
  y_low <- plot_limits[3]
  y_high <- plot_limits[4]
  polygon(
    c(x_low, x_low, y_high, y_low),
    c(y_low, y_high, y_high, y_low),
    col = LOWER_FILL,
    border = NA
  )
  polygon(
    c(y_low, x_high, x_high, y_high),
    c(y_low, y_low, y_high, y_high),
    col = HIGHER_FILL,
    border = NA
  )
  nz <- which(z > 0, arr.ind = TRUE)
  intensity <- sqrt(z[nz] / common_max_share)
  fills <- vapply(
    intensity,
    function(a) adjustcolor(color, alpha.f = 0.18 + 0.72 * a),
    character(1)
  )
  rect(
    edges[nz[, 1]], edges[nz[, 2]],
    edges[nz[, 1] + 1], edges[nz[, 2] + 1],
    col = fills, border = NA
  )
  abline(a = 0, b = 1, lty = 2, lwd = 2.4, col = BLACK)
  axis(1, cex.axis = 1.00, font = 2)
  axis(2, cex.axis = 1.00, font = 2)
  box(bty = "l", lwd = 1.2)
  mtext(label, side = 3, adj = 0, line = 0.45, font = 2, cex = 0.95)
  span <- diff(limits)
  text(
    limits[1] + 0.19 * span, limits[2] - 0.075 * span,
    "Self-report\nlower", cex = 0.84, font = 2, col = BLACK
  )
  text(
    limits[2] - 0.19 * span, limits[1] + 0.075 * span,
    "Self-report\nhigher", cex = 0.84, font = 2, col = BLACK
  )
}

png(
  file.path(OUT, "Figure_1_two_regime_density.png"),
  width = 3300, height = 3000, res = 300, type = "cairo", bg = "white"
)
layout(
  matrix(c(1, 2, 5, 3, 4, 5), nrow = 2, byrow = TRUE),
  widths = c(1, 1, 0.14),
  heights = c(1, 1)
)
draw_density_panel(
  panels[[1]], MALE, "A1  Male | Hip, 2003-2006"
)
draw_density_panel(
  panels[[2]], FEMALE, "A2  Female | Hip, 2003-2006"
)
draw_density_panel(
  panels[[3]], MALE, "B1  Male | Wrist, 2011-2014"
)
draw_density_panel(
  panels[[4]], FEMALE, "B2  Female | Wrist, 2011-2014"
)
par(mar = c(4.2, 0.3, 2.2, 2.7))
plot.new()
levels <- seq(0, 1, length.out = 101)
for (i in seq_len(100)) {
  rect(
    0.15, levels[i], 0.55, levels[i + 1],
    col = adjustcolor(BLACK, alpha.f = 0.18 + 0.72 * levels[i]),
    border = NA
  )
}
axis(
  4, at = c(0, 0.5, 1), labels = c("Low", "Mid", "High"),
  las = 1, tick = FALSE, cex.axis = 0.72, font.axis = 2
)
mtext(
  "Relative\nweighted\ndensity",
  side = 3, line = -0.7, cex = 0.68, font = 2
)
dev.off()

message(
  "Figure 1 completed: hip n=", nrow(hip),
  "; wrist n=", nrow(wrist)
)
