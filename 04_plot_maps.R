# Figure 4: Burn severity maps
# Code developed by M. Buonanduci

# Load packages
library(here)
library(tidyverse)
library(terra)
library(tidyterra)
library(cowplot)
library(ggspatial)


# Load burn severity raster --------
rast_sev <- here("20150814_479980_1161530.tif") %>% rast()

# Calculate uncertainty range (upper - lower)
rast_sev <- rast_sev %>%
  mutate(BA_range = BA_upr - BA_lwr) %>%
  mutate(stems_range = stems_upr - stems_lwr) %>%
  mutate(CC_range = CC_upr - CC_lwr) %>%
  mutate(needle_range = needle_upr - needle_lwr) %>%
  mutate(bchar_circ_range = bchar_circ_upr - bchar_circ_lwr) %>%
  mutate(bchar_ht_range = bchar_ht_upr - bchar_ht_lwr) %>%
  mutate(schar_range = schar_upr - schar_lwr)

# Apply forest mask
rast_sev_forest <- rast_sev * rast_sev$forest


# Plotting function --------

plot_severity_fun <- function(rast, variable, viridis_option, legend_title){
  ggplot() +
    geom_spatraster(data = rast, mapping = aes(fill = {{ variable }})) +
    scale_fill_viridis_c(limits = c(0,1), option = viridis_option, name = legend_title, na.value = "transparent") +
    coord_sf(xlim = xlims, ylim = ylims) + 
    theme_bw() + 
    theme(plot.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "pt"),
          panel.grid = element_blank(), 
          legend.background = element_blank(),
          legend.position = "inside",
          legend.position.inside = c(0.5, 0.4),
          legend.justification = c(0,1),
          legend.key.height = unit(5, "points"),
          legend.key.width = unit(4, "points"),
          legend.title = element_text(size = 7),
          legend.text = element_text(size = 6),
          axis.ticks = element_blank(), 
          axis.text = element_blank())
}



# Figure 4: Example maps -------

# Set ylim
ylims <- c(ext(rast_sev)[3] - 150, ext(rast_sev)[4] + 150)

# Buffer xlim (in meters) to widen plot
xlims <- c(ext(rast_sev)[1] - 200, ext(rast_sev)[2] + 200)

# Basal area mortality 
p11 = plot_severity_fun(rast_sev_forest, BA_lwr, "magma", "Proportion")
p12 = plot_severity_fun(rast_sev_forest, BA, "magma", "Proportion") + theme(legend.position = "none")
p13 = plot_severity_fun(rast_sev_forest, BA_upr, "magma", "Proportion") + theme(legend.position = "none")
p14 = plot_severity_fun(rast_sev_forest, BA_range, "cividis", "Range")
p15 = plot_severity_fun(rast_sev_forest, BA_gt75, "mako", "Probability")

# Dead needle index
p21 = plot_severity_fun(rast_sev_forest, needle_lwr, "magma", "Proportion") + theme(legend.position = "none")
p22 = plot_severity_fun(rast_sev_forest, needle, "magma", "Proportion") + theme(legend.position = "none")
p23 = plot_severity_fun(rast_sev_forest, needle_upr, "magma", "Proportion") + theme(legend.position = "none")
p24 = plot_severity_fun(rast_sev_forest, needle_range, "cividis", "Range") + theme(legend.position = "none")
p25 = plot_severity_fun(rast_sev_forest, needle_gt75, "mako", "Probability") + theme(legend.position = "none")

# Surface char
p31 = plot_severity_fun(rast_sev_forest, schar_lwr, "magma", "Proportion") + theme(legend.position = "none") +
  annotation_scale(location = "br", width_hint = 0.25, line_width = 0.5, style = "ticks", text_cex = 0.5)
p32 = plot_severity_fun(rast_sev_forest, schar, "magma", "Proportion") + theme(legend.position = "none")
p33 = plot_severity_fun(rast_sev_forest, schar_upr, "magma", "Proportion") + theme(legend.position = "none")
p34 = plot_severity_fun(rast_sev_forest, schar_range, "cividis", "Range") + theme(legend.position = "none")
p35 = plot_severity_fun(rast_sev_forest, schar_gt75, "mako", "Probability") + theme(legend.position = "none")

# Create titles
c0 <- ggdraw() + draw_label(" ") 
c1 <- ggdraw() + draw_label("Lower bound\n(quantile 0.05)", 
                            hjust = 0.5, vjust = 0.5, size = 8)
c2 <- ggdraw() + draw_label("Best estimate", 
                            hjust = 0.5, vjust = 0.5, size = 8)
c3 <- ggdraw() + draw_label("Upper bound\n(quantile 0.95)", 
                            hjust = 0.5, vjust = 0.5, size = 8)
c4 <- ggdraw() + draw_label("Uncertainty range\n(upper-lower)", 
                            hjust = 0.5, vjust = 0.5, size = 8)
c5 <- ggdraw() + draw_label("P(severity > 0.75)", 
                            hjust = 0.5, vjust = 0.5, size = 8)

r1 <- ggdraw() + draw_label("Basal area mortality", 
                            hjust = 0.5, vjust = 0.5, angle = 90, size = 8)
r2 <- ggdraw() + draw_label("Dead needle index", 
                            hjust = 0.5, vjust = 0.5, angle = 90, size = 8)
r3 <- ggdraw() + draw_label("Surface char", 
                            hjust = 0.5, vjust = 0.5, angle = 90, size = 8)


plot_grid(c0, c1,  c2,  c3,  c4,  c5,
          r1, p11, p12, p13, p14, p15,
          r2, p21, p22, p23, p24, p25,
          r3, p31, p32, p33, p34, p35,
          nrow = 4,
          rel_widths = c(0.1, 1, 1, 1, 1, 1),
          rel_heights = c(0.2, 1, 1, 1))

