# Figures 2 & 3: Marginal effects plots
# Code developed by M. Buonanduci

# Load packages
library(here)
library(tidyverse)
library(posterior)
library(cowplot)

# Load custom functions
source(here("00_zoib_functions.R"))


# Load model fits --------

m_BA <- here("m_west_BA.rds") %>% readRDS()
m_bchar_circ <- here("m_west_bchar_circ.rds") %>% readRDS()
m_bchar_ht <- here("m_west_bchar_ht.rds") %>% readRDS()
m_CC <- here("m_west_CC.rds") %>% readRDS()
m_needle <- here("m_west_needle.rds") %>% readRDS()
m_schar <- here("m_west_schar.rds") %>% readRDS()
m_stems <- here("m_west_stems.rds") %>% readRDS()


# Plotting functions ------

# Plotting function for post-fire sampling interval
plot_interaction_interval <- function(m_exp, m_yrep, m, ylabel = "Proportion"){
  
  ggplot(m_exp, aes(x, y = median(est),
                       ymin = t(quantile(est, 0.05)), ymax = t(quantile(est, 0.95)), 
                       group = level, fill = level)) +
    geom_point(data = m$data, mapping = aes(x, y), color = "gray", alpha = 0.2, size = 0.5, inherit.aes = FALSE) +
    geom_ribbon(data = m_yrep, mapping = aes(x, ymin = t(quantile(est, 0.05)),  ymax = t(quantile(est, 0.95)), 
                                                group = level, fill = level), alpha = 0.15) +
    geom_ribbon(alpha = 0.4) +
    geom_line(aes(color = level)) +
    scale_fill_manual(values = c("#d45e00", "gray50"), labels = c("1", "2")) + 
    scale_color_manual(values = c("#d45e00", "gray50"), labels = c("1", "2")) + 
    ylim(0, 1) +
    ylab(ylabel) + xlab("RdNBR") +
    labs(color = "Sampling interval\n(years post-fire)", fill = "Sampling interval\n(years post-fire)") +
    coord_cartesian(expand = TRUE) +
    theme_bw() +
    theme(panel.grid = element_blank(),
          title = element_text(size = 7),
          legend.title = element_text(size = 8),
          legend.text = element_text(size = 8),
          axis.title = element_text(size = 7),
          axis.text = element_text(size = 6))
}

# Plotting function for marginal effects of other covariates
plot_interaction_covar <- function(m_exp, m_yrep, m, ylabel = "Proportion", legend_title = "Quantile"){
  
  ggplot(m_exp, aes(x, y = median(est),
                    ymin = t(quantile(est, 0.05)), ymax = t(quantile(est, 0.95)), 
                    group = level, fill = level)) +
    geom_ribbon(data = m_yrep, mapping = aes(x, ymin = t(quantile(est, 0.05)),  ymax = t(quantile(est, 0.95)), 
                                             group = level, fill = level), alpha = 0.15) +
    geom_ribbon(alpha = 0.4) +
    geom_line(aes(color = level)) +
    scale_fill_viridis_d(option = 'mako', direction = -1, begin = 0.2, end = 0.8) +
    scale_color_viridis_d(option = 'mako', direction = -1, begin = 0.2, end = 0.8) +
    guides(color = guide_legend(reverse = T),
           fill = guide_legend(reverse = T)) +
    ylim(0, 1) +
    ylab(ylabel) + xlab("RdNBR") +
    labs(color = legend_title, fill = legend_title) +
    coord_cartesian(expand = TRUE) +
    theme_bw() +
    theme(panel.grid = element_blank(),
          legend.background = element_blank(),
          title = element_text(size = 7),
          legend.title = element_text(size = 7, margin = margin(b = 2)),
          legend.text = element_text(size = 7),
          axis.title = element_text(size = 7),
          axis.text = element_text(size = 6))
}


# Figure 2: Marginal effect of sampling interval --------

m_BA_exp = plot_interaction(m_BA$data, m_BA, int_var = "postfire2", int_var_cont = FALSE, type = "exp", return_data = TRUE)
m_BA_yrep = plot_interaction(m_BA$data, m_BA, int_var = "postfire2", int_var_cont = FALSE, type = "yrep", return_data = TRUE)

m_stems_exp = plot_interaction(m_stems$data, m_stems, int_var = "postfire2", int_var_cont = FALSE, type = "exp", return_data = TRUE)
m_stems_yrep = plot_interaction(m_stems$data, m_stems, int_var = "postfire2", int_var_cont = FALSE, type = "yrep", return_data = TRUE)

m_CC_exp = plot_interaction(m_CC$data, m_CC, int_var = "postfire2", int_var_cont = FALSE, type = "exp", return_data = TRUE)
m_CC_yrep = plot_interaction(m_CC$data, m_CC, int_var = "postfire2", int_var_cont = FALSE, type = "yrep", return_data = TRUE)

m_needle_exp = plot_post_med(m_needle$data, m_needle, type = "exp", return_data = TRUE) %>% mutate(level = "0")
m_needle_yrep = plot_post_med(m_needle$data, m_needle, type = "yrep", return_data = TRUE) %>% mutate(level = "0")

m_bchar_circ_exp = plot_post_med(m_bchar_circ$data, m_bchar_circ, type = "exp", return_data = TRUE) %>% mutate(level = "0")
m_bchar_circ_yrep = plot_post_med(m_bchar_circ$data, m_bchar_circ, type = "yrep", return_data = TRUE) %>% mutate(level = "0")

m_bchar_ht_exp = plot_post_med(m_bchar_ht$data, m_bchar_ht, type = "exp", return_data = TRUE) %>% mutate(level = "0")
m_bchar_ht_yrep = plot_post_med(m_bchar_ht$data, m_bchar_ht, type = "yrep", return_data = TRUE) %>% mutate(level = "0")

m_schar_exp = plot_interaction(m_schar$data, m_schar, int_var = "postfire2", int_var_cont = FALSE, type = "exp", return_data = TRUE)
m_schar_yrep = plot_interaction(m_schar$data, m_schar, int_var = "postfire2", int_var_cont = FALSE, type = "yrep", return_data = TRUE)


p1 = plot_interaction_interval(m_BA_exp, m_BA_yrep, m_BA) + 
  xlim(range(m_BA$data$x)) + # make x-axis range consistent across panels
  ggtitle("Basal area mortality") +
  theme(legend.position = "none",
        axis.title.x = element_blank())
p2 = plot_interaction_interval(m_stems_exp, m_stems_yrep, m_stems) + 
  xlim(range(m_BA$data$x)) + # make x-axis range consistent across panels
  ggtitle("Stem mortality") +
  theme(legend.position = "none",
        axis.title.x = element_blank(),
        axis.title.y = element_blank())
p3 = plot_interaction_interval(m_CC_exp, m_CC_yrep, m_CC) + 
  xlim(range(m_BA$data$x)) + # make x-axis range consistent across panels
  ggtitle("Canopy cover mortality") +
  theme(legend.position = "none",
        axis.title.x = element_blank(),
        axis.title.y = element_blank())
p4 = plot_interaction_interval(m_needle_exp, m_needle_yrep, m_needle) + 
  xlim(range(m_BA$data$x)) + # make x-axis range consistent across panels
  ggtitle("Dead needle index") +
  theme(legend.position = "none")
p5 = plot_interaction_interval(m_bchar_circ_exp, m_bchar_circ_yrep, m_bchar_circ) + 
  xlim(range(m_BA$data$x)) + # make x-axis range consistent across panels
  ggtitle("Bole char (circumf.)") +
  theme(legend.position = "none",
        axis.title.y = element_blank())
p6 = plot_interaction_interval(m_bchar_ht_exp, m_bchar_ht_yrep, m_bchar_ht) + 
  xlim(range(m_BA$data$x)) + # make x-axis range consistent across panels
  ggtitle("Bole char (height)") +
  theme(legend.position = "none",
        axis.title.y = element_blank())
p7 = plot_interaction_interval(m_schar_exp, m_schar_yrep, m_schar) + 
  xlim(range(m_BA$data$x)) + # make x-axis range consistent across panels
  ggtitle("Ground surface char") +
  theme(legend.position = "none",
        axis.title.y = element_blank())

plegend = get_legend( plot_interaction_interval(m_schar_exp, m_schar_yrep, m_schar) + 
                        theme(legend.box.margin = margin(0, 20, 0, 0),
                              legend.key.size = unit(10, "point")) )

plot_grid(p1, p2, p3, plegend, p4, p5, p6, p7, 
          labels = c("(a)", "(b)", "(c)", "", "(d)", "(e)", "(f)", "(g)"), label_size = 10,
          nrow = 2, rel_widths = c(1.07, 1, 1, 1))


# Figure 3: Marginal effects of other covariates --------

m_BA_tcov_exp = plot_interaction(m_BA$data, m_BA, int_var = "rap_tcov", 
                                 .quant = c(0.95, 0.05), int_var_cont = TRUE, type = "exp", return_data = TRUE)
m_BA_tcov_yrep = plot_interaction(m_BA$data, m_BA, int_var = "rap_tcov", 
                                  .quant = c(0.95, 0.05), int_var_cont = TRUE, type = "yrep", return_data = TRUE)

m_BA_cwd_exp = plot_interaction(m_BA$data, m_BA, int_var = "tt_cwd", 
                                .quant = c(0.95, 0.05), int_var_cont = TRUE, type = "exp", return_data = TRUE)
m_BA_cwd_yrep = plot_interaction(m_BA$data, m_BA, int_var = "tt_cwd", 
                                 .quant = c(0.95, 0.05), int_var_cont = TRUE, type = "yrep", return_data = TRUE)

m_CC_tcov_exp = plot_interaction(m_CC$data, m_CC, int_var = "rap_tcov", 
                                 .quant = c(0.95, 0.05), int_var_cont = TRUE, type = "exp", return_data = TRUE)
m_CC_tcov_yrep = plot_interaction(m_CC$data, m_CC, int_var = "rap_tcov", 
                                  .quant = c(0.95, 0.05), int_var_cont = TRUE, type = "yrep", return_data = TRUE)

m_needle_cwd_exp = plot_interaction(m_needle$data, m_needle, int_var = "tt_cwd", 
                                    .quant = c(0.95, 0.05), int_var_cont = TRUE, type = "exp", return_data = TRUE)
m_needle_cwd_yrep = plot_interaction(m_needle$data, m_needle, int_var = "tt_cwd", 
                                     .quant = c(0.95, 0.05), int_var_cont = TRUE, type = "yrep", return_data = TRUE)

m_needle_slp_exp = plot_interaction(m_needle$data, m_needle, int_var = "lf_slp", 
                                    .quant = c(0.95, 0.05), int_var_cont = TRUE, type = "exp", return_data = TRUE)
m_needle_slp_yrep = plot_interaction(m_needle$data, m_needle, int_var = "lf_slp", 
                                     .quant = c(0.95, 0.05), int_var_cont = TRUE, type = "yrep", return_data = TRUE)

m_bchar_circ_cwd_exp = plot_interaction(m_bchar_circ$data, m_bchar_circ, int_var = "tt_cwd", 
                                        .quant = c(0.95, 0.05), int_var_cont = TRUE, type = "exp", return_data = TRUE)
m_bchar_circ_cwd_yrep = plot_interaction(m_bchar_circ$data, m_bchar_circ, int_var = "tt_cwd", 
                                         .quant = c(0.95, 0.05), int_var_cont = TRUE, type = "yrep", return_data = TRUE)

m_bchar_ht_psr_exp = plot_interaction(m_bchar_ht$data, m_bchar_ht, int_var = "lf_psr", 
                                      .quant = c(0.95, 0.05), int_var_cont = TRUE, type = "exp", return_data = TRUE)
m_bchar_ht_psr_yrep = plot_interaction(m_bchar_ht$data, m_bchar_ht, int_var = "lf_psr", 
                                       .quant = c(0.95, 0.05), int_var_cont = TRUE, type = "yrep", return_data = TRUE)

m_schar_cwd_exp = plot_interaction(m_schar$data, m_schar, int_var = "tt_cwd", 
                                   .quant = c(0.95, 0.05), int_var_cont = TRUE, type = "exp", return_data = TRUE)
m_schar_cwd_yrep = plot_interaction(m_schar$data, m_schar, int_var = "tt_cwd", 
                                    .quant = c(0.95, 0.05), int_var_cont = TRUE, type = "yrep", return_data = TRUE)

m_schar_slp_exp = plot_interaction(m_schar$data, m_schar, int_var = "lf_slp", 
                                   .quant = c(0.95, 0.05), int_var_cont = TRUE, type = "exp", return_data = TRUE)
m_schar_slp_yrep = plot_interaction(m_schar$data, m_schar, int_var = "lf_slp", 
                                    .quant = c(0.95, 0.05), int_var_cont = TRUE, type = "yrep", return_data = TRUE)

m_schar_lat_exp = plot_interaction(m_schar$data, m_schar, int_var = "latitude_wgs84", 
                                   .quant = c(0.95, 0.05), int_var_cont = TRUE, type = "exp", return_data = TRUE)
m_schar_lat_yrep = plot_interaction(m_schar$data, m_schar, int_var = "latitude_wgs84", 
                                    .quant = c(0.95, 0.05), int_var_cont = TRUE, type = "yrep", return_data = TRUE)

p1 = plot_interaction_covar(m_BA_tcov_exp, m_BA_tcov_yrep, m_BA, legend_title = "Pre-fire\ntree cover") + 
  ggtitle("Basal area\nmortality") + 
  ylab("Proportion") + 
  xlim(c(-300, 1400)) +
  theme(legend.position = "inside",
        legend.position.inside = c(0.45, 0.46),
        legend.justification = c(0,1),
        legend.key.size = unit(8, "point"))

p2 = plot_interaction_covar(m_BA_cwd_exp, m_BA_cwd_yrep, m_BA, legend_title = "Climatic\nwater deficit") + 
  ggtitle("Basal area\nmortality") + 
  ylab("Proportion basal area mortality") +
  xlim(c(-300, 1400)) +
  theme(legend.position = "inside",
        legend.position.inside = c(0.4, 0.46),
        legend.justification = c(0,1),
        legend.key.size = unit(8, "point"),
        axis.text.y = element_blank(),
        axis.title.y = element_blank())

p3 = plot_interaction_covar(m_CC_tcov_exp, m_CC_tcov_yrep, m_CC, legend_title = "Pre-fire\ntree cover") + 
  ggtitle("Canopy cover\nmortality") + 
  ylab("Proportion canopy cover mortality") + 
  xlim(c(-300, 1400)) +
  theme(legend.position = "inside",
        legend.position.inside = c(0.45, 0.46),
        legend.justification = c(0,1),
        legend.key.size = unit(8, "point"),
        axis.text.y = element_blank(),
        axis.title.y = element_blank())

p4 = plot_interaction_covar(m_needle_cwd_exp, m_needle_cwd_yrep, m_needle, legend_title = "Climatic\nwater deficit") + 
  ggtitle("Dead needle\nindex") + 
  ylab("Proportional dead needle index") + 
  xlim(c(-300, 1400)) +
  theme(legend.position = "inside",
        legend.position.inside = c(0.4, 0.46),
        legend.justification = c(0,1),
        legend.key.size = unit(8, "point"),
        axis.text.y = element_blank(),
        axis.title.y = element_blank())

p5 = plot_interaction_covar(m_needle_slp_exp, m_needle_slp_yrep, m_needle, legend_title = "Slope\nsteepness") + 
  ggtitle("Dead needle\nindex") + 
  ylab("Proportional dead needle index") + 
  xlim(c(-300, 1400)) +
  theme(legend.position = "inside",
        legend.position.inside = c(0.45, 0.46),
        legend.justification = c(0,1),
        legend.key.size = unit(8, "point"),
        axis.text.y = element_blank(),
        axis.title.y = element_blank())

p6 = plot_interaction_covar(m_bchar_circ_cwd_exp, m_bchar_circ_cwd_yrep, m_bchar_circ, legend_title = "Climatic\nwater deficit") + 
  ggtitle("Bole char\n(circumference)") + 
  ylab("Proportion") + 
  xlim(c(-300, 1400)) +
  theme(legend.position = "inside",
        legend.position.inside = c(0.35, 0.46),
        legend.justification = c(0,1),
        legend.key.size = unit(8, "point"))

p7 = plot_interaction_covar(m_bchar_ht_psr_exp, m_bchar_ht_psr_yrep, m_bchar_circ, legend_title = "Hist. percent\nstand-replacing") + 
  ggtitle("Bole char\n(height)") + 
  ylab("Proportion bole char (height)") + 
  xlim(c(-300, 1400)) +
  theme(legend.position = "inside",
        legend.position.inside = c(0.01, 0.99),
        legend.justification = c(0,1),
        legend.key.size = unit(8, "point"),
        axis.text.y = element_blank(),
        axis.title.y = element_blank())

p8 = plot_interaction_covar(m_schar_cwd_exp, m_schar_cwd_yrep, m_schar, legend_title = "Climatic\nwater deficit") + 
  ggtitle("Ground surface\nchar") + 
  ylab("Proportion ground surface char") + 
  xlim(c(-300, 1400)) +
  theme(legend.position = "inside",
        legend.position.inside = c(0.01, 0.99),
        legend.justification = c(0,1),
        legend.key.size = unit(8, "point"),
        axis.text.y = element_blank(),
        axis.title.y = element_blank())

p9 = plot_interaction_covar(m_schar_slp_exp, m_schar_slp_yrep, m_schar, legend_title = "Slope\nsteepness") + 
  ggtitle("Ground surface\nchar") + 
  ylab("Proportion ground surface char") +
  xlim(c(-300, 1400)) +
  theme(legend.position = "inside",
        legend.position.inside = c(0.01, 0.99),
        legend.justification = c(0,1),
        legend.key.size = unit(8, "point"),
        axis.text.y = element_blank(),
        axis.title.y = element_blank())

p10 = plot_interaction_covar(m_schar_lat_exp, m_schar_lat_yrep, m_schar, legend_title = "Latitude") + 
  ggtitle("Ground surface\nchar") + 
  ylab("Proportion ground surface char") +
  xlim(c(-300, 1400)) +
  theme(legend.position = "inside",
        legend.position.inside = c(0.01, 0.99),
        legend.justification = c(0,1),
        legend.key.size = unit(8, "point"),
        axis.text.y = element_blank(),
        axis.title.y = element_blank())

# Simple layout
plot_grid(p1, p2, p3, p4, p5, 
          p6, p7, p8, p9, p10,
          labels = c("(a)", "(b)", "(c)", "(d)", "(e)", 
                     "(f)", "(g)", "(h)", "(i)", "(j)"), 
          label_x = c(0.1, -0.15, -0.15, -0.15, -0.15,
                      0.1, -0.15, -0.15, -0.15, -0.15),
          label_size = 10, nrow = 2, rel_widths = c(1.2, 1, 1, 1, 1))


