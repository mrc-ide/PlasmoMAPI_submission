
# plot_false_positives.R
#
# Author: Bob Verity
# Date: 2020-09-19
#
# Purpose:
# Read in results of simulations and produce plot summarising the spatial
# distribution of false-positive results.
#
# ------------------------------------------------------------------

library(ggplot2)
library(epitools)

# load PlasmoMAPI package (TODO - replace with instructions to install from github)
#devtools::load_all("/Users/rverity/Dropbox/Bob/Work/My Programs/Barriers to gene flow/PlasmoMAPI")

# function for getting exat binomial confidence intervals from list of extreme
# values
get_binom_exact <- function(w) {
  any_pos <- sum(mapply(length, w) > 0)
  ret <- epitools::binom.exact(any_pos, length(w))
  ret$proportion <- ret$proportion * 100
  ret$lower <- ret$lower * 100
  ret$upper <- ret$upper * 100
  ret
}

# read in simulation results
sim_spatial_random <- readRDS("FigureS1_false_positives_spatial/sim_false_positives_spatial_random.rds")
sim_spatial_regular <- readRDS("FigureS1_false_positives_spatial/sim_false_positives_spatial_regular.rds")
reps <- length(sim_spatial_random$which_upper)
coverage <- do.call(rbind, sim_spatial_random$coverage)
hex_values <- do.call(rbind, sim_spatial_random$hex_values)
hex_values2 <- do.call(rbind, sim_spatial_random$hex_values2)
p <- sim_spatial_random$p
n_hex <- length(p$map$hex)

p2 <- p
p2$output$hex_values <- colMeans(hex_values, na.rm = TRUE)
plot_map(p2, min_hex_coverage = 0)





sim_spatial_random$which_lower <- sim_spatial_random$which_lower
sim_spatial_random$which_upper <- sim_spatial_random$which_upper

# get overall familywise error rate for spatial sims
which_any <- mapply(c, sim_spatial_random$which_lower, sim_spatial_random$which_upper, SIMPLIFY = FALSE)
get_binom_exact(which_any)

which_any <- mapply(c, sim_spatial_regular$which_lower, sim_spatial_regular$which_upper, SIMPLIFY = FALSE)
get_binom_exact(which_any)

# count significant hexes and coverage over all sims
list_upper <- sim_spatial_random$which_upper
tally_upper <- rowSums(mapply(function(x) tabulate(x, nbins = n_hex), list_upper))

list_lower <- sim_spatial_random$which_lower
tally_lower <- rowSums(mapply(function(x) tabulate(x, nbins = n_hex), list_lower))

tally_coverage <- colSums(coverage > 10)


# create base plot
plot_base <- ggplot() + theme_bw()# + theme_bw(base_size = 12)

# plot false positives per hex for random data
p$output$hex_values <- tally_upper / tally_coverage
plot1 <- plot_map(p, min_hex_coverage = 0,
                  plot_significance = FALSE, plot_sampling_points = FALSE,
                  base_plot = plot_base) +
  scale_fill_gradient(low = "black", high = "white", name = "false\npositives (%)") +
  ggtitle("A) random upper")

p$output$hex_values <- tally_lower / tally_coverage
plot2 <- plot_map(p, min_hex_coverage = 0,
                  plot_significance = FALSE, plot_sampling_points = FALSE,
                  base_plot = plot_base) +
  scale_fill_gradient(low = "black", high = "white", name = "false\npositives (%)") +
  ggtitle("B) random lower")


# plot false positives per hex for regular data
p <- sim_spatial_regular$p
n_hex <- length(p$map$hex)

unlist_upper <- unlist(sim_spatial_regular$which_upper)
p$output$hex_values <- tabulate(unlist_upper, nbins = n_hex) / reps * 100
plot3 <- plot_map(p, min_hex_coverage = 0,
                  plot_significance = FALSE, plot_sampling_points = TRUE,
                  point_fill = "green", point_colour = "NA", point_size = 1,
                  base_plot = plot_base) +
  scale_fill_gradient(low = "black", high = "white", name = "false\npositives (%)") +
  ggtitle("C) regular upper")

unlist_lower <- unlist(sim_spatial_regular$which_lower)
p$output$hex_values <- tabulate(unlist_lower, nbins = n_hex) / reps * 100
plot4 <- plot_map(p, min_hex_coverage = 0,
                  plot_significance = FALSE, plot_sampling_points = TRUE,
                  point_fill = "green", point_colour = "NA", point_size = 1,
                  base_plot = plot_base) +
  scale_fill_gradient(low = "black", high = "white", name = "false\npositives (%)") +
  ggtitle("D) regular lower")


# produce combined plot
plot_c <- cowplot::plot_grid(plot1, plot3, plot2, plot4)
plot_c

ggsave("FigureS1_false_positives_spatial/false_positives_spatial.pdf")
ggsave("FigureS1_false_positives_spatial/false_positives_spatial.png", dpi = 100)

