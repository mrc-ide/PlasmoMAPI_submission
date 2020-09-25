
# make_figure1.R
#
# Author: Bob Verity
# Date: 2020-09-17
#
# Purpose:
# Produces Figure1 plot and saves to file.
#
# ------------------------------------------------------------------

library(ggplot2)
library(cowplot)

# load PlasmoMAPI package (TODO - replace with instructions to install from github)
#devtools::load_all("/Users/rverity/Dropbox/Bob/Work/My Programs/Barriers to gene flow/PlasmoMAPI")

set.seed(7)

# define node distribution
L <- 5
n_deme <- 30
node_df <- data.frame(long = runif(n_deme, -L, L),
                      lat = runif(n_deme, -L, L))

# define barriers
barrier_list <- list(data.frame(long = c(-1,1,1,-1,-1)*0.5 + 2, lat = c(1,1,-1,-1,1)*2 + 1),
                     data.frame(long = c(-1,1,1,-1,-1)*0.5 - 2, lat = c(1,1,-1,-1,1)*2 - 1))
nb <- length(barrier_list)

# get distance between nodes, taking barriers into account
barrier_penalty <- 100*c(-1, 1)
sim_dist <- get_barrier_intersect(node_df$long, node_df$lat,
                                  barrier_list = barrier_list,
                                  barrier_penalty = barrier_penalty,
                                  barrier_method = 2,
                                  max_barrier_range = Inf)


# convert to statistical distance
stat_distance <- 5e-3*sim_dist + rnorm(length(sim_dist), sd = 0.5)
stat_distance <- 1.5*exp(-stat_distance)/(1 + exp(-stat_distance))

# make symmetric matrix
stat_distance <- as.matrix(as.dist(stat_distance))

# ------------------------------------------------------------------

# create new PlasmoMAPI project
p <- pm_project()

# load node coordinates
p <- load_coords(p, node_df$long, node_df$lat)

# set up map
p <- create_map(p)

# assign edges to hexes
p <- assign_map(p, eccentricity = 0.95)

# load data
p <- load_data(p, stat_distance, check_delete_output = FALSE)

# run analysis
p <- pm_analysis(p, n_perms = 1e3, n_breaks = 50,
                 min_dist = 0, max_dist = Inf,
                 min_group_size = 5)

# plot coverage
#plot_coverage(p)

# plot edge network
#col_scale <- rev(col_hotcold(20))[c(1,5,10,12,14,16,17,18,19,20)]
#col_scale <- viridisLite::magma(10)
col_scale <- rev(col_hotcold())
plot1 <- plot_network(p, node_size = 2, edge_size = 0.7, zlim = c(0, 1), col_scale = col_scale) +
  ggtitle("A)")

# plot spatial vs. statistical distance
plot2 <- plot_dist(p, col = "black") +
  ggtitle("B)")

# create spinoff project with a subsample of edges shown
set.seed(14)
stat_distance2 <- matrix(0, n_deme, n_deme)
w <- sample(n_deme^2, 20)
stat_distance2[w] <- stat_distance[w]
stat_distance2 <- stat_distance2 + t(stat_distance2)
stat_distance2[stat_distance2 == 0] <- NA
p2 <- load_data(p, stat_distance2, check_delete_output = FALSE)

# plot ellipses from spinoff project
plot3 <- plot_ellipses(p2, node_size = 2, alpha = 0.3, n = 50, eccentricity = 0.95,
                       zlim = c(0, 1), col_scale = col_scale) +
  ggtitle("C)")

# plot final hex map
plot4 <- plot_map(p, min_hex_coverage = 10, poly_list = barrier_list) +
  ggtitle("D)")

# produce combined plot
plot_c <- cowplot::plot_grid(plot1, plot2, plot3, plot4)
plot_c

# save plot to file
ggsave("Figure1_schematic/figure1.pdf", width = 10, height = 7)
ggsave("Figure1_schematic/figure1.png", width = 10, height = 7, dpi = 100)
