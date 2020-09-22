
# plot_bias.R
#
# Author: Bob Verity
# Date: 2020-09-22
#
# Purpose:
# Reads in the results of sim_bias.R and produces grid of map plots.
#
# ------------------------------------------------------------------

# load PlasmoMAPI package (TODO - replace with instructions to install from github)
#devtools::load_all("/Users/rverity/Dropbox/Bob/Work/My Programs/Barriers to gene flow/PlasmoMAPI")

# load simulation results from file
x <- readRDS("FigureS3_bias/sim_bias.rds")
list_sim <- x$list_sim
list_p <- x$list_p
node_vec <- node_vec
ecc_vec <- ecc_vec
barrier_list <- x$barrier_list


# tally up significant hexes onto project maps
proj_list <- list()
for (i_arrange in seq_along(list_sim)) {
  proj_list[[i_arrange]] <- list()
  for (i_nodes in seq_along(list_sim[[i_arrange]])) {
    proj_list[[i_arrange]][[i_nodes]] <- list()
    for (i_ecc in seq_along(list_sim[[i_arrange]][[i_nodes]])) {
      
      # extract project relevant to these simulations
      p <- list_p[[i_arrange]][[i_nodes]]
      n_hex <- length(p$map$hex)
      
      # store project for producing upper and lower significance maps
      proj_list[[i_arrange]][[i_nodes]][[i_ecc]] <- p
      
      # extract z-score
      z <- rowMeans(mapply(function(x) x$z, list_sim[[i_arrange]][[i_nodes]][[i_ecc]]), na.rm = TRUE)
      
      # create plotting object
      proj_list[[i_arrange]][[i_nodes]][[i_ecc]]$output$hex_values <- z
      
    }
  }
}

# produce list of plots
plot_list <- list()
for (i_arrange in seq_along(proj_list)) {
  plot_nodes <- list()
  
  # parameters for this arrangement
  plot_points <- (i_arrange == 2)
  
  # create list of plots
  for (i_nodes in seq_along(proj_list[[i_arrange]])) {
    plot_ecc <- list()
    
    for (i_ecc in seq_along(proj_list[[i_arrange]][[i_nodes]])) {
      
      # store map plots
      plot_ecc[[i_ecc]] <- plot_map(proj_list[[i_arrange]][[i_nodes]][[i_ecc]],
                                    min_hex_coverage = 0, poly_list = barrier_list,
                                    plot_significance = FALSE, plot_sampling_points = plot_points,
                                    point_size = 0.5)
      
    }
    plot_nodes[[i_nodes]] <- cowplot::plot_grid(plotlist = plot_ecc, nrow = 1)
  }
  
  # save final plot
  plot_list[[i_arrange]] <- cowplot::plot_grid(plotlist = plot_nodes, ncol = 1)
}
plot_c <- cowplot::plot_grid(plotlist = plot_list)
plot_c

# ----------------------------------------------------------------
# add row and column titles

# parameters
heading_size <- 6
subheading_size <- 5

# create eccentricity names as plot objects
ecc_names <- list()
for (i in seq_along(ecc_vec)) {
  ecc_names[[i]] <- ggplot() + theme_void() + annotate("text", x = 0, y = 0, size = subheading_size,
                                                       label = sprintf("eccentricity = %s", ecc_vec[i]))
}

# add eccentricity names
tmp <- list(cowplot::plot_grid(plotlist = rep(ecc_names, 2), nrow = 1),
            plot_c)
plot_c <- cowplot::plot_grid(plotlist = tmp, ncol = 1, rel_heights = c(1, 15))

# add arrangement names
arrange_names <- list()
arrange_names[[1]] <- ggplot() + theme_void() + annotate("text", x = 0, y = 0, size = heading_size,
                                                         label = "random node arrangement")
arrange_names[[2]] <- ggplot() + theme_void() + annotate("text", x = 0, y = 0, size = heading_size,
                                                         label = "regular node arrangement")
tmp <- list(cowplot::plot_grid(plotlist = arrange_names, nrow = 1),
            plot_c)
plot_c <- cowplot::plot_grid(plotlist = tmp, ncol = 1, rel_heights = c(1, 15))

# create node names as plot objects
node_names <- list()
for (i in seq_along(node_vec)) {
  node_names[[i]] <- ggplot() + theme_void() + annotate("text", x = 0, y = 0, size = subheading_size, angle = 90,
                                                        label = sprintf("nodes = %s", node_vec[i]))
}

# connect node name plots and add buffer
tmp <- cowplot::plot_grid(plotlist = node_names, ncol = 1)
plot_empty <- ggplot() + theme_void()
tmp <- cowplot::plot_grid(plotlist = list(plot_empty, tmp), ncol = 1, rel_heights = c(1, 15))
tmp <- cowplot::plot_grid(plotlist = list(plot_empty, tmp), ncol = 1, rel_heights = c(1, 15))

# add node names
plot_c <- cowplot::plot_grid(plotlist = list(tmp, plot_c), nrow = 1, rel_widths = c(1, 15))

# save plot to file
ggsave("FigureS3_bias/sim_bias.pdf", plot = plot_c, width = 14, height = 8)
ggsave("FigureS3_bias/sim_bias.png", plot = plot_c, width = 14, height = 8, dpi = 100)
