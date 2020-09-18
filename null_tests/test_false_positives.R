
# test_false_posities.R
#
# Author: Bob Verity
# Date: 2020-09-17
#
# Purpose:
# TODO
#
# ------------------------------------------------------------------

library(ggplot2)
library(cowplot)

# load PlasmoMAPI package (TODO - replace with instructions to install from github)
#devtools::load_all("/Users/rverity/Dropbox/Bob/Work/My Programs/Barriers to gene flow/PlasmoMAPI")

set.seed(2)

# define parameters
reps <- 1e2
L <- 5
n_deme <- 30
border_coords <- data.frame(long = c(-L, -L, L, L, -L),
                            lat = c(-L, L, L, -L, -L))

# create new PlasmoMAPI project
p <- pm_project()

# ------------------------------------------------------------------

t0 <- Sys.time()

# repeat analysis multiple times
which_lower <- which_upper <- list()
hex_values <- hex_values2 <- hex_values3 <- hex_coverage <- y_obs <- list()
sig <- rep(NA, reps)
for (i in seq_len(reps)) {
  message(sprintf("rep %s", i))
  
  # define node distribution
  node_df <- data.frame(long = runif(n_deme, -L, L),
                        lat = runif(n_deme, -L, L))
  
  # make symmetric matrix of statistical distance (random noise)
  stat_distance <- matrix(rnorm(n_deme^2), n_deme, n_deme)
  stat_distance <- as.matrix(as.dist(stat_distance))
  
  # load node coordinates
  p <- load_coords(p, node_df$long, node_df$lat, check_delete_output = FALSE)
  
  # set up map (first time only)
  if (is.null(p$map)) {
    p <- create_map(p, border_coords = border_coords, hex_width = 0.5)
  }
  
  # assign edges to hexes
  p <- assign_map(p, eccentricity = 0.9)
  
  # load data
  p <- load_data(p, stat_distance, check_delete_output = FALSE)
  
  # run analysis
  p <- pm_analysis(p, n_perms = 1e3, n_breaks = 1,
                   min_dist = 0, max_dist = Inf,
                   min_group_size = 5)
  
  
  #plot(p$output$hex_values)
  #abline(h = c(-1.96, 1.96))
  sig[i] <- mean(abs(p$output$hex_values) > 1.96, na.rm = TRUE)
  
  #p$output$hex_values <- p$output$hex_values2
  p$output$n_eff <- 100
  
  # get significant hexes
  s <- get_significant_hexes(p)
  which_lower[[i]] <- s$which_lower
  which_upper[[i]] <- s$which_upper
  hex_values[[i]] <- p$output$hex_values
  hex_values2[[i]] <- p$output$hex_values2
  hex_values3[[i]] <- p$output$hex_values3
  hex_coverage[[i]] <- p$output$hex_coverage
  y_obs[[i]] <- p$output$y_obs
}
Sys.time() - t0

# count up false positives overall
any_pos <- mapply(length, which_lower) + mapply(length, which_upper)
sum(any_pos > 0)
mean(any_pos > 0) * 100

hex_values <- do.call(rbind, hex_values)
hex_values2 <- do.call(rbind, hex_values2)
hex_values3 <- do.call(rbind, hex_values3)
hex_coverage <- do.call(rbind, hex_coverage)

z <- hex_values[hex_coverage > 20]

sd(z)

plot(density(as.vector(z)[!is.na(as.vector(z))]))
xv <- seq(-5, 5, l = 1001)
lines(xv, dnorm(xv), col = 2)
lines(xv, dnorm(xv, sd = sd(z, na.rm = T)), col = 3)
#lines(xv, dt(xv, df = 3))

#plot(rowMeans(z > 1.96 | z < -1.96, na.rm = TRUE))
#mean(z > 1.96 | z < -1.96, na.rm = TRUE) * 100

#plot(sig)
mean(sig)
