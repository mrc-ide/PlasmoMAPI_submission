
# sim_false_positives_spatial_regular.R
#
# Author: Bob Verity
# Date: 2020-09-17
#
# Purpose:
# Equivalent to sim_false_positives_spatial_random.R, except nodes are arranged
# in a regular grid. This grid stays the same between simulations, but the
# pairwise statistical values are drawn each time from the standard normal
# distribution.
#
# ------------------------------------------------------------------

# load PlasmoMAPI package (TODO - replace with instructions to install from github)
#devtools::load_all("/Users/rverity/Dropbox/Bob/Work/My Programs/Barriers to gene flow/PlasmoMAPI")

set.seed(1)

# define parameters
reps <- 1e4
L <- 5
border_coords <- data.frame(long = c(-L, -L, L, L, -L),
                            lat = c(-L, L, L, -L, -L))

# define node distribution
node_df <- expand.grid(long = -4:4, lat = -4:4)
n_deme <- nrow(node_df)

# create new PlasmoMAPI project
p <- pm_project()

# ------------------------------------------------------------------

t0 <- Sys.time()

# repeat analysis multiple times
which_lower <- which_upper <- list()
for (i in seq_len(reps)) {
  message(sprintf("rep %s", i))
  
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
  p <- assign_map(p, eccentricity = 0.95)
  
  # load data
  p <- load_data(p, stat_distance, check_delete_output = FALSE)
  
  # run analysis
  p <- pm_analysis(p, n_perms = 1e3, n_breaks = 100,
                   min_dist = 0, max_dist = Inf,
                   min_group_size = 5)
  
  
  # store results
  s <- get_significant_hexes(p, FDR = 0.05)
  which_lower[[i]] <- s$which_lower
  which_upper[[i]] <- s$which_upper
}
Sys.time() - t0

# make return list and save to file
ret <- list(which_lower = which_lower,
            which_upper = which_upper,
            p = p)
saveRDS(ret, "FigureS1_false_positives_spatial/sim_false_positives_spatial_regular.rds")

