
# sim_bias.R
#
# Author: Bob Verity
# Date: 2020-09-22
#
# Purpose:
# Explores how accurately the PlasmoMAPI method can identify known
# barriers/corridors. A series of setups are looped over, including random node
# positions vs. a regular grid, number of nodes and degree of eccentricity.
# Simulations are repeated multiple times where needed. Results are stored in a
# big list and saved as an RDS object for reading into a separate plotting
# script.
#
# ------------------------------------------------------------------

# load PlasmoMAPI package (TODO - replace with instructions to install from github)
#devtools::load_all("/Users/rverity/Dropbox/Bob/Work/My Programs/Barriers to gene flow/PlasmoMAPI")

set.seed(1)

# define parameters
reps <- 1e2
L <- 5
node_vec <- c(289, 81, 25)
ecc_vec <- c(0.9, 0.99)
border_coords <- data.frame(long = c(-L, -L, L, L, -L),
                            lat = c(-L, L, L, -L, -L))

# define spatial barriers/corridors
#barrier_list <- list(data.frame(long = c(-1,1,1,-1,-1)*0.5, lat = c(1,1,-1,-1,1)*5 + 2),
#                     data.frame(long = c(-1,1,1,-1,-1)*2, lat = c(1,1,-1,-1,1)*0.5 - 4))
barrier_list <- list(data.frame(long = c(-4,-3,-3,2,2,-4,-4), lat = c(4,4,-2,-2,-3,-3,4)),
                     data.frame(long = c(2,3,3,2,2), lat = c(-1,-1,-4,-4,-1)))
barrier_penalty <- c(5, -5)

# create new PlasmoMAPI project
p <- pm_project()

# ------------------------------------------------------------------

t0 <- Sys.time()

# loop through setups
list_p <- list_sim <- list()
for (i_arrange in 1:2) {
  list_p[[i_arrange]] <- list_sim[[i_arrange]] <- list()
  for (i_nodes in seq_along(node_vec)) {
    list_sim[[i_arrange]][[i_nodes]] <- list()
    for (i_ecc in seq_along(ecc_vec)) {
      list_sim[[i_arrange]][[i_nodes]][[i_ecc]] <- list()
      
      # loop through reps
      reps_actual <- ifelse(i_arrange == 1, reps, 1)
      for (i in seq_len(reps_actual)) {
        message(sprintf("sampling %s/2, nodes %s/%s, ecc %s/%s, rep %s/%s",
                        i_arrange, i_nodes, length(node_vec), i_ecc, length(ecc_vec), i, reps_actual))
        
        # define node distribution
        if (i_arrange == 1) {
          n_deme <- node_vec[i_nodes]
          node_df <- data.frame(long = runif(n_deme, -L, L),
                                lat = runif(n_deme, -L, L))
        } else {
          if (i_nodes == 1) {
            node_df <- expand.grid(long = seq(-4, 4, 0.5), lat = seq(-4, 4, 0.5))
            n_deme <- nrow(node_df)
          } else if (i_nodes == 2) {
            node_df <- expand.grid(long = seq(-4, 4, 1), lat = seq(-4, 4, 1))
            n_deme <- nrow(node_df)
          } else {
            node_df <- expand.grid(long = seq(-4, 4, 2), lat = seq(-4, 4, 2))
            n_deme <- nrow(node_df)
          }
        }
        
        # draw pairwise distances
        sim_dist <- get_barrier_intersect(node_df$long, node_df$lat,
                                          barrier_list = barrier_list,
                                          barrier_penalty = barrier_penalty,
                                          barrier_method = 2,
                                          max_barrier_range = Inf)
        #sim_dist <- sim_dist + rnorm(length(sim_dist), sd = 1)
        stat_distance <- as.matrix(as.dist(sim_dist))
        
        # load node coordinates
        p <- load_coords(p, node_df$long, node_df$lat, check_delete_output = FALSE)
        
        # set up map (first time only)
        if (is.null(p$map)) {
          p <- create_map(p, border_coords = border_coords, hex_width = 0.5)
        }
        
        # assign edges to hexes
        p <- assign_map(p, eccentricity = ecc_vec[i_ecc])
        
        # load data
        p <- load_data(p, stat_distance, check_delete_output = FALSE)
        #plot_dist(p)
        
        # run analysis
        p <- pm_analysis(p, n_perms = 1e3, n_breaks = 50,
                         min_dist = 0, max_dist = Inf,
                         min_group_size = 5)
        #plot(p)
        
        # store results
        list_sim[[i_arrange]][[i_nodes]][[i_ecc]][[i]] <- list(z = p$output$hex_values,
                                                               signif = get_significant_hexes(p, FDR = 0.05),
                                                               coverage = p$output$hex_coverage)
      }
      
    }
    list_p[[i_arrange]][[i_nodes]] <- p
  }
}
Sys.time() - t0

# save results to file
ret <- list(list_sim = list_sim,
            list_p = list_p,
            node_vec = node_vec,
            ecc_vec = ecc_vec,
            barrier_list = barrier_list)
saveRDS(ret, "FigureS3_bias/sim_bias.rds")
