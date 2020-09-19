
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

# load PlasmoMAPI package (TODO - replace with instructions to install from github)
#devtools::load_all("/Users/rverity/Dropbox/Bob/Work/My Programs/Barriers to gene flow/PlasmoMAPI")

# read in simulation results
sims <- readRDS("null_tests/sim_false_positives_random.rds")
which_lower <- sims$which_lower
which_upper <- sims$which_upper
p <- sims$p
reps <- length(which_lower)
n_hex <- length(p$map$hex)

# count up total false positives
any_pos <- mapply(length, which_lower) + mapply(length, which_upper)
mean(any_pos > 0) * 100

# tally up false positives per hex
unlist_lower <- unlist(which_lower)
hex_lower <- tabulate(unlist_lower, nbins = n_hex)

unlist_upper <- unlist(which_upper)
hex_upper <- tabulate(unlist_upper, nbins = n_hex)


# plot hex map of false positives
col_lower <- grey(seq(0, 1, l = max(hex_lower)))
plot(p$map$hex, col = col_lower[hex_lower])

col_upper <- grey(seq(0, 1, l = max(hex_upper)))
plot(p$map$hex, col = col_upper[hex_upper])
