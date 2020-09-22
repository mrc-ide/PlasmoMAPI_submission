
# plot_false_positives_ecc.R
#
# Author: Bob Verity
# Date: 2020-09-19
#
# Purpose:
# Read in results of simulations and produce plot summarising the proportion of
# false-positive results as a function of eccentricity.
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
sim_random <- readRDS("FigureS2_false_positives_ecc/sim_false_positives_ecc_random.rds")
sim_regular <- readRDS("FigureS2_false_positives_ecc/sim_false_positives_ecc_regular.rds")
ecc_vec <- sim_random$ecc_vec

# get results into dataframes
df_random_lower <- do.call(rbind, mapply(get_binom_exact, sim_random$which_lower, SIMPLIFY = FALSE))
df_random_upper <- do.call(rbind, mapply(get_binom_exact, sim_random$which_upper, SIMPLIFY = FALSE))
df_random <- rbind(cbind(df_random_lower, direction = "lower"),
                   cbind(df_random_upper, direction = "upper"))

df_regular_lower <- do.call(rbind, mapply(get_binom_exact, sim_regular$which_lower, SIMPLIFY = FALSE))
df_regular_upper <- do.call(rbind, mapply(get_binom_exact, sim_regular$which_upper, SIMPLIFY = FALSE))
df_regular <- rbind(cbind(df_regular_lower, direction = "lower"),
                    cbind(df_regular_upper, direction = "upper"))

# produce plots
plot1 <- ggplot(cbind(df_random, ecc = ecc_vec)) + theme_bw() +
  geom_pointrange(aes(x = ecc, y = proportion, ymin = lower, ymax = upper,
                      color = direction), position = position_dodge(width = 0.03)) +
  xlim(0, 1) + ylim(0, 5) + xlab("eccentricity") + ylab("familywise error-rate (%)") +
  ggtitle("A) random")
plot1

plot2 <- ggplot(cbind(df_regular, ecc = ecc_vec)) + theme_bw() +
  geom_pointrange(aes(x = ecc, y = proportion, ymin = lower, ymax = upper,
                      color = direction), position = position_dodge(width = 0.03)) +
  xlim(0, 1) + ylim(0, 5) + xlab("eccentricity") + ylab("familywise error-rate (%)") +
  ggtitle("B) regular")
plot2

# produce combined plot
plot_c <- cowplot::plot_grid(plot1, plot2)
plot_c

# save plot to file
ggsave("FigureS2_false_positives_ecc/false_positives_ecc.pdf", width = 10, height = 4)
ggsave("FigureS2_false_positives_ecc/false_positives_ecc.png", width = 10, height = 4, dpi = 100)
