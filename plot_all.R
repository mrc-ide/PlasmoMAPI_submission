
# plot_all.R
#
# Author: Bob Verity
# Date: 2020-09-25
#
# Purpose:
# Runs all scripts that produce final plots from simulated output.
#
# ------------------------------------------------------------------

# main paper scripts
source("Figure1_schematic/make_figure1.R")

source("Figure3_DRC/plot_DRC.R")


# supplementary material scripts
source("FigureS1_false_positives_spatial/plot_false_positives.R")

source("FigureS2_false_positives_ecc/plot_false_positives_ecc.R")

source("FigureS3_bias/plot_bias.R")
