
# run_all.R
#
# Author: Bob Verity
# Date: 2020-09-25
#
# Purpose:
# Runs all scripts that are involved in generating data for plots. Some of these
# can take a very long time.
#
# ------------------------------------------------------------------

# main paper scripts
source("Figure3_DRC/run_DRC.R")


# supplementary material scripts
source("FigureS1_false_positives_spatial/sim_false_positives_random.R")
source("FigureS1_false_positives_spatial/sim_false_positives_regular.R")

source("FigureS2_false_positives_ecc/sim_false_positives_ecc_random.R")
source("FigureS2_false_positives_ecc/sim_false_positives_ecc_regular.R")

source("FigureS3_bias/sim_bias.R")
