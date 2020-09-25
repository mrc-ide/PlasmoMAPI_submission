
# run_DRC
#
# Author: Bob Verity
# Date: 2020-09-19
#
# Purpose:
# Read in DRC pairwise IBD measures, produce maps at different distance scales
# and save results to file.
#
# ------------------------------------------------------------------

library(magrittr)
library(rworldmap)

# load PlasmoMAPI package (TODO - replace with instructions to install from github)
#devtools::load_all("/Users/rverity/Dropbox/Bob/Work/My Programs/Barriers to gene flow/PlasmoMAPI")

set.seed(1)

# get DRC border shape
world_map <- rworldmap::getMap(resolution = "coarse")
drc_map <- subset(world_map, ISO3 == "COD")
drc_coords <- as.data.frame(drc_map@polygons[[1]]@Polygons[[1]]@coords)
names(drc_coords) <- c("long", "lat")

# load data
raw_data <- pm_file("DRC_relatedness.rds")

# create new project and load coords
p <- pm_project() %>%
  load_coords(raw_data$coords$long, raw_data$coords$lat)

# create and assign map
p <- create_map(p, border_coords = drc_coords, hex_width = 0.3) %>%
  assign_map(eccentricity = 0.99)
#plot(p)

# load data
p <- load_data(p, raw_data$relatedness, check_delete_output = FALSE)
#plot_dist(p)

# run analysis
p <- pm_analysis(p, n_perms = 1e3, n_breaks = 100,
                 min_dist = 0, max_dist = Inf,
                 min_group_size = 5)

# plot coverage
#plot_coverage(p)

# run short-range analysis
p2 <- pm_analysis(p, n_perms = 1e3, n_breaks = 100,
                  min_dist = 0, max_dist = 500,
                  min_group_size = 5)

# plot coverage
#plot_coverage(p)

# run long-rage analysis
p3 <- pm_analysis(p, n_perms = 1e3, n_breaks = 100,
                  min_dist = 1500, max_dist = 2000,
                  min_group_size = 5)

# save results to file
ret <- list(p1 = p, p2 = p2, p3 = p3)
saveRDS(ret, "Figure3_DRC/run_DRC.rds")
