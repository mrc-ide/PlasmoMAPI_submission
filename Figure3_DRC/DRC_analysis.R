
# DRC_analysis.R
#
# Author: Bob Verity
# Date: 2020-09-19
#
# Purpose:
# Read in DRC pairwise IBD measures, produce maps at different distance scales,
# produce combined ggplot and print to file.
#
# ------------------------------------------------------------------

library(rworldmap)
library(magrittr)
library(cowplot)

# load PlasmoMAPI package (TODO - replace with instructions to install from github)
#devtools::load_all("/Users/rverity/Dropbox/Bob/Work/My Programs/Barriers to gene flow/PlasmoMAPI")

set.seed(1)

# plotting parameters
point_stroke <- 0.2

# get DRC map
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

# create base map
plot_base <- ggplot() + theme_bw(base_size = 12) +
  theme(panel.background = element_rect(fill = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank()) +
  geom_polygon(aes(long, lat, group = group),
               size = 0.6, color = grey(0.4),
               fill = grey(0.3), data = world_map) +
  theme(legend.position = "bottom")

# plot map
plot1 <- plot_map(p, min_hex_coverage = 10, base_plot = plot_base, point_stroke = point_stroke) +
  coord_sf(xlim = c(12, 32), ylim = c(-13, 6)) +
  ggplot2::ggtitle("A) All ranges")
plot1

# ------------------------------------------------

# run analysis
p2 <- pm_analysis(p, n_perms = 1e3, n_breaks = 100,
                  min_dist = 0, max_dist = 500,
                  min_group_size = 5)

# plot coverage
#plot_coverage(p)

# plot map
plot2 <- plot_map(p2, min_hex_coverage = 10, base_plot = plot_base, point_stroke = point_stroke) +
  coord_sf(xlim = c(12, 32), ylim = c(-13, 6)) +
  ggplot2::ggtitle("B) Short range (0-500km)")
plot2

# ------------------------------------------------

# run analysis
p3 <- pm_analysis(p, n_perms = 1e3, n_breaks = 100,
                  min_dist = 1500, max_dist = 2000,
                  min_group_size = 5)

# plot coverage
#plot_coverage(p)

# plot map
plot3 <- plot_map(p3, min_hex_coverage = 10, base_plot = plot_base, point_stroke = point_stroke) +
  coord_sf(xlim = c(12, 32), ylim = c(-13, 6)) +
  ggplot2::ggtitle("C) Long range (1500-2000km)")
plot3

# ------------------------------------------------

# produce combined plot
plot_c <- cowplot::plot_grid(plot1, plot2, plot3, nrow = 1)
plot_c

# save to file
ggsave("Figure3_DRC/DRC_analysis.pdf", width = 15, height = 6)
ggsave("Figure3_DRC/DRC_analysis.png", width = 15, height = 6, dpi = 100)
