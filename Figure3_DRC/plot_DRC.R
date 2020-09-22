
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

library(rworldmap)
library(cowplot)

# load PlasmoMAPI package (TODO - replace with instructions to install from github)
#devtools::load_all("/Users/rverity/Dropbox/Bob/Work/My Programs/Barriers to gene flow/PlasmoMAPI")

# plotting parameters
point_stroke <- 0.2

# read in data from file
x <- readRDS("Figure3_DRC/run_DRC.rds")

# get DRC map
world_map <- rworldmap::getMap(resolution = "coarse")
drc_map <- subset(world_map, ISO3 == "COD")
drc_coords <- as.data.frame(drc_map@polygons[[1]]@Polygons[[1]]@coords)
names(drc_coords) <- c("long", "lat")

# create base map
plot_base <- ggplot() + theme_bw(base_size = 14) +
  theme(panel.background = element_rect(fill = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank()) +
  geom_polygon(aes(long, lat, group = group),
               size = 0.6, color = grey(0.4),
               fill = grey(0.3), data = world_map) +
  theme(plot.margin = unit(c(0.1, 0, 0, 0), "cm"))

# produce plots
plot1 <- plot_map(x$p1, min_hex_coverage = 10, base_plot = plot_base, point_stroke = point_stroke) +
  coord_sf(xlim = c(12, 32), ylim = c(-13, 6)) +
  ggplot2::ggtitle("A) all ranges") +
  theme(legend.position = "bottom")
plot1

plot2 <- plot_map(x$p2, min_hex_coverage = 10, base_plot = plot_base, point_stroke = point_stroke) +
  coord_sf(xlim = c(12, 32), ylim = c(-13, 6)) +
  ggplot2::ggtitle("B) short range (0-500km)")
plot2

plot3 <- plot_map(x$p3, min_hex_coverage = 10, base_plot = plot_base, point_stroke = point_stroke) +
  coord_sf(xlim = c(12, 32), ylim = c(-13, 6)) +
  ggplot2::ggtitle("C) long range (1500-2000km)")
plot3


# produce combined plot
plot_c <- cowplot::plot_grid(plot1, cowplot::plot_grid(plot2, plot3, nrow = 2), nrow = 1, rel_heights = c(2,1))
plot_c

# save to file
ggsave("Figure3_DRC/DRC_analysis.pdf", plot = plot_c, width = 13, height = 7)
ggsave("Figure3_DRC/DRC_analysis.png", plot = plot_c, width = 13, height = 7, dpi = 100)
