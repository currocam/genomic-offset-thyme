library(tidyverse)
library(ggpubr)
theme_set(theme_pubr())
options(digits = 3)        # number of digits printed by R default (vectors, data.frames, lists)
options(pillar.sigfig = 3) # number of digits printed by tibbles default.

text_base_size   <- 16   # in pt
fig.witdh        <- 180  # in mm
fig.height       <- 125  # in mm

# Set all text in plots to same size
ggplot2::theme_set(
  cowplot::theme_cowplot(
    font_size = text_base_size, rel_small = 1,
    rel_tiny = 1, rel_large = 1
  )
)

# Setting text size inside plots (geom_text, geom_label etc.)
ggplot_text_size <- text_base_size / ggplot2::.pt
# Set locale if you want danish month names etc.
options(ggplot2.continuous.colour = ggplot2::scale_colour_viridis_c)
