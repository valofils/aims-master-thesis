library(sf)
library(ggplot2)
library(rnaturalearth)
library(rnaturalearthdata)
library(rnaturalearthhires)
library(cowplot)
library(dplyr)

# --- Brazil states (main panel) ---
brazil <- ne_states(country = "brazil", returnclass = "sf")

p_brazil <- ggplot(brazil) +
  geom_sf(fill = "#ADD8E6", color = "white", linewidth = 0.3) +
  theme_void()

# --- South America inset ---
world <- ne_countries(scale = "medium", returnclass = "sf")

p_world <- ggplot(world) +
  geom_sf(fill = "grey85", color = "grey60", linewidth = 0.2) +
  geom_sf(data = filter(world, name == "Brazil"),
          fill = "#00008B", color = "white", linewidth = 0.3) +
  coord_sf(xlim = c(-85, -30), ylim = c(-60, 15)) +
  theme_void() +
  theme(panel.border = element_rect(color = "grey40", fill = NA))

# --- Combine ---
ggdraw(p_brazil) +
  draw_plot(p_world, x = 0.6, y = 0.55, width = 0.38, height = 0.38) +
  theme(plot.caption = element_text(hjust = 0.5, family = "serif"))

ggsave("brazil_locator_map.pdf", width = 8, height = 6)