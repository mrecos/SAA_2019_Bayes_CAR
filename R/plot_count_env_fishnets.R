ggplot() + 
  geom_sf(data = sites_fishnet3, color = "black", fill="white") + 
  theme_void() + 
  theme(
    panel.grid.major = element_line(colour = "white"),
    legend.position = "none"
  )

ggsave("./plots/clean.png", width = 6, height = 6)


ggplot() + 
  geom_sf(data = sites_fishnet3, color = "black", aes(fill = log(count+0.5))) + 
  theme_void() + 
  scale_fill_distiller(palette = "GnBu", direction = 1) +
  # scale_fill_viridis_d(option="A", direction = -1) +
  theme(
    panel.grid.major = element_line(colour = "white"),
    legend.position = "none"
    )

ggsave("./plots/count.png", width = 6, height = 6)


ggplot() + 
  geom_sf(data = sites_fishnet3, color = "black", aes(fill = mean_ed_h2)) + 
  theme_void() + 
  # scale_fill_viridis_c(option="A", direction = -1)+
  scale_fill_distiller(palette = "Spectral", direction = 1) +
  theme(
    panel.grid.major = element_line(colour = "white"),
    legend.position = "none"
  )
ggsave("./plots/env1.jpg", width = 6, height = 6)

ggplot() + 
  geom_sf(data = sites_fishnet3, color = "black", aes(fill = mean_cd_h7)) + 
  theme_void() + 
  # scale_fill_viridis_c(option="B", direction = -1)+
  scale_fill_distiller(palette = "PRGn", direction = -1) +
  theme(
    panel.grid.major = element_line(colour = "white"),
    legend.position = "none"
  )
ggsave("./plots/env2.jpg", width = 6, height = 6)

ggplot() + 
  geom_sf(data = sites_fishnet3, color = "black", aes(fill = mean_std_32c)) + 
  theme_void() + 
  # scale_fill_viridis_c(option="C", direction = -1)+
  scale_fill_distiller(palette = "PiYG", direction = 1) +
  theme(
    panel.grid.major = element_line(colour = "white"),
    legend.position = "none"
  )
ggsave("./plots/env3.jpg", width = 6, height = 6)


ggplot() + 
  geom_sf(data = sites_fishnet3, color = "black", aes(fill = mean_elev_2_drainh)) + 
  theme_void() + 
  # scale_fill_viridis_c(option="D", direction = -1)+
  scale_fill_distiller(palette = "PiYG", direction = 1) +
  theme(
    panel.grid.major = element_line(colour = "white"),
    legend.position = "none"
  )
ggsave("env4.svg", width = 6, height = 6)
