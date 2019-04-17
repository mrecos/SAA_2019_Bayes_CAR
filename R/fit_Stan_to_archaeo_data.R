library("maptools");
library("spdep");
library("rgdal")
library("rstan");
library("bayesplot")

options(mc.cores = 6);

# source functions
source(file.path("R","nb_data_funs.R"))
source(file.path("R","archaeo_BYM_Functions.R"))

input_fishnet = sites_fishnet_stan

# convert fishnet to spatial
fishnet_sp <- as(input_fishnet, "Spatial")

#### MAKE SPATIAL COMPONENT
# make Neighbor list object
nb_fishnet = poly2nb(fishnet_sp);
# cast neighbors to graph (nodes and edges)
graph_fishnet = nb2graph(nb_fishnet);
N = graph_fishnet$N;
node1 = graph_fishnet$node1;
node2 = graph_fishnet$node2;
N_edges = graph_fishnet$N_edges;
scaling_factor = scale_nb_components(nb_fishnet)[1];

#### MAKE INDEPENDENT VAR COMPONENT
x = st_drop_geometry(input_fishnet) %>% 
  dplyr::select(mean_val1, mean_val2, mean_val3, mean_val4) %>% 
  as.matrix()
K = ncol(x)

#### MAKE DEPENDENT VAR COMPONENT
y = input_fishnet$count
# Set expected value (Should deal with this better)
E = rep(0,length(y));
# set pop > 0 so we can use log(pop) as offset
E[E < 1] = 0.01;

#### Compile Stan model
bym2_stan = stan_model(file.path("STAN","BYM2.stan"), verbose = FALSE);

#### Fit Stan model
bym2_fit = sampling(bym2_stan, data=list(N = N,
                                         N_edges = N_edges,
                                         node1 = node1,
                                         node2 = node2,
                                         y = y,
                                         E = E,
                                         K = K,
                                         x = x,
                                         scaling_factor = scaling_factor), 
                    control = list(adapt_delta = 0.97), 
                    chains=4, warmup=2000, iter=4000, save_warmup=FALSE, verbose = TRUE);


#### Inspect Stan results
print(bym2_fit, digits=3, pars=c("beta0", "rho", "logit_rho", "sigma", "betas[1]", "betas[2]",
"betas[3]", "betas[4]", "mu[1]", "mu[2]", "mu[3]", "phi[1]", "phi[2]", "phi[3]", "theta[1]", "theta[2]", "theta[3]"), probs=c(0.025, 0.5, 0.975))

#### Plot
# traceplot(bym2_fit, pars = c("beta0"), inc_warmup = FALSE, nrow = 2)

posterior <- as.matrix(bym2_fit)

# plot of parameter distributions
plot_title <- ggtitle("Posterior distributions",
                      "with medians and 80% intervals")
# betas
mcmc_areas(posterior,
           pars = c("betas[1]", "betas[2]", "betas[3]", "betas[4]"),
           prob = 0.8) + plot_title

# log odds of rho (logit)
color_scheme_set("blue")
mcmc_areas(posterior,
           pars = c("logit_rho"),
           prob = 0.8,
           prob_outer = 0.99, # 99%
           point_est = "median") + plot_title

# predictions of cell counts
color_scheme_set("red")
mcmc_areas(posterior,
           pars = c("mu[5]","mu[13]","mu[33]"),
           prob = 0.8) + plot_title

color_scheme_set("darkgray")
mcmc_scatter(
  as.matrix(bym2_fit),
  pars = c("rho", "sigma"), 
  np = nuts_params(bym2_fit), 
  np_style = scatter_style_np(div_color = "green", div_alpha = 0.8)
)

## plot predictions
mu <- rstan::extract(bym2_fit, "mu")
mu_plot <- data.frame(mu) %>% 
  gather(parameter, estimate) %>% 
  group_by(parameter) %>% 
  summarise(mean = quantile(estimate, probs=0.5),
            high = quantile(estimate, probs=0.9),
            low  = quantile(estimate, probs=0.1)) %>% 
  mutate(parameter = as.numeric(str_remove(parameter, "mu."))) %>% 
  dplyr::left_join(., data.frame("parameter" = seq(1:length(y)), obs = y), by="parameter")

ggplot(mu_plot) +
  geom_point(aes(x = obs, y = mean)) +
  geom_point(aes(x = obs, y = high), color = "blue") +
  geom_point(aes(x = obs, y = low), color = "red") +
  scale_x_continuous(limits = c(0,max(mu_plot$high))) +
  scale_y_continuous(limits = c(0,max(mu_plot$high))) +
  labs(x = "Observed", y = "Predicted") +
  geom_hline(yintercept = mean(y)) +
  geom_abline() +
  coord_fixed() +
  theme_bw() +
  theme(aspect.ratio=1)

# colmeans of phi
hist(rstan::extract(bym2_fit, "phi")[[1]])
# colnames of theta
hist(rstan::extract(bym2_fit, "theta")[[1]])


## map mcmc results
mu_plot_map <- mu_plot %>% 
  arrange(parameter) %>% 
  as.data.frame() %>% 
  bind_cols(.,input_fishnet) %>% 
  st_as_sf() %>% 
  arrange(desc(count))

# mapview(mu_plot_map, zcol = "mean") + mapview(sites)
bym2_mu <- tm_shape(mu_plot_map) +
  tm_fill("mean", midpoint=NA, breaks = seq(0,max(mu_plot_map$mean),0.25),
          palette = tmaptools::get_brewer_pal("Greens", n = 12)) +
  tm_borders(col = "gray50", lwd = 0.5, alpha = 0.5) +
  tm_shape(sites) +
  tm_borders(col = "black") + 
  tm_legend(show=FALSE)


# Phi plot
phi <- rstan::extract(bym2_fit, "phi")
phi_plot_map <- data.frame(phi) %>% 
  gather(parameter, estimate) %>% 
  group_by(parameter) %>% 
  summarise(mean = quantile(estimate, probs=0.5),
            high = quantile(estimate, probs=0.9),
            low  = quantile(estimate, probs=0.1)) %>% 
  mutate(parameter = as.numeric(str_remove(parameter, "phi."))) %>% 
  dplyr::left_join(., data.frame("parameter" = seq(1:length(y)), obs = y), by="parameter") %>% 
  arrange(parameter) %>% 
  as.data.frame() %>% 
  bind_cols(.,input_fishnet) %>% 
  st_as_sf() %>% 
  arrange(desc(count)) 

# mapview(phi_plot_map, zcol = "mean") + mapview(sites)
bym2_phi <- tm_shape(phi_plot_map) +
  tm_fill("mean", midpoint=0, breaks = seq(min(phi_plot_map$mean),max(phi_plot_map$mean),0.25),
          palette = tmaptools::get_brewer_pal("PiYG", n = 12)) +
  tm_borders(col = "gray50", lwd = 0.5, alpha = 0.5) +
  tm_shape(sites) +
  tm_borders(col = "black") + 
  tm_legend(show=FALSE)
tmap_save(bym2_phi, "bym2_phi.svg", width = 5, height = 5.5)


## Theta plot
theta <- rstan::extract(bym2_fit, "theta")
theta_plot_map <- data.frame(theta) %>% 
  gather(parameter, estimate) %>% 
  group_by(parameter) %>% 
  summarise(mean = quantile(estimate, probs=0.5),
            high = quantile(estimate, probs=0.9),
            low  = quantile(estimate, probs=0.1)) %>% 
  mutate(parameter = as.numeric(str_remove(parameter, "theta."))) %>% 
  dplyr::left_join(., data.frame("parameter" = seq(1:length(y)), obs = y), by="parameter") %>% 
  arrange(parameter) %>% 
  as.data.frame() %>% 
  bind_cols(.,input_fishnet) %>% 
  st_as_sf() %>% 
  arrange(desc(count))

# mapview(theta_plot_map, zcol = "mean") + mapview(sites)
bym2_theta <- tm_shape(theta_plot_map) +
  tm_fill("mean", midpoint=0,
          palette = tmaptools::get_brewer_pal("BrBG", n = 12)) +
  tm_shape(sites) +
  tm_borders(col = "black") + 
  tm_legend(show=FALSE)
tmap_save(bym2_theta, "bym2_theta.svg", width = 5, height = 5.5)

## convolved_re plot
convolved_re <- rstan::extract(bym2_fit, "convolved_re")
convolved_re_plot_map <- data.frame(convolved_re) %>% 
  gather(parameter, estimate) %>% 
  group_by(parameter) %>% 
  summarise(mean = quantile(estimate, probs=0.5),
            high = quantile(estimate, probs=0.9),
            low  = quantile(estimate, probs=0.1)) %>% 
  mutate(parameter = as.numeric(str_remove(parameter, "convolved_re."))) %>% 
  dplyr::left_join(., data.frame("parameter" = seq(1:length(y)), obs = y), by="parameter") %>% 
  arrange(parameter) %>% 
  as.data.frame() %>% 
  bind_cols(.,input_fishnet) %>% 
  st_as_sf() %>% 
  arrange(desc(count))

# mapview(theta_plot_map, zcol = "mean") + mapview(sites)
bym2_convolved_re <- tm_shape(convolved_re_plot_map) +
  tm_fill("mean", midpoint=0, breaks = seq(min(convolved_re_plot_map$mean),
                                           max(convolved_re_plot_map$mean),0.25),
          palette = tmaptools::get_brewer_pal("BrBG", n = 12)) +
  tm_borders(col = "gray50", lwd = 0.5, alpha = 0.5) +
  tm_shape(sites) +
  tm_borders(col = "black") + 
  tm_legend(show=FALSE)
tmap_save(bym2_convolved_re, "bym2_convolved_re.svg", width = 5, height = 5.5)

###
bym2_conved_mu  <- convolved_re_plot_map %>% 
  dplyr::left_join(., st_drop_geometry(mu_plot_map), by = "parameter") %>% 
  mutate(conved_mu_mean = mean.x * mean.y,
         conved_mu_low  = low.x  * low.y,
         conved_mu_high = high.x * high.y) %>% 
  arrange(parameter) %>% 
  as.data.frame() %>% 
  bind_cols(.,input_fishnet) %>% 
  st_as_sf()

tm_shape(bym2_conved_mean) +
  tm_fill("conved_mu_mean", midpoint=NA, 
          breaks = seq(0,max(bym2_conved_mu$conved_mu_mean),0.25),
          palette = tmaptools::get_brewer_pal("BuPu", n = max(bym2_conved_mu$conved_mu_mean)+1)) +
  tm_borders(col = "gray50", lwd = 0.5, alpha = 0.5) +
  tm_shape(sites_facet) +
  tm_borders(col = "black")

mu_conv_mu_compare <- bym2_conved_mu %>% 
  st_drop_geometry() %>% 
  select(conved_mu_mean, observed = obs.y, bym2_mu = mean.y) %>% 
  gather(parameter, value, -observed)

ggplot(mu_conv_mu_compare, aes(x = value, y = observed, group = parameter, color = parameter)) +
  geom_point() +
  geom_abline() +
  geom_smooth(method = "lm") +
  coord_equal() +
  theme_bw()



