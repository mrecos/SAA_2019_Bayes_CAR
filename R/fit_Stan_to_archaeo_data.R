library("maptools");
library("spdep");
library("rgdal")
library("rstan");
library("bayesplot")
# library("tidyr")
# library("stringr")

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
# Set expected value (NOT SURE HOW TO DEAL WITH THIS!)
E = rep(0,length(y));
# set pop > 0 so we can use log(pop) as offset
E[E < 1] = 0.01;

#### MAKE TEST DATA COMPONENT
# N_test <- 30
# x_test <- x[1:N_test,]
# E_test <- rep(0.01,N_test)

# N_test <- nrow(test_fishnet)
# x_test <-  st_drop_geometry(test_fishnet) %>% 
#   dplyr::select(mean_val1, mean_val2, mean_val3, mean_val4) %>% 
#   as.matrix()
# E_test <- rep(0.01,N_test)

#### Compile Stan model
bym2_stan = stan_model(file.path("STAN","BYM2.stan"), verbose = FALSE);


# prior_dist <- stan(file=file.path("STAN","BYM2.stan"),
#                    iter=1000, warmup=0, chains=1,
#                    seed=4838282, algorithm="Fixed_param")

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

mapview(mu_plot_map, zcol = "mean") + mapview(sites)

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

mapview(phi_plot_map, zcol = "mean") + mapview(sites)

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

mapview(theta_plot_map, zcol = "mean") + mapview(sites)




