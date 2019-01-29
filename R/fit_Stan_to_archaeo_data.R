library(maptools);
library(spdep);
library(rgdal)
library(rstan);
library("bayesplot")
options(mc.cores = 6);

# source functions
source(file.path("R","nb_data_funs.R"))
source(file.path("R","archaeo_BYM_Functions.R"))

# convert fishnet to spatial
fishnet.sp <- as(pass_fishnet2, "Spatial")

#### MAKE SPATIAL COMPONENT
# make Neighbor list object
nb_fishnet = poly2nb(fishnet.sp);
# cast neighbors to graph (nodes and edges)
graph_fishnet = nb2graph(nb_fishnet);
N = graph_fishnet$N;
node1 = graph_fishnet$node1;
node2 = graph_fishnet$node2;
N_edges = graph_fishnet$N_edges;
scaling_factor = scale_nb_components(nb_fishnet)[1];

#### MAKE INDEPENDENT VAR COMPONENT
x = st_drop_geometry(pass_fishnet2) %>% 
  dplyr::select(med_elev) %>% 
  as.matrix()
K = ncol(x)

#### MAKE DEPENDENT VAR COMPONENT
y = pass_fishnet2$count
# Set expected value (NOT SURE HOW TO DEAL WITH THIS!)
E = rep(0,length(y));
# set pop > 0 so we can use log(pop) as offset
E[E < 1] = 1;

#### Compile Stan model
bym2_stan = stan_model(file.path("STAN","BYM2.stan"), verbose = FALSE);

#### Fit Stan model
bym2_fit = sampling(bym2_stan, data=list(N,
                                         N_edges,
                                         node1,
                                         node2,
                                         y,
                                         x,
                                         K,
                                         E,
                                         scaling_factor), 
                    control = list(adapt_delta = 0.97), 
                    chains=3, warmup=2000, iter=3000, save_warmup=FALSE, verbose = TRUE);

#### Inspect Stan results
print(bym2_fit, digits=3, pars=c("beta0", "rho", "sigma", "mu[1]", "mu[2]", "mu[3]", "mu[50]", "phi[1]", "phi[2]", "phi[3]", "phi[50]", "theta[1]", "theta[2]", "theta[3]", "theta[50]"), probs=c(0.025, 0.5, 0.975))

mu <- rstan::extract(bym2_fit, "mu")
colMeans(mu[[1]])

#### Plot
posterior <- as.matrix(bym2_fit)

# plot of parameter distributions
plot_title <- ggtitle("Posterior distributions",
                      "with medians and 80% intervals")
mcmc_areas(posterior,
           pars = c("rho", "sigma"),
           prob = 0.8) + plot_title

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
  labs(x = "Predicted", y = "Observed") +
  geom_abline() +
  coord_fixed() +
  theme_bw() +
  theme(aspect.ratio=1)





