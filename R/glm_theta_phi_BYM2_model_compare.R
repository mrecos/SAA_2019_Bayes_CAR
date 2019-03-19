library("maptools");
library("spdep");
library("rgdal")
library("rstan");
library("bayesplot")
# library("tidyr")
# library("stringr")

rmse <- function(obs,pred){
  sqrt(mean((obs-pred)^2,na.rm=TRUE))
}

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

#### Compile Stan model
glm_stan   = stan_model(file.path("STAN","pois_glm_MDH.stan"), verbose = FALSE);
theta_stan = stan_model(file.path("STAN","pois_glm_theta_MDH.stan"), verbose = FALSE);
phi_stan   = stan_model(file.path("STAN","pois_glm_phi_MDH.stan"), verbose = FALSE);
bym2_stan  = stan_model(file.path("STAN","BYM2.stan"), verbose = FALSE);
 
#### Fit Stan model
glm_stan_fit = sampling(glm_stan, data=list(N = N,
                                         y = y,
                                         E = E,
                                         K = K,
                                         x = x), 
                    control = list(adapt_delta = 0.97), 
                    chains=4, warmup=2000, iter=4000, save_warmup=FALSE, verbose = TRUE);

theta_stan_fit = sampling(theta_stan, data=list(N = N,
                                            y = y,
                                            E = E,
                                            K = K,
                                            x = x), 
                        control = list(adapt_delta = 0.97), 
                        chains=4, warmup=2000, iter=4000, save_warmup=FALSE, verbose = TRUE);

phi_stan_fit = sampling(theta_stan, data=list(N = N,
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

print(glm_stan_fit, digits=3, pars=c("beta0", "sigma", "betas[1]", "betas[2]", "betas[3]", "betas[4]",
                                 "mu[1]", "mu[10]", "mu[20]", "mu[30]", "mu[40]", "mu[50]"),
      probs=c(0.025, 0.5, 0.975))
print(theta_stan_fit, digits=3, pars=c("beta0", "sigma", "betas[1]", "betas[2]", "betas[3]", "betas[4]",
                                     "mu[1]", "mu[10]", "mu[20]", "mu[30]", "mu[40]", "mu[50]"),
      probs=c(0.025, 0.5, 0.975))
print(phi_stan_fit, digits=3, pars=c("beta0", "sigma", "betas[1]", "betas[2]", "betas[3]", "betas[4]",
                                           "mu[1]", "mu[10]", "mu[20]", "mu[30]", "mu[40]", "mu[50]"),
      probs=c(0.025, 0.5, 0.975))
print(bym2_fit, digits=3, pars=c("beta0", "sigma", "betas[1]", "betas[2]", "betas[3]", "betas[4]",
                                           "mu[1]", "mu[10]", "mu[20]", "mu[30]", "mu[40]", "mu[50]"),
      probs=c(0.025, 0.5, 0.975))

glm_mu   <- rstan::extract(glm_stan_fit, "mu") %>% 
  data.frame(.) %>% 
  gather(parameter, estimate) %>% 
  group_by(parameter) %>% 
  summarise(mean = quantile(estimate, probs=0.5),
            high = quantile(estimate, probs=0.9),
            low  = quantile(estimate, probs=0.1)) %>% 
  mutate(parameter = as.numeric(str_remove(parameter, "mu."))) %>% 
  dplyr::left_join(., data.frame("parameter" = seq(1:length(y)), obs = y), by="parameter")
theta_mu <- rstan::extract(theta_stan_fit, "mu")%>% 
  data.frame(.) %>% 
  gather(parameter, estimate) %>% 
  group_by(parameter) %>% 
  summarise(mean = quantile(estimate, probs=0.5),
            high = quantile(estimate, probs=0.9),
            low  = quantile(estimate, probs=0.1)) %>% 
  mutate(parameter = as.numeric(str_remove(parameter, "mu."))) %>% 
  dplyr::left_join(., data.frame("parameter" = seq(1:length(y)), obs = y), by="parameter")
phi_mu   <- rstan::extract(phi_stan_fit, "mu")%>% 
  data.frame(.) %>% 
  gather(parameter, estimate) %>% 
  group_by(parameter) %>% 
  summarise(mean = quantile(estimate, probs=0.5),
            high = quantile(estimate, probs=0.9),
            low  = quantile(estimate, probs=0.1)) %>% 
  mutate(parameter = as.numeric(str_remove(parameter, "mu."))) %>% 
  dplyr::left_join(., data.frame("parameter" = seq(1:length(y)), obs = y), by="parameter")
bym2_mu  <- rstan::extract(bym2_fit, "mu")%>% 
  data.frame(.) %>% 
  gather(parameter, estimate) %>% 
  group_by(parameter) %>% 
  summarise(mean = quantile(estimate, probs=0.5),
            high = quantile(estimate, probs=0.9),
            low  = quantile(estimate, probs=0.1)) %>% 
  mutate(parameter = as.numeric(str_remove(parameter, "mu."))) %>% 
  dplyr::left_join(., data.frame("parameter" = seq(1:length(y)), obs = y), by="parameter")

mu_compare <- data.frame(id = glm_mu$parameter,
                         glm_mean   = round(glm_mu$mean,3),
                         theta_mean = round(theta_mu$mean,3),
                         phi_mean   = round(phi_mu$mean,3),
                         bym2_mean  = round(bym2_mu$mean,3),
                         observed   = round(bym2_mu$obs,3)) %>% 
  rowwise() %>% 
  mutate(glm_RMSE   = rmse(observed, glm_mean),
         theta_RMSE = rmse(observed, theta_mean),
         phi_RMSE   = rmse(observed, phi_mean),
         bym2_RMSE  = rmse(observed, bym2_mean)) %>% 
  arrange(desc(observed))
  
round(colMeans(mu_compare),3)

ggplot(gather(mu_compare, model, rmse, -id, -observed) %>% filter(str_detect(model,"RMSE")),
       aes(x = model, y = rmse, color = model)) +
  geom_jitter( width = 0.25) +
  scale_y_log10() +
  theme_bw() +
  theme(
    legend.position = "none"
  )

ggplot(gather(mu_compare, model, mean, -id, -observed) %>% filter(str_detect(model,"mean")),
       aes(x = mean, y = observed, group = model, color = model)) +
  geom_point() +
  facet_wrap(~model) +
  scale_x_continuous(limits = c(0,max(mu_compare$observed))) +
  labs(x = "Prediction Mean") +
  coord_equal() +
  geom_abline(linetype = "dashed") +
  theme_bw() +
  theme(
    legend.position = "none"
  )

###
mu_plot_map <- mu_compare %>% 
  arrange(id) %>% 
  as.data.frame() %>% 
  bind_cols(.,input_fishnet) %>% 
  st_as_sf() %>% 
  arrange(desc(count))

# mapview(mu_plot_map, zcol = "glm_mean") + mapview(sites)
library(tmap)

data(World, NLD_muni, NLD_prov, land, metro)

current.mode <- tmap_mode("plot")

# CASE 1: Facets defined by constant values
xx <- gather(mu_plot_map, model, mean, -id, -observed, -geometry) %>% filter(str_detect(model,"_mean"))
sites_facet <- rbind(sites,sites,sites,sites) %>% 
  mutate(model = rep(c("glm_mean","theta_mean","phi_mean","bym2_mean"),each=nrow(sites)))
tm_shape(xx) +
  tm_fill("mean", midpoint=NA, breaks = seq(0,max(xx$mean),0.5),
          palette = tmaptools::get_brewer_pal("RdPu", n = max(y)+1)) +
  tm_borders(col = "gray50", lwd = 0.5, alpha = 0.5) +
  # tm_shape(sites_facet) +
  # tm_borders(col = "black") +
  tm_facets(by = "model", nrow = 2)


#### convoluting times the mean
mean(-log(dpois(mu_compare$observed, lambda=mu_compare$glm_mean)))
mean(-log(dpois(mu_compare$observed, lambda=mu_compare$phi_mean)))
mean(-log(dpois(mu_compare$observed, lambda=mu_compare$theta_mean)))
mean(-log(dpois(mu_compare$observed, lambda=mu_compare$bym2_mean)))
mean(-log(dpois(bym2_conved_mu$obs.y, lambda=bym2_conved_mu$mean.y)))
mean(-log(dpois(bym2_conved_mu$obs.y, lambda=ifelse(bym2_conved_mu$conved_mu_mean<=0,0.01,bym2_conved_mu$conved_mu_mean))))

bym2_conved_mu  <- rstan::extract(bym2_fit, "convolved_re")%>% 
  data.frame(.) %>% 
  gather(parameter, estimate) %>% 
  group_by(parameter) %>% 
  summarise(mean = quantile(estimate, probs=0.5),
            high = quantile(estimate, probs=0.9),
            low  = quantile(estimate, probs=0.1)) %>% 
  mutate(parameter = as.numeric(str_remove(parameter, "convolved_re."))) %>% 
  dplyr::left_join(., data.frame("parameter" = seq(1:length(y)), obs = y), by="parameter") %>% 
  dplyr::left_join(., bym2_mu, by = "parameter") %>% 
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
