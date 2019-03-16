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

####
siteI_envI <- sites_fishnet_stan

# count 1, vars 1
siteC_envC <- siteI_envI %>%
  mutate(count = 1,
         mean_val1 = 1,
         mean_val2 = 1,
         mean_val3 = 1,
         mean_val4 = 1)

# count 1, vars same
siteC_envI <- siteI_envI %>%
  mutate(count = 1)

# count 1, vars runif (0,1)
siteC_envR <- siteI_envI %>%
  mutate(count = 1,
         mean_val1 = runif(nrow(siteI_envI),0,1),
         mean_val2 = runif(nrow(siteI_envI),0,1),
         mean_val3 = runif(nrow(siteI_envI),0,1),
         mean_val4 = runif(nrow(siteI_envI),0,1))

# count 1, vars uniform
siteC_envU <- siteI_envI %>%
  mutate(count = 1, 
         mean_val1 = runif(nrow(siteI_envI),min(mean_val1),max(mean_val1)),
         mean_val2 = runif(nrow(siteI_envI),min(mean_val2),max(mean_val2)),
         mean_val3 = runif(nrow(siteI_envI),min(mean_val3),max(mean_val3)),
         mean_val4 = runif(nrow(siteI_envI),min(mean_val4),max(mean_val4)))

# count same, vars 1
siteI_envC <- siteI_envI %>%
  mutate(mean_val1 = 1,
         mean_val2 = 1,
         mean_val3 = 1,
         mean_val4 = 1)

# count same, vars runif (0,1)
siteI_envR <- siteI_envI %>%
  mutate(mean_val1 = runif(nrow(siteI_envI),0,1),
         mean_val2 = runif(nrow(siteI_envI),0,1),
         mean_val3 = runif(nrow(siteI_envI),0,1),
         mean_val4 = runif(nrow(siteI_envI),0,1))

# count same, vars runif in range
siteI_envU <- siteI_envI %>%
  mutate(mean_val1 = runif(nrow(siteI_envI),min(mean_val1),max(mean_val1)),
         mean_val2 = runif(nrow(siteI_envI),min(mean_val2),max(mean_val2)),
         mean_val3 = runif(nrow(siteI_envI),min(mean_val3),max(mean_val3)),
         mean_val4 = runif(nrow(siteI_envI),min(mean_val4),max(mean_val4)))

# count random, vars 1
siteR_envC <- siteI_envI %>%
  mutate(count = sample(c(0,1),nrow(siteI_envI),replace = TRUE),
         mean_val1 = 1,
         mean_val2 = 1,
         mean_val3 = 1,
         mean_val4 = 1)

# count random, vars same
siteR_envI <- siteI_envI %>%
  mutate(count = sample(c(0,1),nrow(siteI_envI),replace = TRUE))

# count sample (0,1), vars runif (0,1)
siteR_envR <- siteI_envI %>%
  mutate(count = sample(c(0,1),nrow(siteI_envI),replace = TRUE),
         mean_val1 = runif(nrow(siteI_envI),0,1),
         mean_val2 = runif(nrow(siteI_envI),0,1),
         mean_val3 = runif(nrow(siteI_envI),0,1),
         mean_val4 = runif(nrow(siteI_envI),0,1))

# count random, vars runif
siteR_envU <- siteI_envI %>%
  mutate(count = sample(c(0,1),nrow(siteI_envI),replace = TRUE),
         mean_val1 = runif(nrow(siteI_envI),min(mean_val1),max(mean_val1)),
         mean_val2 = runif(nrow(siteI_envI),min(mean_val2),max(mean_val2)),
         mean_val3 = runif(nrow(siteI_envI),min(mean_val3),max(mean_val3)),
         mean_val4 = runif(nrow(siteI_envI),min(mean_val4),max(mean_val4)))

# count uniform, vars 1
siteU_envC <- siteI_envI %>%
  mutate(count = sample(0:max(count), nrow(siteI_envI),replace = TRUE),
         mean_val1 = 1,
         mean_val2 = 1,
         mean_val3 = 1,
         mean_val4 = 1)

# count range, vars same
siteU_envI <- siteI_envI %>%
  mutate(count = sample(0:max(count), nrow(siteI_envI),replace = TRUE)) 

# count range, vars runif (0,1)
siteU_envR <- siteI_envI %>%
  mutate(count = sample(0:max(count), nrow(siteI_envI),replace = TRUE), 
         mean_val1 = runif(nrow(siteI_envI),0,1),
         mean_val2 = runif(nrow(siteI_envI),0,1),
         mean_val3 = runif(nrow(siteI_envI),0,1),
         mean_val4 = runif(nrow(siteI_envI),0,1))

# count range, vars runif range
siteU_envU <- siteI_envI %>%
  mutate(count = sample(0:max(count), nrow(siteI_envI),replace = TRUE),
         mean_val1 = runif(nrow(siteI_envI),min(mean_val1),max(mean_val1)),
         mean_val2 = runif(nrow(siteI_envI),min(mean_val2),max(mean_val2)),
         mean_val3 = runif(nrow(siteI_envI),min(mean_val3),max(mean_val3)),
         mean_val4 = runif(nrow(siteI_envI),min(mean_val4),max(mean_val4)))

####


# 
# input_fishnet_random = siteI_envI_random5
# 
# # convert fishnet to spatial
# fishnet_sp_random <- as(input_fishnet_random, "Spatial")
# 
# #### MAKE SPATIAL COMPONENT
# # make Neighbor list object
# nb_fishnet_random = poly2nb(fishnet_sp_random);
# # cast neighbors to graph (nodes and edges)
# graph_fishnet_random = nb2graph(nb_fishnet_random);
# N = graph_fishnet_random$N;
# node1 = graph_fishnet_random$node1;
# node2 = graph_fishnet_random$node2;
# N_edges = graph_fishnet_random$N_edges;
# scaling_factor = scale_nb_components(nb_fishnet_random)[1];
# 
# #### MAKE INDEPENDENT VAR COMPONENT
# x_random = st_drop_geometry(input_fishnet_random) %>% 
#   dplyr::select(mean_val1, mean_val2, mean_val3, mean_val4) %>% 
#   as.matrix()
# K = ncol(x_random)
# 
# #### MAKE DEPENDENT VAR COMPONENT
# y_random = input_fishnet_random$count
# # Set expected value (NOT SURE HOW TO DEAL WITH THIS!)
# E = rep(0,length(y_random));
# # set pop > 0 so we can use log(pop) as offset
# E[E < 1] = 0.01;
# 
# #### Compile Stan model
# bym2_stan_random = stan_model(file.path("STAN","BYM2.stan"), verbose = FALSE);
# 
# # prior_dist <- stan(file=file.path("STAN","BYM2.stan"),
# #                    iter=1000, warmup=0, chains=1,
# #                    seed=4838282, algorithm="Fixed_param")
# 
# #### Fit Stan model
# bym2_fit_random = sampling(bym2_stan_random, data=list(N = N,
#                                          N_edges = N_edges,
#                                          node1 = node1,
#                                          node2 = node2,
#                                          y = y_random,
#                                          E = E,
#                                          K = K,
#                                          x = x_random,
#                                          scaling_factor = scaling_factor), 
#                     control = list(adapt_delta = 0.97), 
#                     chains=4, warmup=2000, iter=4000, save_warmup=FALSE, verbose = TRUE);
# 
# #### Inspect Stan results
# print(bym2_fit_random, digits=3, pars=c("beta0", "rho", "logit_rho", "sigma", "betas[1]", "betas[2]",
#                                  "betas[3]", "betas[4]", "mu[1]", "mu[2]", "mu[3]", "phi[1]", "phi[2]", "phi[3]", "theta[1]", "theta[2]", "theta[3]"), probs=c(0.025, 0.5, 0.975))
# 
# posterior_random <- as.matrix(bym2_fit_random)
# 
# # betas
# mcmc_areas(posterior_random,
#            pars = c("betas[1]", "betas[2]", "betas[3]", "betas[4]"),
#            prob = 0.8)
# 
# ## plot predictions
# mu_random <- rstan::extract(bym2_fit_random, "mu")
# mu_plot_random <- data.frame(mu_random) %>% 
#   gather(parameter, estimate) %>% 
#   group_by(parameter) %>% 
#   summarise(mean = quantile(estimate, probs=0.5),
#             high = quantile(estimate, probs=0.9),
#             low  = quantile(estimate, probs=0.1)) %>% 
#   mutate(parameter = as.numeric(str_remove(parameter, "mu."))) %>% 
#   dplyr::left_join(., data.frame("parameter" = seq(1:length(y_random)), 
#                                  obs = y_random), by="parameter")
# 
# ggplot(mu_plot_random) +
#   geom_point(aes(x = obs, y = mean)) +
#   geom_point(aes(x = obs, y = high), color = "blue") +
#   geom_point(aes(x = obs, y = low), color = "red") +
#   scale_x_continuous(limits = c(0,max(mu_plot_random$high))) +
#   scale_y_continuous(limits = c(0,max(mu_plot_random$high))) +
#   labs(x = "Observed", y = "Predicted") +
#   geom_hline(yintercept = mean(y)) +
#   geom_abline() +
#   coord_fixed() +
#   theme_bw() +
#   theme(aspect.ratio=1)
# 
# 
# xx <- data.frame(posterior_random)
# xy <- data.frame(model = "xx",
#            rho = xx$rho,
#            logit_rho = xx$logit_rho,
#            betas_1 = xx$betas.1.,
#            betas_2 = xx$betas.2.,
#            betas_3 = xx$betas.3.,
#            betas_4 = xx$betas.4.,
           # sigma = xx$sigma)




