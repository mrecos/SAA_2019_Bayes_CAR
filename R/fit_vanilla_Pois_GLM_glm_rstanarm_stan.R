library("maptools");
library("spdep");
library("rgdal")
library("rstan");
library("bayesplot")
library("pdp")

options(mc.cores = 6);

# source functions
source(file.path("R","nb_data_funs.R"))
source(file.path("R","archaeo_BYM_Functions.R"))

input_fishnet = sites_fishnet_stan

#### MAKE INDEPENDENT VAR COMPONENT
x = st_drop_geometry(input_fishnet) %>% 
  dplyr::select(mean_val1, mean_val2, mean_val3, mean_val4) %>% 
  as.matrix()
K = ncol(x)

#### MAKE DEPENDENT VAR COMPONENT
y = input_fishnet$count
E = rep(0,length(y));
# set pop > 0 so we can use log(pop) as offset
E[E < 1] = 0.01;


#### Compile Stan model
library(rstanarm)

# Estimate original model
glm_dat <- as.data.frame(cbind(y,x))
pairs(glm_dat)
glm1 <- glm(y ~ ., data = glm_dat, offset = log(E), family = poisson)
glm_pred <- predict(glm1,  type = "response") %>% 
  as.data.frame() %>% 
  mutate(obs = y) %>% 
  arrange(desc(obs))

glm1 %>% # the %>% operator is read as "and then"
  partial(pred.var = "mean_val4", inv.link = exp) %>%
  plotPartial(smooth = FALSE, lwd = 2)

pd <- partial(glm1, pred.var = c("mean_val1", "mean_val2"), , inv.link = exp)
pdp1 <- plotPartial(pd)
rwb <- colorRampPalette(c("red", "white", "blue"))
pdp2 <- plotPartial(pd, contour = TRUE, col.regions = rwb)
pdp3 <- plotPartial(pd, levelplot = FALSE, zlab = "cmedv", drape = TRUE,
                    colorkey = TRUE, screen = list(z = -20, x = -60))
# Figure 3
grid.arrange(pdp1, pdp2, pdp3, ncol = 3)


# Estimate Bayesian version with stan_glm
stan_glm1 <- stan_glm(y ~ ., offset = log(E),
                      data = glm_dat, family = poisson, 
                      prior = normal(0,1), prior_intercept = normal(0,1),
                      chains = 4, cores = 4, seed = 717)

#### Inspect Stan results
print(stan_glm1, digits=3, pars=c("beta0", "sigma", "betas[1]", "betas[2]",
                                 "betas[3]", "betas[4]", "mu[1]", "mu[2]", "mu[3]"),
      probs=c(0.025, 0.5, 0.975))

#### Plot
# traceplot(stan_glm1, pars = c("beta0"), inc_warmup = FALSE, nrow = 2)

glm_posterior <- as.matrix(stan_glm1)

# plot of parameter distributions
plot_title <- ggtitle("Posterior distributions",
                      "with medians and 80% intervals")
# betas
mcmc_areas(glm_posterior,
           pars = c("mean_val1", "mean_val2", "mean_val3", "mean_val4"),
           prob = 0.8) + 
  plot_title +
  scale_x_continuous(breaks = seq(-2,1,0.1))


# from stan BYM2 fit
post_compare <- cbind(glm_posterior, posterior[1:nrow(glm_posterior),
                                               c("betas[1]","betas[2]","betas[3]","betas[4]")])

mcmc_areas(post_compare,
           pars = c("betas[1]", "mean_val1",
                    "betas[2]", "mean_val2", 
                    "betas[3]", "mean_val3", 
                    "betas[4]", "mean_val4"),
           prob = 0.8) + plot_title +
  scale_x_continuous(breaks = seq(-2,1,0.1))


### rstanarm predict
glm_mu <- rstanarm_pred <- posterior_predict(stan_glm1, type = "response")

glm_mu_plot <- data.frame(glm_mu) %>%
  gather(parameter, estimate) %>%
  group_by(parameter) %>%
  summarise(mean = quantile(estimate, probs=0.5),
            high = quantile(estimate, probs=0.9),
            low  = quantile(estimate, probs=0.1)) %>%
  mutate(parameter = as.numeric(str_remove(parameter, "X"))) %>%
  dplyr::left_join(., data.frame("parameter" = seq(1:length(y)), obs = y), by="parameter") %>% 
  arrange(desc(obs))

ggplot(glm_mu_plot) +
  geom_jitter(aes(x = obs, y = mean), width = 0.2, height = 0.1) +
  # geom_point(aes(x = obs, y = high), color = "blue") +
  # geom_point(aes(x = obs, y = low), color = "red") +
  # scale_x_continuous(limits = c(0,max(mu_plot$high))) +
  # scale_y_continuous(limits = c(0,max(mu_plot$high))) +
  labs(x = "Observed", y = "Predicted") +
  geom_hline(yintercept = mean(y)) +
  geom_abline() +
  coord_fixed() +
  theme_bw() +
  theme(aspect.ratio=1)

#### pure stan version (betas distribution matches output from rstan)
stan_glm_stan = stan_model(file.path("STAN","pois_glm_MDH.stan"), verbose = FALSE);
stan_glm_fit = sampling(stan_glm_stan, data=list(N = N,
                                         y = y,
                                         E = E,
                                         K = K,
                                         x = x), 
                    control = list(adapt_delta = 0.97), 
                    chains=4, warmup=2000, iter=4000, save_warmup=FALSE, verbose = TRUE);

print(stan_glm_fit, digits=3, pars=c("beta0", "betas[1]", "betas[2]",
                                 "betas[3]", "betas[4]", "mu[1]", "mu[2]", "mu[3]"), probs=c(0.025, 0.5, 0.975))

stan_glm_posterior <- as.matrix(stan_glm_fit)

## plot predictions
mu <- rstan::extract(stan_glm_fit, "mu")
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


## map mcmc results
mu_plot_map <- mu_plot %>% 
  arrange(parameter) %>% 
  as.data.frame() %>% 
  bind_cols(.,input_fishnet) %>% 
  st_as_sf() %>% 
  arrange(desc(count))

# mapview(mu_plot_map, zcol = "mean") + mapview(sites)


