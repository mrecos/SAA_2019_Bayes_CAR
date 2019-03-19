library(tidyverse)

data.frame(values = rnorm(1000000,0,2.5),
           param  = "betas") %>% 
  ggplot(. ,aes(x = values)) +
  geom_density(fill = "skyblue") +
  geom_vline(xintercept = 0, color = "firebrick") +
  theme_minimal() +
  labs(x='',y='') +
  theme(
    axis.title.y=element_blank(),
    axis.text.y=element_blank(),
    axis.ticks.y=element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),
    strip.text.x = element_blank()
  )
ggsave("c:/TEMP/betas.svg", width = 5, height = 5)


data.frame(values = rnorm(1000000,0,1),
           param  = "theta") %>% 
  ggplot(. ,aes(x = values)) +
  geom_density(fill = "skyblue") +
  geom_vline(xintercept = 0, color = "firebrick") +
  theme_minimal() +
  labs(x='',y='') +
  theme(
    axis.title.y=element_blank(),
    axis.text.y=element_blank(),
    axis.ticks.y=element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),
    strip.text.x = element_blank()
  )
ggsave("c:/TEMP/theta.svg", width = 5, height = 5)

data.frame(values = rnorm(1000000,0,5),
           param  = "sigma") %>% 
  ggplot(. ,aes(x = values)) +
  geom_density(fill = "skyblue") +
  geom_vline(xintercept = 0, color = "firebrick") +
  theme_minimal() +
  labs(x='',y='') +
  theme(
    axis.title.y=element_blank(),
    axis.text.y=element_blank(),
    axis.ticks.y=element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),
    strip.text.x = element_blank()
  )
ggsave("c:/TEMP/sigma.svg", width = 5, height = 5)


data.frame(values = rbeta(1000000,0.5,0.5),
           param  = "rho") %>% 
  ggplot(. ,aes(x = values)) +
  geom_density(fill = "skyblue") +
  geom_vline(xintercept = 0.5, color = "firebrick") +
  theme_minimal() +
  labs(x='',y='') +
  theme(
    axis.title.y=element_blank(),
    axis.text.y=element_blank(),
    axis.ticks.y=element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),
    strip.text.x = element_blank()
  )
ggsave("c:/TEMP/rho.svg", width = 5, height = 5)

### or facets for betas, theta, sigma
betas <- data.frame(values = rnorm(100000,0,2.5),
                    param  = "betas") 

theta <- data.frame(values = rnorm(100000,0,1),
                    param  = "theta")

sigma <- data.frame(values = rnorm(100000,0,5),
                    param  = "sigma") 

param_data <- rbind(betas, theta, sigma)

ggplot(param_data ,aes(x = values, group = param)) +
  geom_density(fill = "skyblue") +
  geom_vline(xintercept = 0, color = "firebrick") +
  theme_minimal() +
  facet_wrap(~param, ncol = 1) +
  labs(x='',y='') +
  theme(
    axis.title.y=element_blank(),
    axis.text.y=element_blank(),
    axis.ticks.y=element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),
    strip.text.x = element_blank()
  )

ggsave("c:/TEMP/betas_theta_sigma_params.svg", width = 3, height = 7)

