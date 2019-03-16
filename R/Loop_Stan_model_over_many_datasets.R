library("LaCroixColoR")
library("pals")

models_list <- list(siteC_envC = siteC_envC,
                    siteC_envI = siteC_envI,
                    siteC_envR = siteC_envR,
                    siteC_envU = siteC_envU,
                    siteI_envC = siteI_envC,
                    siteI_envI = siteI_envI,
                    siteI_envR = siteI_envR,
                    siteI_envU = siteI_envU,
                    siteR_envC = siteR_envC,
                    siteR_envI = siteR_envI,
                    siteR_envR = siteR_envR,
                    siteR_envU = siteR_envU,
                    siteU_envC = siteU_envC,
                    siteU_envI = siteU_envI,
                    siteU_envR = siteU_envR,
                    siteU_envU = siteU_envU)
                    

rho_results <- data.frame()

for(i in seq_along(models_list)){
  
  cat("Running model",i,"of",length(models_list),"-",names(models_list[i]),"\n")
  input_fishent = models_list[[i]]
  
  # convert fishnet to spatial
  fishnet.sp <- as(input_fishent, "Spatial")
  
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
  x = st_drop_geometry(input_fishent) %>% 
    dplyr::select(mean_val1, mean_val2, mean_val3, mean_val4) %>% 
    as.matrix()
  K = ncol(x)
  
  #### MAKE DEPENDENT VAR COMPONENT
  y = input_fishent$count
  # Set expected value (NOT SURE HOW TO DEAL WITH THIS!)
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
                      chains=4, warmup=200, iter=4000, save_warmup=FALSE, verbose = FALSE);
  
  
  #### Plot
  posterior <- as.matrix(bym2_fit)
  
  post_df <- data.frame(posterior)
  rho_results <- rbind(rho_results, data.frame(model = names(models_list[i]),
                                               rho = post_df$rho,
                                               logit_rho = post_df$logit_rho,
                                               betas_1   = post_df$betas.1.,
                                               betas_2   = post_df$betas.2.,
                                               betas_3   = post_df$betas.3.,
                                               betas_4   = post_df$betas.4.,
                                               sigma     = post_df$sigma))
  
}

rho_results <- rho_results %>% group_by(model) %>%
  mutate(med = median(logit_rho),
         site = str_sub(model, 1, 5),
         env  = str_sub(model, 7, 10)) 

# saveRDS(rho_results, "rho_results.rdata")

ggplot(rho_results, aes(x=logit_rho)) +
  # geom_density(fill = "skyblue", color = "skyblue4") +
  geom_density(aes(fill=model)) +
  geom_vline(xintercept = 0, colour = 'gray30', linetype = "dashed") +
  geom_vline(aes(xintercept = med, group = model), colour = 'black') +
  scale_x_continuous(limits = c(-10,10))+
  facet_grid(env~site) +
  labs(x="") +
  # scale_fill_manual(values = lacroix_palette("KiwiSandia",n = 4)) +
  scale_fill_manual(values = stepped(n = 16)) +
  theme_bw() +
  theme(
    legend.position = "none",
    axis.title.y=element_blank(),
    axis.text.y=element_blank(),
    axis.ticks.y=element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    strip.background = element_blank()
    # strip.text.x = element_blank(),
    # strip.text.y = element_blank()
  )


# ggplot(rho_results, aes(x=rho)) +
#   geom_density(fill = "skyblue", color = "skyblue4") +
#   geom_vline(xintercept = 0.5, colour = 'gray30') +
#   scale_x_continuous(limits = c(0,1))+
#   facet_grid(site~env) +
#   theme_bw()
# 
# 
# ggplot(rho_results, aes(x=sigma)) +
#   geom_density(fill = "skyblue", color = "skyblue4") +
#   # geom_vline(xintercept = 0.5, colour = 'gray30') +
#   # scale_x_continuous(limits = c(0,1))+
#   facet_grid(site~env) +
#   theme_bw()
# 
# 
# ggplot(rho_results, aes(x=betas_1)) +
#   geom_density(fill = "skyblue", color = "skyblue4") +
#   # geom_vline(xintercept = 0.5, colour = 'gray30') +
#   scale_x_continuous(limits = c(-2,2))+
#   facet_grid(site~env) +
#   theme_bw()
# 

