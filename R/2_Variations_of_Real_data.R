library("maptools");
library("spdep");
library("rgdal")
library("rstan");
library("bayesplot")

# source functions
source(file.path("R","Functions_nb_data.R"))
source(file.path("R","Functions_archaeo_BYM.R"))

#### From PASS_data_prep.R
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
