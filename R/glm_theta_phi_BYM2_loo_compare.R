library(loo)

log_lik_glm <- extract_log_lik(glm_stan_fit, merge_chains = FALSE)
rel_n_eff_glm <- relative_eff(exp(log_lik_glm))
glm_loo <- loo(log_lik_glm, r_eff = rel_n_eff_glm, cores = 2)

log_lik_theta <- extract_log_lik(theta_stan_fit, merge_chains = FALSE)
rel_n_eff_theta <- relative_eff(exp(log_lik_theta))
theta_loo <- loo(log_lik_theta, r_eff = rel_n_eff_theta, cores = 2)

log_lik_phi <- extract_log_lik(phi_stan_fit, merge_chains = FALSE)
rel_n_eff_phi <- relative_eff(exp(log_lik_phi))
phi_loo <- loo(log_lik_phi, r_eff = rel_n_eff_phi, cores = 2)

log_lik_bym2 <- extract_log_lik(bym2_fit, merge_chains = FALSE)
rel_n_eff_bym2 <- relative_eff(exp(log_lik_bym2))
bym2_loo <- loo(log_lik_bym2, r_eff = rel_n_eff_bym2, cores = 2)

loo::compare(glm_loo, theta_loo, phi_loo, bym2_loo)
