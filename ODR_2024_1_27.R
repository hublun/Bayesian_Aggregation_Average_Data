rm(list=ls())
library(rstan)
options(mc.cores = parallel::detectCores())
#--------------------------------------------
# color palette
c_light <- c("#fee0d2")
c_mid <- c("#fc9272")
c_dark <- c("#de2d26")
c_simu <- "chartreuse3"
c_posterior = "orange"
c_prior = "aquamarine2"

set.seed(1) # for reproductibility
#------------------------------------------

getwd()
list.dirs()
setwd("~/Bayesian_Aggregation_Average_Data")
#--------------------------------------------
library(tidybayes)
library(gridExtra)
library(tidyverse)
df_swiss <- read_csv("~/Bayesian_Aggregation_Average_Data/disease_transmission_workflow/data/swiss_agg_data.csv")
#-------------------------------------------------------------------------------------------------------------------
df_swiss %>% 
  ggplot() + 
  geom_bar(mapping = aes(x = date, y = report_dt), fill = c_mid, color = c_dark, stat = "identity") +
  labs(y="Number of reported cases")
#--------------------------------------------------------------------------------------------------------------------
# Swiss population
N <- 8.57E6;

#initial conditions
i0 <- 1
s0 <- N - i0
r0 <- 0

y0 = c(S = s0, I = i0, R = r0) # initial state

#---------------------------------------------------------------------------------------------------------
sir_model <- stan_model("./disease_transmission_workflow/stan_models/models_swiss/seir_ode.stan")
#---------------------------------------------------------------------------------------------------------
cases <- df_swiss$report_dt

# times
n_days <- length(cases)
t <- seq(1, n_days, by = 1)
t0 = 0
ts <- t

data_seir <- list(n_days = n_days, t0 = t0, ts = ts, N = N, cases = cases)
#-----------------------------------------------------------------------------
fit_seir <- sampling(sir_model, 
                    data_seir, 
                    iter=1000,
                    seed = 0)
#-----------------------------------------------------------------------------
check_hmc_diagnostics(fit_sir)
pars = c("beta", "gamma", "phi", "alpha", "p_rep", "lp__")
summary(fit_sir, pars=pars)
traceplot(fit_sir, pars=pars)
#----------------------------------------------------------------------------------------------------------------------------------
dim(as.data.frame(summary(fit_seir, pars = "pred_cases", probs = c(0.025, 0.05, 0.1, 0.5, 0.9, 0.95, 0.975))$summary))

smr_pred <- cbind(as.data.frame(summary(fit_seir, pars = "pred_cases", probs = c(0.025, 0.05, 0.1, 0.5, 0.9, 0.95, 0.975))$summary), 
                  t=1:(n_days-1), 
                  cases = cases[1:length(cases)-1]
              )

colnames(smr_pred)

colnames(smr_pred) <- make.names(colnames(smr_pred)) # to remove % in the col names

ggplot(smr_pred, mapping = aes(x = t)) +
  #geom_ribbon(aes(ymin = X2.5., ymax = X97.5.), fill = c_dark, ) +
  geom_ribbon(aes(ymin = X5., ymax = X95.), fill = c_posterior, alpha=0.35) +
  #geom_ribbon(aes(ymin = X10., ymax = X90.), fill = c_light) +
  geom_line(mapping = aes(x = t, y = X50.), color = c_posterior) +
  geom_point(mapping = aes(y = cases)) +
  labs(x = "Day", y = "Incidence")
#----------------------------------------------------------------------------------------------------------------------------------------------
fit_seir %>% 
  spread_draws(pred_cases[n_days]) %>% 
  left_join(tibble(cases = cases, n_days = 1:length(cases))) %>% 
  group_by(n_days, .chain) %>% 
  summarise(cases = mean(cases), pred_mean = mean(pred_cases), pred_9 = quantile(pred_cases, 0.95), pred_1 = quantile(pred_cases, 0.05)) %>% 
  ggplot(aes(x = n_days)) +
  geom_ribbon(aes(ymin = pred_1, ymax = pred_9), fill = c_posterior, alpha=0.35)+
  geom_line(mapping = aes(y=pred_mean), color = c_posterior)+
  geom_point(mapping = aes(y=cases), size=0.1)+
  facet_wrap(~.chain)
#-----------------------------------------------------------------------------------------------------------------------------------------------