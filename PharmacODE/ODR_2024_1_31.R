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
#list.dirs()
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
seir_model <- stan_model("./PharmacODE/seir_ode_final.stan")
#---------------------------------------------------------------------------------------------------------
cases <- df_swiss$report_dt

# times
n_days <- length(cases)
t <- seq(1, n_days, by = 1)
t0 = 0
ts <- t

date_switch <- "2020-03-13" # date of introduction of control measures
tswitch <- df_swiss %>% 
  filter(date < date_switch) %>% 
  nrow() + 1 # convert time to number

data_seir <- list(n_days = n_days, t0 = t0, ts = ts, tswitch = tswitch, N = N, cases = cases)

date_survey_left <- "2020-05-04"
date_survey_right <- "2020-05-07"
t_survey_start <- df_swiss %>% filter(date < date_survey_left) %>% nrow() + 1 # convert time to number
t_survey_end <- df_swiss %>% filter(date < date_survey_right) %>% nrow() + 1 # convert time to number
n_infected_survey <-  83
n_tested_survey <-  775
# add these data to the data given to stan
data_seir_survey <- c(data_seir, list(t_survey_start = t_survey_start, 
                                      t_survey_end = t_survey_end,
                                      n_infected_survey = n_infected_survey,
                                      n_tested_survey = n_tested_survey))
#-----------------------------------------------------------------------------
fit_seir <- sampling(seir_model, 
                    data_seir_survey, 
                    iter=1000,
                    seed = 0)
#-----------------------------------------------------------------------------
check_hmc_diagnostics(fit_seir)
pars = c("beta", "gamma", "phi", "alpha", "p_rep", "eta","nu","xi","lp__")
summary(fit_seir, pars=pars)
traceplot(fit_seir, pars=pars)
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
fit_seir %>% 
  spread_draws(alpha) %>% 
  filter(.chain == 2) %>% 
  ggplot() +
  geom_histogram(mapping = aes(x = alpha), fill=c_posterior, color=c_dark)
#-------------------------------------------------------------------------
