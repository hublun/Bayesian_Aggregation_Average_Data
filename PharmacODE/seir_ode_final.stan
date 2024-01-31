functions {
  //================================================================
    real switch_eta(real t, 
                    real t1, 
                    real eta, 
                    real nu, 
                    real xi) {
    return(eta + (1 - eta) / (1 + exp(xi * (t - t1 - nu))));
  }
  //---------------------------------
  vector seir( real t, 
              vector y, 
              vector theta, 
              array[] real x_r, 
              array[] int x_i) {
      
      vector[4] dy_dt;
      
      real N = x_i[1];
      real tswitch = x_r[1];

      real beta = theta[1];
      real gamma = theta[2];
      real alpha = theta[3];

      real eta = theta[4];
      real nu = theta[5];
      real xi = theta[6];

      real i0 = theta[7];
      real e0 = theta[8];
      real forcing = switch_eta(t,tswitch,eta,nu,xi); // switch function
      real beta_eff = beta * forcing; // effective beta taking into consideraton of control policy

      array[4] real init = {N - i0 - e0, e0, i0, 0};
      
      real S = y[1] + init[1];
      real E = y[2] + init[2];
      real I = y[3] + init[3];
      real R = y[4] + init[4];      
      
      
      dy_dt[1] = -beta_eff * I * S / N;
      dy_dt[2] = beta_eff * I * S / N - alpha * E;
      dy_dt[3] = alpha * E - gamma * I;
      dy_dt[4] =  gamma * I;
      
      return dy_dt;
  }
}
data {
  int<lower=1> n_days;
  //vector[3] y0;
  real t0;
  real tswitch; // date of introduction of control measures
  array[n_days] real ts;
  int N;
  array[n_days] int cases;
  
  int t_survey_start; // antibody survey data
  int t_survey_end;
  int n_infected_survey;
  int n_tested_survey;
}
transformed data {
  real x_r[1] = { tswitch };
  int x_i[1] = { N };
}
parameters {
  real<lower=0> gamma;
  real<lower=0> beta;
  real<lower=0> alpha; // rate of exposure

  real<lower=0, upper=1> eta; // reduction in transmission due to control measures (in proportion of beta)
  real<lower=0> nu; // shift of quarantine implementation (strictly positive as it can only occur after tswitch)
  real<lower=0,upper=1> xi_raw; // slope of quarantine implementation (strictly positive as the logistic must be downward)  

  real<lower=0> phi_inv;
  real<lower=0, upper=1> p_rep; # proportion of infected reported
  real<lower=0> e0; // number of initialll exposed
  real<lower=0> i0; // number of initially infected
}
transformed parameters{
  array[n_days] vector[4] y;
  array[n_days-1] real incidence;
  real phi = 1. / phi_inv;
  real xi = xi_raw + 0.5;
  real<lower=0, upper=1> p_infected_survey; //proportion of people having been infected at week 5 (between 4 and 7 may)
  
  {
    vector[8] theta;

    theta[1] = beta;
    theta[2] = gamma;
    theta[3] = alpha;
    theta[4] = eta;
    theta[5] = nu;
    theta[6] = xi;
    theta[7] = i0;
    theta[8] = e0;

    y = ode_rk45(seir, rep_vector(0.0, 4), t0, ts, theta, x_r, x_i);
  }
  for (i in 1:n_days-1){
    incidence[i] =  -(y[i+1, 2] - y[i, 2] + y[i+1,1] - y[i,1]) * p_rep; //S(t+1) - S(t) + E(t+1) - S(t)
  }
  // mean number of infected + recovered people during week 5
  p_infected_survey = mean(to_vector(y[t_survey_start:t_survey_end, 4])) / N;
  //=============================
}
model {
  //priors
  beta ~ normal(2, 1);
  gamma ~ normal(0.4, 0.5);
  phi_inv ~ exponential(5);
  alpha ~ normal(0.4, 0.5);
  
  eta ~ beta(2.5, 4);
  nu ~ exponential(1./5);
  xi_raw ~ beta(1, 1);

  p_rep ~ beta(1, 2);
  i0 ~ normal(0, 10);
  e0 ~ normal(0, 10);
  
 //likelihood
  n_infected_survey ~ binomial(n_tested_survey, p_infected_survey); // we fit the survey data to our latent parameter
  cases[1:(n_days-1)] ~ neg_binomial_2(incidence, phi);
  //=======================================================
}
generated quantities {
  real R0 = beta / gamma;
  real Reff[n_days]; // R0 but taking into account environmental changes
  real recovery_time = 1 / gamma;
  array[n_days-1] real pred_cases;
  pred_cases = neg_binomial_2_rng(incidence, phi);
  
  for (i in 1:n_days)
    Reff[i] = switch_eta(i, tswitch, eta, nu, xi) * beta / gamma;
  //-----------------------------------------
}

