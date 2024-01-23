
   /* ... function declarations and definitions ... */
functions {
  //---------------------------------------------  
      vector sir(   real t, 
                vector y, 
                vector theta, 
                array[] real x_r, 
                array[] int x_i
            ) {

      vector[3] dy_dt;

      real S = y[1];
      real I = y[2];
      real R = y[3];   // constant real
      real N = x_i[1];  //constant int
      
      real beta = theta[1];
      real gamma = theta[2];
      
      dy_dt[1] = -beta * I * S / N;
      dy_dt[2] =  beta * I * S / N - gamma * I;
      dy_dt[3] =  gamma * I;
      
      return dy_dt;
  }
  //-----------------------------------------------
}

data {
  int<lower=1> n_days;
  vector[3] y0;
  real t0;
  array[n_days] real ts;
  int N;
  array[n_days] int cases;
}

transformed data {
  real x_r[0];
  int x_i[1] = { N };
}

parameters {
  real<lower=0> gamma;
  real<lower=0> beta;
  real<lower=0> phi_inv;
}

transformed parameters{
  array[n_days] vector[3] y;
  matrix[n_days, 3] y_matrix;
  real phi = 1. / phi_inv;
  {
    vector[2] theta;
    theta[1] = beta;
    theta[2] = gamma;

    y = ode_rk45(sir, y0, t0, ts, theta, x_r, x_i);
  }

  for (i in 1:n_days){
    y_matrix[i] = y[i]';
  }
}

model {
  //priors
  beta ~ normal(2, 1); //truncated at 0
  gamma ~ normal(0.4, 0.5); //truncated at 0
  phi_inv ~ exponential(5);
  
  //sampling distribution
  //col(matrix x, int n) - The n-th column of matrix x. Here the number of infected people
  cases ~ neg_binomial_2( col(y_matrix,2), phi);
}

generated quantities {
  real R0 = beta / gamma;
  real recovery_time = 1 / gamma;
  real pred_cases[n_days];
  pred_cases = neg_binomial_2_rng(col(y_matrix, 2) + 1e-5, phi);
}
