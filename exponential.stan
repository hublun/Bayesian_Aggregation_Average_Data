data{
    int LEN;
    vector[LEN] y;
}

parameters{
    real<lower=0> lambda;
}

model{
    real alpha;
    real beta;
    alpha = 1.0;
    beta = 1.0;

    lambda ~ gamma(alpha, beta);
    y ~ exponential(lambda);
}

generated quantities {
   /* ... declarations ... statements ... */
   real pred;
   pred = exponential_rng(lambda);
}
