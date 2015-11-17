data{
  int N;
  int M;
  matrix[N,M] X;
  real<lower=0,upper=1> r[N];
  int idx[N]; // survey index
  int basen[N]; // number requested to survey
}
transformed data{
  real<lower=0> V[N];
  for(i in 1:N)
    V[i] <- (r[i] * (1 - r[i] ))/basen[i];
}
parameters{
  vector[M] betas;
  real<lower=0> sigma;
  real<lower=0> tau;
  vector[39] alpha;
}
model{
  alpha ~ normal(0,tau);

  for(i in 1:N)
    r[i] ~ normal(X[i] * betas + alpha[idx[i]],
      sigma + V[i] );
}