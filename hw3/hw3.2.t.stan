data {
  int N;
  int N_e[N];
  int N_s[N];
  real mean_e[N];
  real mean_s[N];
  real<lower=0> se_e[N];
  real<lower=0> se_s[N];
  real freq[N];
}
transformed data{
  real log_y[N];
  for (i in 1:N)
    log_y[i] <- log(mean_e[i]);
}
parameters {
  real freq_effect[N];
}
model {
  real mu;
  real<lower=0> sigma_e;
  mu <- 0;
  sigma_e <- 1;  
  freq_effect ~ normal(mu,sigma_e);
  for (i in 1:N){
    log_y[i] ~ normal(freq_effect[i],se_e[i]);
  }
}