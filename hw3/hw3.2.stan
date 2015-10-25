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
parameters {
  real freq_effect[N];
  real sham_effect[N];
  real mu;
  real eta;
  real<lower=0> sigma_e;
  real<lower=0> sigma_s;
}
model {
  sham_effect ~ normal(eta,sigma_s);
  freq_effect ~ normal(mu,sigma_e);
  for (i in 1:N){
    mean_e[i] ~ cauchy(freq_effect[i],se_e[i]);
    mean_s[i] ~ cauchy(sham_effect[i],se_s[i]);
  }
}
generated quantities{
  real diff;
  diff <- mu - eta;
}