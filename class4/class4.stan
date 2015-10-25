data {
  int N;
  vector[N] x;
  vector[N] y;
}
parameters{
  real a;
  real b;
  real<lower=0> sigma_y;
}
model{
  vector[N] y_pred;
  for (n in 1:N){
    y_pred[n] <- a + b*x[n];
  }
  y ~ normal(y_pred,sigma_y);
}