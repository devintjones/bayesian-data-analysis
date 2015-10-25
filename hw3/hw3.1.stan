data {
  int N;
  real y[N];
}
parameters {
  real<lower=0,upper=1> theta;
}
model {
  y ~ cauchy(theta,1);
}
