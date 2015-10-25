data{
  int J;
}
parameters{
  real alpha;
  real beta;
  real<lower=0,upper=1> x[J];
}
transformed parameters{
  real theta[J];
  for (i in 1:J)
    theta[i] <- inv_logit(alpha + beta * x[i]);
}
model {
  int y[J];
  int n[J];

  for (i in 1:J)
    n[i]  ~ poisson(5);

  beta ~ student_t(4,0,1);
  alpha ~ student_t(4,0,2);

  for (i in 1:J)
    y[i] ~ binomial(n[i],theta[i]);	
}