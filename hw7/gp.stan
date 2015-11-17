data{
  int N;
  int M;
  int<lower=0, upper=1> gender[N];
  int<lower=0> numKnown[N];
  int<lower=0> total[N];
  real age[N];
}
transformed data{
  cov_matrix[N] sigma;
  vector[N] mu;
  real sq_diff;
  for(i in 1:N) mu[i] <- 0;
  for(i in 1:N){
    for(j in 1:N){
      sigma[i,j] <- exp(-pow(age[i]-age[j],2)) + 
                      exp(-pow(gender[i]-gender[j],2)) + 
                      if_else(i==j,0.1,0.0);
    }
  }
}
parameters{
  vector[N] theta;
}
model{
  theta ~ multi_normal(mu,sigma);
  for(i in 1:N)
    numKnown[i] ~ binomial_logit(total[i],theta[i]);
}
generated quantities{
  vector[N] prob;
  for(i in 1:N)
    prob[i] <- inv_logit(theta[i]);
}