data{
  int N;
  int M;
	matrix[N,M] X;
	int y[N];
}
parameters{
  vector[M] betas;
}
model{
  for(i in 1:M)
    betas[i] ~ cauchy(0,2.5);
  for(i in 1:N)
    y[i] ~ poisson(exp(X[i]*betas));
}