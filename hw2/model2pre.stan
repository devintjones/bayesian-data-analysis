data {
	int J;
	vector[J] x_trials;
	int<lower=0,upper=1> y_trials[J];
}

parameters {
	real<lower=-5,upper=10> alpha;
	real<lower=-10,upper=40> beta;
	vector[2] b_prior;
}

model {
	vector[2] b_pre;
	matrix[2,2] sigma;

	// prior on alpha, beta
	alpha ~ normal(0,4);
	beta ~ normal(10,100);

	sigma[1,1] <- 4;
	sigma[2,2] <- 100;
	sigma[1,2] <- 10;
	sigma[2,1] <- 10;


	b_pre[1] <- alpha;
	b_pre[2] <- beta;
	b_prior ~ multi_normal(b_pre,sigma);

	// posterior
	//y_trials ~ bernoulli_logit(b_prior[1] + b_prior[2] * x_trials);
}
