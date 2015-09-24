data {
	int N;
	vector[N] cPost;
	vector[N] cPre;
	vector[N] ePost;
	vector[N] ePre;
}
transformed data {
	vector[N] cdiff;
	vector[N] ediff;

	cdiff <- cPost - cPre;
	ediff <- ePost - ePre;
}
parameters {
	real meanC;
	real meanE;
	real<lower=0> seC;
	real<lower=0> seE;
}
model {
	for(i in 1:N)
		cdiff[i] ~ student_t(N, meanC, seC);
	for(i in 1:N)
		ediff[i] ~ student_t(N, meanE, seE);
}
generated quantities {
	vector[N] tdiff;
	tdiff <- ediff - cdiff;
}
