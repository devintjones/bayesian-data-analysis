data {
	int N;
	int K;
	vector[N] cPost;
	vector[N] cPre;
	vector[N] ePost;
	vector[N] ePre;
	vector[N] city;
	vector[N] ones;
	matrix[N,K] grade;

}
transformed data {

	matrix[N,K+1] temp1;
	matrix[N,K+1] temp2;

	matrix[N,K+2] eData;
	matrix[N,K+2] cData;

	vector[N] cdiff;
	vector[N] ediff;
	
	cdiff <- cPost - cPre;
	ediff <- ePost - ePre;

	temp1 <- append_col(ones,grade);
	cData <- append_col(temp1,city);
	
	temp2 <- append_col(ones,grade);
	eData <- append_col(temp2,city);

}
parameters {
	real<lower=0> eVar;
	real<lower=0> cVar;
	real<lower=0> tVar;

	vector[K+2] cBeta;
	vector[K+2] eBeta;

	vector[N] tdiff;
}
model {
	for(i in 1:N)
		cdiff[i] ~ normal(cData[1] * cBeta[1] +
						  cData[2] * cBeta[2] +
						  cData[3] * cBeta[3] +
						  cData[4] * cBeta[4] +
						  cData[5] * cBeta[5], cVar);
	for(i in 1:N)
		ediff[i] ~ normal(eData[1] * eBeta[1] +
						  eData[2] * eBeta[2] +
						  eData[3] * eBeta[3] +
						  eData[4] * eBeta[4] +
						  eData[5] * eBeta[5], eVar);
	for(i in 1:N)
		tdiff[i] ~ normal(ediff[i] - cdiff[i], tVar);

}
