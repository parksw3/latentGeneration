data {
	int t_censor; // when was the data censored
	
	int npatient;
	
	int from[npatient];
	
	int<lower=0> dt_start_exp[npatient]; // first day of exposure
	int<lower=0> dt_end_exp[npatient]; // last day of exposure
	int<lower=0> dt_onset[npatient]; // day of symptom onset
	
	int prior_only;
}

parameters {
	real<lower=0> popR;
	real<lower=0> sigmaR;
	vector<lower=0>[npatient] indR;
	
	real<lower=0> shape_gen;
	real<lower=0> rate_gen;
	real<lower=0> shape_inf;
	real<lower=0> rate_inf;
	real<lower=0> phi;
	
	vector<lower=0, upper=1>[npatient] inf; // infectious period
}

transformed parameters {
	vector[npatient] dt_infected;
	vector[npatient] infectious_period;
	vector[npatient] gen;

	for (i in 1:npatient) {
		dt_infected[i] = dt_start_exp[i] + inf[i] * (dt_end_exp[i] - dt_start_exp[i]);
		infectious_period[i] = dt_onset[i] - dt_infected[i];
	}
	
	for (i in 1:npatient) {
		if (from[i] != 0) {
			gen[i] = dt_infected[i] - dt_infected[from[i]];
		} else {
			gen[i] = 0;
		}
	}
}

model {
	popR ~ gamma(2, 2);
	sigmaR ~ gamma(2, 4);
	
	shape_gen ~ gamma(4, 2);
	rate_gen ~ gamma(4, 20);
	
	shape_inf ~ gamma(4, 2);
	rate_inf ~ gamma(4, 10);
	
	inf ~ uniform(0, 1);
	
	phi ~ gamma(4, 20);
	
	indR ~ normal(popR, sigmaR);
	
	if (prior_only==0) {
		
		// infectious period
		for (i in 1:npatient) {
			target += gamma_lpdf(infectious_period[i] | shape_inf, rate_inf);
		}
	
		// generation time
		for (i in 1:npatient) {
			if (from[i] != 0) {
				target += gamma_lpdf(gen[i] | shape_gen, rate_gen); 
			}
		}
	
		// big R
		for (i in 1:npatient) {
			if (from[i] != 0) {
				target += log(indR[from[i]]) - phi * dt_infected[i];
			} else {
				target += log(popR) - phi * dt_infected[i];
			}
		}
	
		// discounting
		for (i in 1:npatient) {
			target += - indR[i] * (1/phi) * (exp(- phi * dt_infected[i]) - exp(- phi * t_censor)) * gamma_cdf(t_censor - dt_infected[i], shape_gen, rate_gen);
		}
	}
}
