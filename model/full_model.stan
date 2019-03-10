data {
	int npatient;
	int n_missing_onset;
	int n_missing_end_exp;
	
	int from[npatient];
	
	int<lower=-1> dt_start_exp[npatient]; // first day of exposure
	int<lower=-1> dt_end_exp[npatient]; // last day of exposure
	int<lower=-1> dt_onset[npatient]; // day of symptom onset
	int<lower=-1> dt_report[npatient];
	int missing_onset[npatient];
	int missing_end_exp[npatient];
	
	
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
	
	real<lower=0> rate_report;
	
	vector<lower=0, upper=1>[npatient] inf; // unscaled infectious period
	vector<lower=0, upper=1>[n_missing_onset] rep; // unscaled reporting period
}

transformed parameters {
	vector[npatient] dt_infected;
	vector[npatient] infectious_period;
	vector[npatient] gen;
	vector[npatient] dt_onset_new;
	
	for (i in 1:npatient) {
		if (missing_onset[i] == 0) {
			dt_onset_new[i] = dt_onset[i];
		} else {
			dt_onset_new[i] = dt_start_exp[i] + rep[missing_onset[i]] * (dt_report[i] - dt_start_exp[i]);
		}
	}

	for (i in 1:npatient) {
		if (missing_onset[i] == 0) {
			dt_infected[i] = dt_start_exp[i] + inf[i] * (dt_end_exp[i] - dt_start_exp[i]);
		} else {
			dt_infected[i] = dt_start_exp[i] + inf[i] * (dt_onset_new[i] - dt_start_exp[i]);
		}
		
		infectious_period[i] = dt_onset_new[i] - dt_infected[i];
	}
	
	for (i in 1:npatient) {
		if (from[i] != -1) {
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
	
	rate_report ~ exponential(1);
	
	inf ~ uniform(0, 1);
	rep ~ uniform(0, 1);
	
	phi ~ gamma(4, 20);
	
	indR ~ normal(popR, sigmaR);
	
	if (prior_only == 0) {
		// reporting kernel
		for (i in 1:npatient) {
			target += exponential_lpdf(dt_report[i] - dt_onset_new[i] | rate_report);
		}
		
		// infectious period
		for (i in 1:npatient) {
			target += gamma_lpdf(infectious_period[i] | shape_inf, rate_inf);
		}
		
		// generation time
		for (i in 1:npatient) {
			if (from[i] != -1) {
				target += gamma_lpdf(gen[i] | shape_gen, rate_gen); 
			}
		}
		
		// big R
		for (i in 1:npatient) {
			if (from[i] != -1) {
				target += log(indR[from[i]]) - phi * dt_infected[i];
			} else {
				target += log(popR) - phi * dt_infected[i];
			}
		}
		
		target += - popR * (1/phi);
		
		// discounting
		for (i in 1:npatient) {
			target += - indR[i] * exp(-phi * dt_infected[i]) * pow(rate_gen, shape_gen)/pow(rate_gen + phi, shape_gen);
		}
	}
}
