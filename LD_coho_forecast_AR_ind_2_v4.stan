data{
	//Number of years (total, includes several missing years for some stocks)
	int <lower=0> n_year;

	//Number of total populations
	int <lower=0> n_pop;

	//Number of populations with return data
	int <lower=0> n_pop_tot;

	//Which populations possess return data
	int <lower=0> pop_tot[n_pop_tot];

	//Length of the return data vectors
	int <lower=0> n_tot;

	//Vectors of all return data across all populations
	real <lower=0> tot_dat [n_tot];

	//Vectors of the indices identifying which years are those with non-NA data for the return data
	int <lower=0> tot_true [n_tot];

	//Paired vectors of slice points indicating the beginning, and end of the data for a particular population
	int <lower=0> slice_tot_start [n_pop_tot];
	int <lower=0> slice_tot_end [n_pop_tot];
	
}
transformed data{
  real log_tot_dat [n_tot]; //log transformed return data
  log_tot_dat[1:n_tot] = log(tot_dat[1:n_tot]);
}
parameters {//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  vector<lower=0> [n_pop] sigma_tot;
  vector<lower=0,upper=1> [n_pop] p_proc;
  vector<lower=0,upper=1> [n_pop] p_obs;
  //pop-specific observation sd
  //vector <lower=0> [n_pop] sigma_obs;
	//Population-specific process error sd
	//vector <lower=0> [n_pop] sigma_proc;
	//autocorrelation coefficients
	vector <lower=-0.99999, upper = 0.99999> [n_pop] phi;
	//process residuals
	matrix [n_year-1,n_pop] proc_dev;
	// log of initial adult abundance residual
  vector [n_pop] init;
  //mean log abundance
  vector [n_pop] mu;
}
transformed parameters {//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// State estimates of log adult returns
  matrix [n_year,n_pop] log_adult_est;
  //pop-specific observation sd
  vector <lower=0> [n_pop] sigma_obs;
	//Population-specific process error sd
	vector <lower=0> [n_pop] sigma_proc;
  // State estimates of adult returns
  matrix <lower=0> [n_year,n_pop] adult_est;
	//Population-specific random walk of adult returns
	for(i in 1:n_pop){
	  sigma_proc[i] = sqrt(square(sigma_tot[i]) * p_proc[i]);
	  sigma_obs[i] = sqrt(square(sigma_tot[i])* p_obs[i]);
		for(y in 1:n_year){
			if(y==1){
				log_adult_est[y,i] = mu[i] + init[i] * sqrt(square(sigma_proc[i])/(1-square(phi[i])));
			}else{
				log_adult_est[y,i] = mu[i] + phi[i] * (log_adult_est[y-1,i] - mu[i]) + proc_dev[y-1,i] * sigma_proc[i];
			}
		}
	}
	//exponentiate adult returns from estimation in log space
	adult_est = exp(log_adult_est);
}
model {//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
 	
 	//Standardized Process deviations drawn from population-specific error distributions
 	for(i in 1:n_pop){
 		proc_dev[,i] ~ std_normal();
 	}

  //residuals from long term mean for first state
  init ~ std_normal();

 	//Vague prior for initial log abundance
 	mu ~ normal(8,4); 
  
  sigma_tot ~ std_normal();
  p_proc ~ beta(1,1);
  p_obs ~ beta(1,1);
    
  // Prior proc error sd
 	//sigma_obs ~ std_normal();

 	// Prior obs error sd
 	//sigma_proc ~ std_normal();
 	
 	//prior on autocorrelation coefficients
 	phi ~ normal(0,10);
 	
	//Likelihood
	for(i in 1:n_pop_tot){
		//tot_dat[slice_tot_start[i]:slice_tot_end[i]] ~ lognormal(log_adult_est[tot_true[slice_tot_start[i]:slice_tot_end[i]], pop_tot[i]], sigma_obs[i]);
		log_tot_dat[slice_tot_start[i]:slice_tot_end[i]] ~ normal(log_adult_est[tot_true[slice_tot_start[i]:slice_tot_end[i]], pop_tot[i]], sigma_obs[i]);
	}

}
generated quantities {//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //Vector of adult return forecasts
  vector[n_pop] adult_pred;

  for(i in 1:n_pop){
		adult_pred[i] = exp(log_adult_est[n_year, i]);
	}
}

