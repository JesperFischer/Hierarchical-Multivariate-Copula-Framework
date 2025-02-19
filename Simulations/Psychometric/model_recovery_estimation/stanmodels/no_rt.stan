functions {
  
  real psychometric(real x,real threshold, real slope,real lapse){
  
    return(lapse+(1-2*lapse)*(1 / (1+exp(-slope * (x - threshold)))));
      
  }
    vector psychometric_vec(vector x,real threshold, real slope,real lapse){
  
    return(lapse+(1-2*lapse).*(1 / (1+exp(-slope .* (x - threshold)))));
      
  }
  
  



}

data {
  int<lower=0> N;
  int<lower=0> S;
  
  array[N] int binom_y;
  vector[N] RT;
  vector[N] x;
  vector[S] minRT;
  array[N] int S_id; // Vector of binary responses
  array[S] int starts;
  array[S] int ends;
  





}
transformed data{
  int P = 3;
}

parameters {
  vector[P] gm;  // Group means 

  vector<lower = 0>[P]  tau_u;   // Between participant scales

  matrix[P, S] z_expo;    // Participant deviation from the group means


  
}

transformed parameters{
  
  matrix[P, S] param;
  
  // adding the participant deviation to the group means
  for(p in 1:P){
    param[p,]= to_row_vector(gm[p] + z_expo[p,] * tau_u[p]);
  }
  
  
  row_vector[S] threshold = (param[1,]);
  row_vector[S] slope = exp(param[2,]);
  row_vector[S] lapse = inv_logit(param[3,]) / 2;
  
  
  vector[size(x)] expect;
    
  // for (n in 1:N) {
  //    expect[n] = psychometric(x[n],threshold[S_id[n]],(slope[S_id[n]]),(lapse[S_id[n]])); 
  // }
  
  for (s in 1:S) {
    expect[starts[s]:ends[s]] = psychometric_vec(x[starts[s]:ends[s]],threshold[s],(slope[s]),(lapse[s])); 
  }

  
}


model {
  
  gm[1] ~ normal(0,5);
  gm[2] ~ normal(-1.5,0.5);
  gm[3] ~ normal(-3,1);

  to_vector(z_expo) ~ std_normal();

  tau_u[1] ~ normal(5,2);
  tau_u[2] ~ normal(0.6,0.2);
  tau_u[3] ~ normal(0.6,0.2);

  
  // for (n in 1:N) {
  //   target += lognormal_lpdf(RT[n] - rt_ndt[S_id[n]]| rt_int[S_id[n]] + rt_beta[S_id[n]] * entropy(expect[n]), rt_sd[S_id[n]]);
  // }

  for(s in 1:S){

    target += bernoulli_lpmf(binom_y[starts[s]:ends[s]] |  expect[starts[s]:ends[s]]);

  }

}

generated quantities {

  real mu_threshold = gm[1];
  real mu_slope = gm[2];
  real mu_lapse = gm[3];

  real tau_threshold = tau_u[1];
  real tau_slope = tau_u[2];
  real tau_lapse = tau_u[3];



  vector[N] log_lik_bin = rep_vector(0,N);
  vector[N] log_lik = rep_vector(0,N);



  for(n in 1:N){


    log_lik_bin[n] = binomial_lpmf(binom_y[n] | 1, expect[n]);

    log_lik[n] = log_lik_bin[n];
  }

}
