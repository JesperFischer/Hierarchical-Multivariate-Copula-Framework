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
  int P = 6;
}

parameters {
  vector[P] gm;  // Group means 

  vector<lower = 0>[P]  tau_u;   // Between participant scales

  matrix[P, S] z_expo;    // Participant deviation from the group means

  vector<lower=0, upper = minRT>[S] rt_ndt;
  
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
  row_vector[S] alpha = exp(param[4,]);
  row_vector[S] beta = inv_logit(param[5,]);
  row_vector[S] delta = (param[6,]);
  
  
  vector[size(x)] expect;
    
  
  for (s in 1:S) {
    expect[starts[s]:ends[s]] = psychometric_vec(x[starts[s]:ends[s]],threshold[s],(slope[s]),(lapse[s])); 
  }

  
}


model {
  
  gm[1] ~ normal(0,5);
  gm[2] ~ normal(1,1);
  gm[3] ~ normal(-3,1);
  gm[4] ~ normal(0.5,1);
  gm[5] ~ normal(0,0.5);
  gm[6] ~ normal(5,1);
  
  to_vector(z_expo) ~ std_normal();

  tau_u[1] ~ normal(5,2);
  tau_u[2] ~ normal(0.6,1);
  tau_u[3] ~ normal(0.6,1);
  tau_u[4] ~ normal(0.4,1);
  tau_u[5] ~ normal(0.4,1);
  tau_u[6] ~ normal(0.5,1);

  
  rt_ndt ~ normal(0.2,0.05);
  // for (n in 1:N) {
  //   target += lognormal_lpdf(RT[n] - rt_ndt[S_id[n]]| rt_int[S_id[n]] + rt_beta[S_id[n]] * entropy(expect[n]), rt_sd[S_id[n]]);
  // }



  for(n in 1:N){
    int c = binom_y[n];
    if(c == 1){
      RT[n] ~ wiener(alpha[S_id[n]],rt_ndt[S_id[n]],beta[S_id[n]],(expect[n]-(1-expect[n]))*delta[S_id[n]]);
    }else if(c == 0){
      RT[n] ~ wiener(alpha[S_id[n]],rt_ndt[S_id[n]],1-beta[S_id[n]],-(expect[n]-(1-expect[n]))*delta[S_id[n]]);
    }
  }

}

generated quantities {

  real mu_threshold = gm[1];
  real mu_slope = gm[2];
  real mu_lapse = gm[3];
  real mu_alpha = gm[4];
  real mu_beta = gm[5];
  real mu_delta = gm[6];

  real tau_threshold = tau_u[1];
  real tau_slope = tau_u[2];
  real tau_lapse = tau_u[3];
  real tau_alpha = tau_u[4];
  real tau_beta = tau_u[5];
  real tau_delta = tau_u[6];



  vector[N] log_lik = rep_vector(0,N);



  for(n in 1:N){
    int c = binom_y[n];
    if(c == 1){
      log_lik[n] = wiener_lpdf(RT[n] | alpha[S_id[n]],rt_ndt[S_id[n]],beta[S_id[n]],(expect[n]-(1-expect[n]))*delta[S_id[n]]);
    }else if(c == 0){
      log_lik[n] = wiener_lpdf(RT[n] | alpha[S_id[n]],rt_ndt[S_id[n]],1-beta[S_id[n]],-(expect[n]-(1-expect[n]))*delta[S_id[n]]);
    }
  }


}
