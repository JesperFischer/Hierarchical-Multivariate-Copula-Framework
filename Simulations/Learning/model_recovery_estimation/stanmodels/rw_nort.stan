functions {
 
  real learning(real x,real expect, real learningrate){
 
    return(expect + learningrate * (x - expect));
     
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
  array[N] int id_ind;
  array[S] int starts;
  array[S] int ends;
 





}
transformed data{
  int P = 2;
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
 
 
  row_vector[S] learningrate = inv_logit(param[1,]);
  row_vector[S] e0 = inv_logit(param[2,]);
 
 
  vector[size(x)+1] expect;


  for (n in 1:N) {
    if(id_ind[n] == 1){
      expect[n] = e0[S_id[n]];
    }else{
      expect[n] = learning(x[n-1],expect[n-1],learningrate[S_id[n]]);  
    }
  }
 
 
}


model {
 
  gm[1] ~ normal(-1,1);
  gm[2] ~ normal(0,0.2);

  to_vector(z_expo) ~ std_normal();

  tau_u[1] ~ normal(0.5,0.1);
  tau_u[2] ~ normal(0.5,0.1);

 

  for(s in 1:S){
    target += bernoulli_lpmf(binom_y[starts[s]:ends[s]] |  expect[starts[s]:ends[s]]);
  }


}

generated quantities {

  real mu_learningrate = gm[1];
  real mu_e0 = gm[2];

  real tau_learningrate = tau_u[1];
  real tau_e0 = tau_u[2];


  vector[N] log_lik_bin = rep_vector(0,N);
  vector[N] log_lik = rep_vector(0,N);
  


  for(n in 1:N){

    log_lik_bin[n] = binomial_lpmf(binom_y[n] | 1, expect[n]);

    log_lik[n] = log_lik_bin[n];
  }
 
}
