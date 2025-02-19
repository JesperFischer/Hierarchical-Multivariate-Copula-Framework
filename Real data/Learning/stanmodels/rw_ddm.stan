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
  array[N] int id_ind;
  array[N] int S_id; // Vector of binary responses
  array[S] int starts;
  array[S] int ends;
  





}
transformed data{
  int P = 5;
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
  
  
  row_vector[S] learningrate = inv_logit(param[1,]);
  row_vector[S] e0 = inv_logit(param[2,]);
  row_vector[S] alpha = exp(param[3,]);
  row_vector[S] beta = inv_logit(param[4,]);
  row_vector[S] delta = (param[5,]);
  
  
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
  gm[2] ~ normal(0,1);
  gm[3] ~ normal(0.5,0.5);
  gm[4] ~ normal(0,0.5);
  gm[5] ~ normal(2,1);
  
  to_vector(z_expo) ~ std_normal();
  
  tau_u[1] ~ normal(0.5,1);
  tau_u[2] ~ normal(0.5,1);
  tau_u[3] ~ normal(0.4,1);
  tau_u[4] ~ normal(0.4,1);
  tau_u[5] ~ normal(0.5,1);

  
  rt_ndt ~ normal(0.2,0.05);



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

  real mu_learningrate = gm[1];
  real mu_e0 = gm[2];
  real mu_alpha = gm[3];
  real mu_beta = gm[4];
  real mu_delta = gm[5];

  real tau_learningrate = tau_u[1];
  real tau_e0 = tau_u[2];
  real tau_alpha = tau_u[3];
  real tau_beta = tau_u[4];
  real tau_delta = tau_u[5];



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
