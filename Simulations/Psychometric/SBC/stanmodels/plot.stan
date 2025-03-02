functions {
  
  real psychometric(real x,real threshold, real slope,real lapse){
  
    return(lapse+(1-2*lapse)*(1 / (1+exp(-slope * (x - threshold)))));
      
  }
    vector psychometric_vec(vector x,real threshold, real slope,real lapse){
  
    return(lapse+(1-2*lapse).*(1 / (1+exp(-slope .* (x - threshold)))));
      
  }
  
  
  real entropy(real p){
    return(-p * log(p) - (1-p) * log(1-p));
  }
  
  vector entropy_vec(vector p){
    return(-p .* log(p) - (1-p) .* log(1-p));
  }
  
  
  real response_times(real unc, real rt_int, real slope){
    return(rt_int + slope * unc);
  }
      
  real gauss_copula_cholesky_lpdf(matrix u, matrix L) {
    array[rows(u)] row_vector[cols(u)] q;
    for (n in 1:rows(u)) {
      q[n] = inv_Phi(u[n]);
    }

    return multi_normal_cholesky_lpdf(q | rep_row_vector(0, cols(L)), L)
            - std_normal_lpdf(to_vector(to_matrix(q)));
  }

  vector gauss_copula_cholesky_pointwise(matrix u, matrix L) {
    int N = rows(u);
    int J = cols(u);
    matrix[J,J] Sigma_inv = chol2inv(L);
    vector[J] inv_sigma_inv = inv(diagonal(Sigma_inv));
    matrix[N, J] log_lik_mat;
    matrix[N, J] G;
    matrix[N, J] q;

    for (n in 1:N) {
      q[n] = inv_Phi(u[n]);
    }

    G = q * Sigma_inv;

    for (n in 1:N) {
      for (j in 1:J) {
        log_lik_mat[n, j] = normal_lpdf(q[n, j] | q[n, j] - G[n, j] * inv_sigma_inv[j], sqrt(inv_sigma_inv[j]))
                              - std_normal_lpdf(q[n, j]);
      }
    }
    return to_vector(log_lik_mat);
  }
  
  matrix gauss_copula_cholesky_pointwise_full(matrix u, matrix L) {
    int N = rows(u);
    int J = cols(u);
    matrix[J,J] Sigma_inv = chol2inv(L);
    vector[J] inv_sigma_inv = inv(diagonal(Sigma_inv));
    matrix[N, J] log_lik_mat;
    matrix[N, J] G;
    matrix[N, J] q;

    for (n in 1:N) {
      q[n] = inv_Phi(u[n]);
    }

    G = q * Sigma_inv;

    for (n in 1:N) {
      for (j in 1:J) {
        log_lik_mat[n, j] = normal_lpdf(q[n, j] | q[n, j] - G[n, j] * inv_sigma_inv[j], sqrt(inv_sigma_inv[j]))
                              - std_normal_lpdf(q[n, j]);
      }
    }
    return (log_lik_mat);
  }

  matrix uvar_bounds(array[] int binom_y, vector gm,vector tau_u, matrix z_expo,array[] int S_id, int S, vector x, array[] int starts, array[] int ends, int is_upper) {
    int N = size(binom_y);
                       
    matrix[N, 1] u_bounds;
    
    
    matrix[3, S] param;
  
    // adding the participant deviation to the group means
    for(p in 1:3){
      param[p,] = to_row_vector(gm[p] + z_expo[p,] * tau_u[p]);
    }
  
  
  row_vector[S] threshold = (param[1,]);
  row_vector[S] slope = exp(param[2,]);
  row_vector[S] lapse = inv_logit(param[3,]) / 2;
    
    
    vector[size(x)] expect;
    
    
    for (s in 1:S) {
      expect[starts[s]:ends[s]] = psychometric_vec(x[starts[s]:ends[s]],threshold[s],(slope[s]),(lapse[s])); 
    }

    for (n in 1:N) {
      if (is_upper == 0) {
        u_bounds[n, 1] = binom_y[n] == 0.0
                          ? 0.0 : binomial_cdf(binom_y[n] - 1 | 1, expect[n]);
      } else {
        u_bounds[n, 1] = binomial_cdf(binom_y[n] | 1, expect[n]);
      }
    }

    return u_bounds;
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

  array[S] cholesky_factor_corr[2] rho_chol;
  
  vector<lower=0, upper = minRT>[S] rt_ndt;

  matrix<
    lower=uvar_bounds(binom_y, gm, tau_u,z_expo,S_id,S, x, starts, ends, 0),
    upper=uvar_bounds(binom_y, gm, tau_u,z_expo,S_id,S, x, starts, ends, 1)
  >[N, 1] u;
  
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
  
  row_vector[S] rt_int = param[4,];
  row_vector[S] rt_beta = param[5,];
  row_vector[S] rt_sd = exp(param[6,]);
  
  
  
  vector[size(x)] expect;
    
  // for (n in 1:N) {
  //    expect[n] = psychometric(x[n],threshold[S_id[n]],(slope[S_id[n]]),(lapse[S_id[n]])); 
  // }
  
  for (s in 1:S) {
    expect[starts[s]:ends[s]] = psychometric_vec(x[starts[s]:ends[s]],threshold[s],(slope[s]),(lapse[s])); 
  }
  matrix[N, 2] u_mix;

  for (n in 1:N) {
    u_mix[n, 1] = u[n,1];
    u_mix[n, 2] = lognormal_cdf(RT[n] - rt_ndt[S_id[n]] | rt_int[S_id[n]] + rt_beta[S_id[n]] * entropy(expect[n]), (rt_sd[S_id[n]]));
  }
  
  
}


model {
  
  gm[1] ~ normal(0,5);
  gm[2] ~ normal(-1.5,0.5);
  gm[3] ~ normal(-3,1);
  gm[4] ~ normal(-1,0.25);
  gm[5] ~ normal(1.5,0.5);
  gm[6] ~ normal(-1,0.25);  

  to_vector(z_expo) ~ std_normal();

  tau_u[1] ~ normal(5,2);
  tau_u[2] ~ normal(0.6,0.2);
  tau_u[3] ~ normal(0.6,0.2);
  tau_u[4] ~ normal(0.6,0.2);
  tau_u[5] ~ normal(0.6,0.2);
  tau_u[6] ~ normal(0.5,0.1);  

  rt_ndt ~ normal(0.3,0.05);
  
  // for (n in 1:N) {
  //   target += lognormal_lpdf(RT[n] - rt_ndt[S_id[n]]| rt_int[S_id[n]] + rt_beta[S_id[n]] * entropy(expect[n]), rt_sd[S_id[n]]);
  // }

  for(s in 1:S){
    target += lognormal_lpdf(RT[starts[s]:ends[s]] - rt_ndt[s]| rt_int[s] + rt_beta[s] * entropy_vec(expect[starts[s]:ends[s]]), rt_sd[s]);
    rho_chol[s] ~ lkj_corr_cholesky(12);
    u_mix[starts[s]:ends[s],] ~ gauss_copula_cholesky(rho_chol[s]);
  }

}

generated quantities {
  // 
  vector[S] rho;
  
  for(s in 1:S){
    rho[s] = multiply_lower_tri_self_transpose(rho_chol[s])[1, 2];
  }
  real mu_threshold = gm[1];
  real mu_slope = gm[2];
  real mu_lapse = gm[3];
  real mu_rt_int = gm[4];
  real mu_rt_beta = gm[5];
  real mu_rt_sd = gm[6];
 
  real tau_threshold = tau_u[1];
  real tau_slope = tau_u[2];
  real tau_lapse = tau_u[3];
  real tau_rt_int = tau_u[4];
  real tau_rt_beta = tau_u[5];
  real tau_rt_sd = tau_u[6];
  
  vector[N] expect_gq = expect;
  vector[N] mu_rts_gq;
  vector[N] ndt_rts_gq;
  vector[N] sd_rts_gq;
  
  
  for(n in 1:N){
    mu_rts_gq[n] = rt_int[S_id[n]] + rt_beta[S_id[n]] * entropy(expect[n]);
    sd_rts_gq[n] = rt_sd[S_id[n]];
    ndt_rts_gq[n] = rt_ndt[S_id[n]];
  }


}
