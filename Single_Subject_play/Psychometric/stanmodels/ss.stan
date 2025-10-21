functions {
  
  
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
  
  real psycho(real x, real alpha, real beta, real lapse){
   return (lapse + (1 - 2 * lapse) * (tanh(beta*(x-alpha)) / 2 + 0.5));
  }

  real entropy(real p){
    return(-p * log(p) - (1-p) * log(1-p));
  }

  real contin_resp(real unc, real rt_int, real slope){
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

  matrix uvar_bounds(array[] int binom_y, vector gm, vector X,
                     int is_upper) {
    int N = size(binom_y);
                       
    matrix[N, 1] u_bounds;
    

    
  real alpha = (gm[1]);
  real beta = exp(gm[2]);
  real lapse = inv_logit(gm[3]) / 2;

  
    for (n in 1:N) {
      real theta = psycho(X[n], alpha, beta, lapse);
      if (is_upper == 0) {
        u_bounds[n, 1] = binom_y[n] == 0.0
                          ? 0.0 : binomial_cdf(binom_y[n] - 1 | 1, theta);
      } else {
        u_bounds[n, 1] = binomial_cdf(binom_y[n] | 1, theta);
      }
    }

    return u_bounds;
  }
  
  
  
  real ord_beta_reg_cdf(real y, real mu, real phi, real cutzero, real cutone) {

    vector[2] thresh;
    thresh[1] = cutzero;
    thresh[2] = cutzero + exp(cutone);

    real p0 = 1-inv_logit(mu - thresh[1]);

    real p_m = (inv_logit(mu - thresh[1])-inv_logit(mu - thresh[2]))  * beta_cdf(y | exp(log_inv_logit(mu) + log(phi)), exp(log1m_inv_logit(mu) + log(phi)));



    if (y < 0) {
      return 0;
    } else if (y == 0) {
      return p0;
    } else if (y == 1) {
      return 1-(1e-12);
    } else {
      return (p0 + p_m);
    }
  }

  // ordered beta function
  real ord_beta_reg_lpdf(real y, real mu, real phi, real cutzero, real cutone) {

    vector[2] thresh;
    thresh[1] = cutzero;
    thresh[2] = cutzero + exp(cutone);

  if(y==0) {
      return log1m_inv_logit(mu - thresh[1]);
    } else if(y==1) {
      return log_inv_logit(mu  - thresh[2]);
    } else {
      return log_diff_exp(log_inv_logit(mu - thresh[1]), log_inv_logit(mu - thresh[2])) +
                beta_lpdf(y|exp(log_inv_logit(mu) + log(phi)),exp(log1m_inv_logit(mu) + log(phi)));
    }
  }

  real induced_dirichlet_lpdf(real nocut, vector alpha, real phi, int cutnum, real cut1, real cut2) {
    int K = num_elements(alpha);
    vector[K-1] c = [cut1, cut1 + exp(cut2)]';
    vector[K - 1] sigma = inv_logit(phi - c);
    vector[K] p;
    matrix[K, K] J = rep_matrix(0, K, K);

    if(cutnum==1) {

    // Induced ordinal probabilities
    p[1] = 1 - sigma[1];
    for (k in 2:(K - 1))
      p[k] = sigma[k - 1] - sigma[k];
    p[K] = sigma[K - 1];

    // Baseline column of Jacobian
    for (k in 1:K) J[k, 1] = 1;

    // Diagonal entries of Jacobian
    for (k in 2:K) {
      real rho = sigma[k - 1] * (1 - sigma[k - 1]);
      J[k, k] = - rho;
      J[k - 1, k] = rho;
    }

    // divide in half for the two cutpoints

    // don't forget the ordered transformation

      return   dirichlet_lpdf(p | alpha)
           + log_determinant(J) + cut2;
    } else {
      return(0);
    }
  }
}
  
  

data {
  int<lower=0> N;
  
  array[N] int binom_y;
  vector[N] RT;
  vector[N] Conf;
  
  vector[N] X;
  
  real minRT;
  
  vector[N] ACC; // Vector of deltaBPM values that match the binary response



}
transformed data{
  int P = 12;
}

parameters {
  vector[P] gm;

  matrix<
    lower=uvar_bounds(binom_y, gm, X, 0),
    upper=uvar_bounds(binom_y, gm, X, 1)
  >[N, 1] u;
  
  // cholesky_factor_corr[2] rho_chol;
  
  cholesky_factor_corr[3] rho_chol;
  
  real c0;
  real c11;
  real <lower=0, upper = minRT> rt_ndt;

}

transformed parameters{
  

  
  real alpha = (gm[1]);
  real beta = exp(gm[2]);
  real lapse = inv_logit(gm[3]) / 2;
  
  real rt_int = gm[4];
  real rt_slope = gm[5];
  real rt_prec = exp(gm[6]);
  real rt_stim = gm[7];


  real conf_int = gm[8];
  real conf_ACC = gm[9];
  real conf_entropy = gm[10];
  real conf_entropy_ACC = gm[11];
  real conf_prec = exp(gm[12]);

  
  
  vector[N] entropy_t;

  vector[N] conf_mu;
  
  profile("likelihood") {

  for (n in 1:N){
    
  entropy_t[n] = entropy(psycho(X[n], alpha, beta, lapse));
  
  conf_mu[n] = conf_int +                                           // intercept
    conf_ACC * ACC[n] +                                  // main effect: ACC
    conf_entropy * entropy_t[n] +                        // main effect: entropy
  
    conf_entropy_ACC * ACC[n] * entropy_t[n];                 // 2-way interaction: ACC × entropy
  }
  
  }
  
  
  matrix[N, 3] u_mix;
  for (n in 1:N) {
    
    u_mix[n, 1] = u[n,1];
    
    u_mix[n, 2] = lognormal_cdf(RT[n] - rt_ndt | rt_int + rt_slope * entropy_t[n] + rt_stim * X[n], rt_prec);
    
    u_mix[n, 3] = ord_beta_reg_cdf(Conf[n] | conf_mu[n], conf_prec, c0, c11);
    
  }
  
  
  
}

model {

  gm[1] ~ normal(0,10); //global mean of beta
  gm[2] ~ normal(-2,3); //global mean of beta
  gm[3] ~ normal(-4,2); //global mean of beta
  gm[4:12] ~ normal(0,5); //global mean of beta
  

  rt_ndt ~ normal(0.3,0.05);

  
  for (n in 1:N) {
    
    target += lognormal_lpdf(RT[n] - rt_ndt | rt_int + rt_slope * entropy_t[n]+ rt_stim * X[n], rt_prec);
    
    target += ord_beta_reg_lpdf(Conf[n] | conf_mu[n], conf_prec, c0, c11);
    
    
  }
    
    c0 ~ induced_dirichlet([1,10,1]', 0, 1, c0, c11);
    c11 ~ induced_dirichlet([1,10,1]', 0, 2, c0, c11);

    rho_chol ~ lkj_corr_cholesky(12);

    u_mix ~ gauss_copula_cholesky(rho_chol);
  
}

generated quantities {
  // 
  real rho_p_rt;
  real rho_p_conf;
  real rho_rt_conf;
  


  rho_p_rt = multiply_lower_tri_self_transpose(rho_chol)[1, 2];
  rho_p_conf = multiply_lower_tri_self_transpose(rho_chol)[1, 3];
  rho_rt_conf = multiply_lower_tri_self_transpose(rho_chol)[2, 3];

  matrix[N,3] log_lik_cop_full;

  vector[N] log_lik_bin = rep_vector(0,N);
  vector[N] log_lik_rt  = rep_vector(0,N);
  vector[N] log_lik_conf  = rep_vector(0,N);
  vector[N] log_lik = rep_vector(0,N);


 log_lik_cop_full = gauss_copula_cholesky_pointwise_full(u_mix,rho_chol);

  for(n in 1:N){

    log_lik_rt[n] = lognormal_lpdf(RT[n] - rt_ndt | rt_int + rt_slope * entropy_t[n]+ rt_stim * X[n], rt_prec);

    log_lik_bin[n] = binomial_lpmf(binom_y[n] | 1, psycho(X[n], alpha, beta, lapse));

    log_lik_conf[n] = ord_beta_reg_lpdf(Conf[n] | conf_mu[n], conf_prec, c0, c11);

    log_lik[n] = log_lik_rt[n] + log_lik_bin[n] + log_lik_conf[n] + log_lik_cop_full[n,1] + log_lik_cop_full[n,2]+ log_lik_cop_full[n,3];
  }
  




}
