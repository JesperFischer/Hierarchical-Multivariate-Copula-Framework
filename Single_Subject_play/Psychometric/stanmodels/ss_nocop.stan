functions {
  
  
  real psycho(real x, real alpha, real beta, real lapse){
   return (lapse + (1 - 2 * lapse) * (tanh(beta*(x-alpha)) / 2 + 0.5));
  }

  real entropy(real p){
    return(-p * log(p) - (1-p) * log(1-p));
  }

  real contin_resp(real unc, real rt_int, real slope){
    return(rt_int + slope * unc);
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
  
  vector[N] prob;

  
  profile("likelihood") {

  for (n in 1:N){
    
    
    prob[n] = psycho(X[n], alpha, beta, lapse);
    
  entropy_t[n] = entropy(psycho(X[n], alpha, beta, lapse));
  
  conf_mu[n] = conf_int +                                           // intercept
    conf_ACC * ACC[n] +                                  // main effect: ACC
    conf_entropy * entropy_t[n] +                        // main effect: entropy
  
    conf_entropy_ACC * ACC[n] * entropy_t[n];                 // 2-way interaction: ACC × entropy
  }
  
  }
}

model {

  gm[1] ~ normal(0,10); //global mean of beta
  gm[2] ~ normal(-2,3); //global mean of beta
  gm[3] ~ normal(-4,2); //global mean of beta
  gm[4:12] ~ normal(0,5); //global mean of beta
  

  rt_ndt ~ normal(0.3,0.05);

  

  matrix[N, 3] u_mix;
  for (n in 1:N) {

    target += lognormal_lpdf(RT[n] - rt_ndt | rt_int + rt_slope * entropy_t[n]+ rt_stim * X[n], rt_prec);
    
    target += ord_beta_reg_lpdf(Conf[n] | conf_mu[n], conf_prec, c0, c11);
    
    
  }
    
    c0 ~ induced_dirichlet([1,10,1]', 0, 1, c0, c11);
    c11 ~ induced_dirichlet([1,10,1]', 0, 2, c0, c11);


  
}


generated quantities{
  
  vector[N] log_lik_p;
  vector[N] log_lik_rt;
  vector[N] log_lik_c;
  vector[N] log_lik;
  
  for(n in 1:N){
    
    log_lik_p[n] = bernoulli_lpmf(binom_y[n] | prob[n]);
    log_lik_rt[n] = lognormal_lpdf(RT[n] - rt_ndt | rt_int + rt_slope * entropy_t[n]+ rt_stim * X[n], rt_prec);
    log_lik_c[n] = ord_beta_reg_lpdf(Conf[n] | conf_mu[n], conf_prec, c0, c11);
    
    log_lik[n] = log_lik_p[n] + log_lik_rt[n] + log_lik_c[n];
    
  }
}
