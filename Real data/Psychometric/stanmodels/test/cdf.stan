//
// This Stan program defines a simple model, with a
// vector of values 'y' modeled as normally distributed
// with mean 'mu' and standard deviation 'sigma'.
//
// Learn more about model development with Stan at:
//
//    http://mc-stan.org/users/interfaces/rstan.html
//    https://github.com/stan-dev/rstan/wiki/RStan-Getting-Started
//

functions {

  real ord_beta_reg_prob(real y, real mu, real phi, real cutzero, real cutone) {

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
      return 1;
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

// The input data is a vector 'y' of length 'N'.
data {
  int<lower=0> N; //number of trials
  vector[N] VAS;
}

// The parameters accepted by the model.
// subject level with normal distribution around the parameters
parameters {

  real mu;
  real phi;
  real c0;
  real c11;
}



// The model to be estimated.
model {


  mu ~ beta(5,10);
  phi ~ lognormal(3,1);
  c0 ~ induced_dirichlet([1,10,1]', 0, 1, c0, c11);
  c11 ~ induced_dirichlet([1,10,1]', 0, 2, c0, c11);

  // estimating phs and rt
  for(n in 1:N){
      target += ord_beta_reg_lpdf(VAS[n] | logit(mu), phi, c0, c11);

  }
}

generated quantities{

  real c1 = c0+exp(c11);

  vector[N] probs;
  vector[N] probs_beta;

  for (n in 1:N){
    probs[n] = ord_beta_reg_prob(VAS[n],logit(mu),phi,c0,c11);
    probs_beta[n] = beta_cdf(VAS[n]|mu*phi,(1-mu)*phi);

  }

}

