//Replication of the Normal/IG semi-conjugate case
data {
  int<lower=0> T;     //size of data
  int<lower=0> p;     //order of AR
  vector[T] y;        //Time series
  real<lower=0> a;    //sig prior
  real<lower=0> b;    //sig prior
  real<lower=0> s;      //prior variance of phi - in full generality should be a matrix, but restricting to iid phi here
}

parameters {
  real<lower=0> sigma2;
  vector[p] phi; //we do not confine to the stationary region here
}

model {
  sigma2 ~ inv_gamma(a, b);
  for(i in 1:p){
    phi[i] ~ normal(0, sqrt(s));
  }
  for(t in (p+1):T) {
    real mu=0;
    for(i in 1:p){
       mu += phi[i] * y[t-i];
    }
    y[t] ~ normal(mu, sqrt(sigma2));
  }
}
