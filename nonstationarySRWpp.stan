//random walk in R on the autoregressive coefficients, not in the stationary region
data {
  int<lower=0> p;       //order
  int<lower=0> T;       //Time
  vector[T] Y;          //Data
  real<lower=0> a;      //prior
  real<lower=0> b;      //prior
  real<lower=0> alpha;  //prior
  real<lower=0> beta;   //prior
  real<lower=0> s;      //initialiser for SRW
  int<lower=0> Tstar; // steps ahead
}

parameters {
  real<lower=0> sigmaR;      //RW variance
  real<lower=0> sigmaA;      //ARp variance
  matrix[T,p] phi; //random walk
}

model { 
  real mu;
  //prior
  sigmaA ~ inv_gamma(a,b);
  sigmaR ~ inv_gamma(alpha,beta);
  for(i in 1:p) {
    phi[1, i] ~ normal(0,s);
  }
  for(t in 2:T){
      phi[t,] ~ normal(phi[(t-1),],sqrt(sigmaR)); //RW
  }
  for(t in (p+1):T){
    mu=0;
    for(i in 1:p){
      mu += phi[t, i] * Y[t-i];
    }
    Y[t] ~ normal(mu, sqrt(sigmaA));
  }
}

generated quantities { //posterior predictive densities for next T* values
  vector[p + Tstar] Ystar; //final p values and Tstar predictions so I'm not there ALL day
  matrix[p + Tstar,p] phi_star; //the next steps in the random walk
  real mu; //construction of AR(p)
  Ystar[1:p] = Y[(T-p+1):T]; //last p values
  for(i in 1:p){
    for(j in 1:p){
      phi_star[i,j] = phi[(T-p+i),j];
    }
  }
  for(i in 1:p){ //first new random walk using estimate for sigmaR
    phi_star[p + 1,i] = normal_rng(phi[T,i], sqrt(sigmaR));
  }
  
  for(t in 2:Tstar) {
    for(j in 1:p) {
     phi_star[(p+t),j] = normal_rng(phi_star[(p+t-1),j], sqrt(sigmaR)); 
    }
  }
  for(t in 1:Tstar){ 
  //handles the first case also, i.e. the block generates the next T* values in the RW
  //then it computes the AR(p) predictions.
    mu=0;
    for(k in 1:p) {
      mu += phi_star[(p+t), k] * Ystar[p + t - k];
    }
    Ystar[p + t] = normal_rng(mu, sigmaA);
  }
}
