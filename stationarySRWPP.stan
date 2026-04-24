functions {
  vector r2phi(vector r) {
    //Transformation from Marriott et al.
    //Make matrix put PAC on diagonal
    //Recursion on each PAC until end
    //Report final column as autocorrelation coefficients
    int p = size(r);
    matrix[p,p] phimat;
    phimat[1,1] = r[1];
    for(k in 2:p){
      phimat[k,k] = r[k];
      for(i in 1:(k-1)){
        phimat[i,k] = phimat[i,(k-1)] - r[k] * phimat[(k-i),(k-1)];
      }
    }
    vector[p] phi = col(phimat, p);
    return phi;
  }
}

data {
  int<lower=0> p;       //order
  int<lower=0> T;       //Time
  vector[T] Y;          //Data
  real<lower=0> a;      //prior
  real<lower=0> b;      //prior
  real<lower=0> alpha;  //prior
  real<lower=0> beta;   //prior
  real<lower=0> s;      //initialiser for SRW
  real<lower=0> kappap;  //kappa prior
  int<lower=0> Tstar; // steps ahead
}

parameters {
  real<lower=0> sigmaR;      //RW variance
  real<lower=0> sigmaA;      //ARp variance
  matrix[T,p] rst; //random walk
  real kappa;
}

transformed parameters {
  matrix[T,p] phi;
  for(t in 1:T){
    phi[t,] = r2phi( ((exp(rst[t,]) - 1 ) ./ (exp(rst[t,]) + 1 ))' )';
  }
}

model { 
  real mu;
  //prior
  kappa ~ normal(0,kappap);
  sigmaA ~ inv_gamma(a,b);
  sigmaR ~ inv_gamma(alpha,beta);
  for(i in 1:p) {
    rst[1, i] ~ normal(0,s);
  }
  for(t in 2:T){
      rst[t,] ~ normal(kappa*rst[(t-1),],sqrt(sigmaR)); //AR1
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
  matrix[p + Tstar, p] rstar; //the next steps in the random walk
  matrix[p + Tstar,p] phi_star; //the AR coefficients computed as a result of the RW
  real mu; //construction of AR(p)
  Ystar[1:p] = Y[(T-p+1):T]; //last p values so that
  for(i in 1:p){
    for(j in 1:p){
      phi_star[i,j] = phi[(T-p+i),j];
    }
  }
  
  //Put all the first p values into the Ystar, rstar, phi_star
  for(t in 1:p){
    rstar[t,] = rst[(T-p+t), ]; //final p values of the RW on the og chain
  }
  for(i in 1:p){ //first new random walk using estimate for sigmaR
    rstar[p + 1,i] = normal_rng(rst[T,i], sqrt(sigmaR));
  }
  
  for(t in 2:Tstar) {
    for(j in 1:p) {
     rstar[(p+t),j] = normal_rng(rstar[(p+t-1),j], sqrt(sigmaR)); 
    }
  }
  for(t in 1:Tstar) { //generate phi
    phi_star[p+t,] = r2phi( ((exp(rstar[t,]) - 1 ) ./ (exp(rstar[t,]) + 1 ))' )';
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
