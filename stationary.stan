functions {
  matrix kprod(matrix A, matrix B){
	int n_A = rows(A);
	int n_B = rows(B);
	int n_Z = n_A * n_B;
	matrix[n_Z, n_Z] Z;
	for (i in 1:n_A) {
		for (j in 1:n_A) {
			for (p in 1:n_B) {
				for (q in 1:n_B) {
					int row = (i - 1) * n_B + p;
					int col = (j - 1) * n_B + q;
					Z[row, col] = A[i, j] * B[p, q];
				}
			}
		}
	}
	return(Z);
  }
  
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
  
  matrix stationaryvar(vector phi, real Var) {
    //state space matrix 
    //Kronecker product of state space with itself
    //VecV = VecSigma * inv(Ip2 - GXG)
    //V = matrix(VecV)
    int p = num_elements(phi);
    matrix[p,p] G = rep_matrix(0,p,p);
    G[1,] = (phi)';
    for(i in 1:(p-1)){
      G[(i+1),i] = 1;
      }
    /*matrix[p*p, p*p] Mat = rep_matrix(0,p*p,p*p); //p blocks
    
    for(i in 1:p){ //i think this construction does work 
      //row assignment
      for(j in 1:p) {
        Mat[i,((j-1)*p+1):(j*p)] = G[i,j] * G[i];
      }
    }*/
    matrix[p*p,p*p] Mat = -1 * add_diag(kprod(G,G),-1);
    
    vector[p*p] Sig = rep_vector(0, p*p);
    Sig[1] = Var;
    
    row_vector[p*p] vecV = Sig' / Mat';
    matrix[p,p] V;
    for(i in 1:p) {
      V[i] = segment(vecV, (i-1)*p+1, p);
      //vecV[((i-1)*p+1):(i*p)];
    } 
    return V;
  }
}

data {
  int<lower=0> T;     //size of data
  int<lower=0> p;     //order of AR
  vector[T] y;      //Time series in total; no need to generate stationary var here
  real<lower=0> a;    //sig prior
  real<lower=0> b;    //sig prior
  vector<lower=0>[p] alpha; //prior on betas
  vector<lower=0>[p] beta;  //prior on betas
}

parameters {
  real<lower=0> sigma2;
  vector<lower=-1,upper=1>[p] r;
}

transformed parameters {
  vector<lower=0,upper=1>[p] rst = 0.5 + 0.5*r;
  vector[p] phi = r2phi(r);
  //matrix[p,p] V = stationaryvar(phi,sigma2);
  cov_matrix[p] V = stationaryvar(phi,sigma2);
}

model {
  sigma2 ~ inv_gamma(a, b); //prior
  for(n in 1:p){
    rst[n] ~ beta(alpha[n],beta[n]);
  }
  y[1:p] ~ multi_normal(rep_vector(0,p),V); //likelihood
  for(t in (p+1):T) {
    real mu = 0;
    for(i in 1:p){
      mu += phi[i] * y[t-i];
    }
    y[t] ~ normal(mu, sqrt(sigma2));
  }
}
