library(invgamma)
library(mvtnorm)

library(rstan)
library(durhamSLR)
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

library(scales) #colours

statespace <- function(phi){
  p=length(phi)
  g=diag(x=1, nrow = p, ncol = p)[-p,]
  G=rbind(phi, g)
  return(G)
}

truephi=c(0.8,0.4,-0.3) ; truesigma2=4
abs(polyroot(c(1,0.8,0.4,-0.3)))>1
set.seed(8)
y <- arima.sim(n=1000, list(ar = truephi), sd = 2)

abs(eigen(statespace(truephi))$values) < 1
plot(ts(y), col=4 , main = 'Simulated AR(3) process')


NormalIGGibbs <- function(N, data, sigma2.p, phi.p){ #doesn't actually use data structure beyond this prior so it's fine
  p <- 3 #so i stay sane
  T <- length(data) - p #must specify p in order to get the prior - will extend later
  mat <- matrix(rep(0,N*(p+1)), ncol = (p+1), nrow = N) #output for diagnostics - row 1 sigmasq, 
  #Prior specification
  a <- spri[1] #prior on a
  b <- spri[2] #prior on b
  mat[1,1] <- rinvgamma(1,a,b) #sample from prior
  mat[1,c(2:(p+1))] <- rmvnorm(1, mean = rep(0,p), sigma = ppri) #we have a problem of low eigenvalues here
  Vinv <- solve(ppri) #assume phi normal(0, V) prior so for simplicity i will take 0 mean
  # response variable and design matrix
  yT <- data[(p+1):T]
  X <- matrix(rep(0, p*(T-p)), nrow = T-p, ncol = p) #i should be pulling from data NOT yT
  for(j in 1:(T-p)){
    xj <- data[j:(j+p-1)]
    X[j,] <- rev(xj)
  }
  
  #gibbs sampler
  for(i in 2:N){
    sigma2old <- mat[i-1,1]
    phiold <- mat[i-1,c(2:(p+1))]
    
    beta <- 0.5*t(yT - X%*%phiold)%*%(yT - X%*%phiold) + b #updated quantities
    SigmaMat <- solve(Vinv + 1/(sigma2old)*t(X)%*%X)
    Mu <- 1/(sigma2old)*SigmaMat%*%t(X)%*%yT
    
    mat[i,1] <- rinvgamma(1, (T-p)/2, beta) #new sigma draw
    mat[i,c(2:(p+1))] <- rmvnorm(1, mean = Mu, sigma = SigmaMat) #new phi draw
    
  }
  return(mat)
}
#priors

spri <- c(2,2)
ppri <- diag(0.1,3)
outputgibbs <- NormalIGGibbs(10000, y, spri, ppri)

plot(ts(outputgibbs[-c(1:1000), 1]), main = "Trace plot", ylab='sigma^2')
plot(ts(outputgibbs[-c(1:1000), 2]), main = "Trace plot", ylab = 'phi_1')
plot(ts(outputgibbs[-c(1:1000), 3]), main = "Trace plot", ylab = 'phi_2')
plot(ts(outputgibbs[-c(1:1000), 4]), main = "Trace plot", ylab = 'phi_3')
plot(acf(outputgibbs[-c(1:1000),1]), main = "sigma^2 ACF")
plot(acf(outputgibbs[-c(1:1000),2]), main = "phi_1 ACF")
plot(acf(outputgibbs[-c(1:1000),3]), main = "phi_2 ACF")
plot(acf(outputgibbs[-c(1:1000),4]), main = "phi_3 ACF")

#variance
#colours - 4,10,cyan


###################################
#  Posterior summary for Gibbs    #
###################################


hist(outputgibbs[-c(1:1000),1], breaks=20, probability = TRUE,
     main='Posterior Summary', xlab=NULL)
lines(density(outputgibbs[-c(1:1000),1]), col = 4, lwd=3)
abline(v=truesigma2,col='red',lwd=2)

hist(outputgibbs[-c(1:1000),2], breaks=20, probability = TRUE,
     main='Posterior Summary', xlab=NULL)
lines(density(outputgibbs[-c(1:1000),2]), col = 4, lwd=3)
abline(v=truephi[1],col='red',lwd=2)

hist(outputgibbs[-c(1:1000),3], breaks=20, probability = TRUE,
     main='Posterior Summary', xlab=NULL)
lines(density(outputgibbs[-c(1:1000),3]), col = 4, lwd=3)
abline(v=truephi[2],col='red',lwd=2)

hist(outputgibbs[-c(1:1000),4], breaks=20, probability = TRUE,
     main='Posterior Summary', xlab=NULL)
lines(density(outputgibbs[-c(1:1000),4]), col = 4, lwd=3)
abline(v=truephi[3],col='red',lwd=2)


#####################################
# Convergence diagnostics for Gibbs #
#####################################

#######################
# stan for comparison #
#######################

data <- list(T = 1000,
             p = 3, 
             y = y,
             a = 2,
             b = 3,
             s = 0.1)
output <- stan('nonstationary.stan', iter = 10000,
               data = data, chains = 4) #this now works
Usefuloutput <- extract(output, c("sigma2", "phi"))


#####################################
# Graphs for comparison stan/gibbs  #
#####################################

hist(outputgibbs[-c(1:1000),1], breaks=20, probability = TRUE,
     main='Posterior Summary', xlab=NULL, col = alpha('purple', 0.4))
hist(Usefuloutput$sigma2, breaks=20, probability = TRUE,
     xlab=NULL, add=TRUE, col = alpha('orange',0.5))
lines(density(outputgibbs[-c(1:1000),1]), col = 'purple', lwd=3)
lines(density(Usefuloutput$sigma2), col = 'red', lwd=3)
abline(v=truesigma2,col='cyan',lwd=3)
legend('topright', legend = c('Gibbs density estimate','NUTS density estimate','True value'),
       fill = c('purple', 10, 'cyan'))

hist(outputgibbs[-c(1:1000),2], breaks=20, probability = TRUE,
     main='Posterior Summary', xlab=NULL, col = alpha('purple', 0.4))
hist(Usefuloutput$phi[,1], breaks=20, probability = TRUE,
     xlab=NULL, add=TRUE, col = alpha('orange',0.5))
lines(density(outputgibbs[-c(1:1000),2]), col = 'purple', lwd=3)
lines(density(Usefuloutput$phi[,1]), col = 'red', lwd=3)
abline(v=truephi[1],col='cyan',lwd=3)
legend('topright', legend = c('Gibbs density estimate','NUTS density estimate','True value'),
       fill = c('purple', 10, 'cyan'))

hist(outputgibbs[-c(1:1000),3], breaks=20, probability = TRUE,
     main='Posterior Summary', xlab=NULL, col = alpha('purple', 0.4))
hist(Usefuloutput$phi[,2], breaks=20, probability = TRUE,
     xlab=NULL, add=TRUE, col = alpha('orange',0.5))
lines(density(outputgibbs[-c(1:1000),3]), col = 'purple', lwd=3)
lines(density(Usefuloutput$phi[,2]), col = 'red', lwd=3)
abline(v=truephi[2],col='cyan',lwd=3)
legend('topright', legend = c('Gibbs density estimate','NUTS density estimate','True value'),
       fill = c('purple', 10, 'cyan'))

hist(outputgibbs[-c(1:1000),4], breaks=20, probability = TRUE,
     main='Posterior Summary', xlab=NULL, col = alpha('purple', 0.4))
hist(Usefuloutput$phi[,3], breaks=20, probability = TRUE,
     xlab=NULL, add=TRUE, col = alpha('orange',0.5))
lines(density(outputgibbs[-c(1:1000),4]), col = 'purple', lwd=3)
lines(density(Usefuloutput$phi[,3]), col = 'red', lwd=3)
abline(v=truephi[3],col='cyan',lwd=3)
legend('topright', legend = c('Gibbs density estimate','NUTS density estimate','True value'),
       fill = c('purple', 10, 'cyan'))

