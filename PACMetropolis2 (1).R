library(invgamma)
library(mvtnorm)
library(scales)

library(rstan)
#library(durhamSLR)
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)
################################
#PRIOR:

#sigma2 has an IG(a,b) prior
#Rho is a length p vector with rho[i] = 2*Y[i] - 1 
#Y[i] ~ Beta(alpha[i], beta[i])

#So would need to specify a, b, alpha and beta
#we will take our test model as an AR(4) process
p <- 4
#T <- 1000
a <- 2
b <- 2
alpha <- c(2,3,3,2)
beta <- c(3,2,2,3)

################################
#TOY MODEL:

#assumed stationary AR(4) process, generate synthetic data
truerho <- c(0.3,-0.2,-0.7,0.4)
truesigma2 <- 4
#this function will be used many times
pac.to.ar <- function(rho){ #type vector length p
  #make matrix
  p <- length(rho)
  mat <- diag(as.numeric(rho), nrow = p)
  k <- 2
  while(k<(p+1)){ #indexing :)
    for(i in 1:(k-1)){
      mat[i,k] <- mat[i,(k-1)]-rho[k]*mat[(k-i),(k-1)]
    }
    k <- k + 1
  }
  return(as.numeric(mat[,p]))
}
#y <- arima.sim(n=1000, list(ar = pac.to.ar(truerho)), sd = 2)
y <- arima.sim(n=1000, list(ar = pac.to.ar(truerho)), sd = sqrt(truesigma2)) # SEH: changed hard-coded std. dev
plot(y)
#################################
#MCMC SAMPLING
#We will sample from the posterior via MCMC,
#where we take a Gibbs sample of sigma2 from an IG distribution conditioned on rho,
#since we have semi-conjugacy
#We will take a Metropolis step on rho, where we sample on R, convert to [-1,1] and then compute acc.prob

#################################
#LIKELIHOOD:
#i will define the log-likelihood here; we require the stationary variance for this

statespace <- function(phi){
  p <- length(phi)
  g <- diag(x=1, nrow = p, ncol = p)[-p,]
  G <- rbind(phi, g)
  #rm(list = c('p', 'g')) # SEH: don't need this
  return(G)
}
stationaryvar <- function(phi, sigma2){ #phi coefficients here --equivalent to lyapunov
  #phi matrix constructed - called G here
  p <- length(phi)
  G <- statespace(phi)
  #Kronecker product
  sig <- matrix(data = c(sigma2, rep(0, p^2-1)), nrow = p)
  
  vec <- solve(diag(1,nrow = p^2, ncol = p^2)-(G%x%G))%*%as.numeric(sig)
  V <- matrix(vec, nrow = p, ncol = p)
  return(V)
}
#we also need design matrix, which is invariant
#X <- matrix(rep(0, p*(length(y)-p)), nrow = length(y)-p, ncol = p) #i should be pulling from data NOT yT
#for(j in 1:(length(y)-p)){
#  xj <- y[j:(j+p-1)]
#  X[j,] <- rev(xj)
#}
# put this in MCMC function

#model: f(yt | rho, sigma2) = f(y[1:p]|rho, sigma2) * f(y[(p+1)]|rho, sigma2)
#llike <- function(y, X, rho, sigma2) {
llike <- function(y, X, rho, sigma2) { # SEH: pass X as argument
  # SEH: extract p from args of function
  p <- length(rho)
  phi <- pac.to.ar(rho)
  #first p values
  V <- stationaryvar(phi, sigma2)
  val <- dmvnorm(y[1:p], rep(0,p), V, log = TRUE)
  #(p+1):T values
  #val <- val + dmvnorm(y[-c(1:p)], X%*%phi, diag(sigma2, T-p), log = TRUE)
  val <- val + dmvnorm(y[-c(1:p)], X%*%phi, diag(sigma2, length(y)-p), log = TRUE) # SEH: avoid using variables defined in global env
  return(val)
}
#llike(y, X, c(0.3,0.4,0.3,0.3), 2)

###############################
#GIBBS STEP
#sigma2Gibbs <- function(y, rho, sigma2, a, b) {
sigma2Gibbs <- function(y, X, rho, sigma2, a, b) { # SEH: as above
  p <- length(rho) # SEH: as above
  T <- length(y) # SEH: as above
  phi <- pac.to.ar(rho)
  omega <- 1/sigma2*stationaryvar(phi, sigma2)
  
  sigma2new <- rinvgamma(1,
                         #a + (T-p)/2 + 1,
                         a + (T-p)/2 - 0.5, # SEH: the contribution from the determinant of the stationary variance is (sigma^2)^{-p/2}
                         b + 0.5*(t(y[1:p])%*%solve(omega)%*%y[1:p] + t((y[(p+1):T]-X%*%phi))%*%(y[(p+1):T]-X%*%phi)))
  return(sigma2new)
}
#sigma2Gibbs(y, X, c(0.2,0.3,0.4,0.4), 0.4, 2,2) ##### here
###############################
#ACCEPTANCE PROBABILITY
#l.acc.prob = pi(rhocan) * q(rho|rhocan) / pi(rho) * q(rhocan|rho)
#pi is prior * likelihood
#proposal is a sample drawn from a random walk

#l.acc.prob <- function(y, rhocur, rhocan, sigma2, lambda) {
l.acc.prob <- function(y, X, rhocur, rhocan, sigma2, lambda) { # SEH: as above
  p <- length(rhocur) # SEH: as above
  #proposal density
  rhocur.st <- log(1+rhocur) - log(1-rhocur)
  rhocan.st <- log(1+rhocan) - log(1-rhocan)
  #density in rho.st space i.e. R, then the Jacobian takes that back down to [-1,1]
  #val <- dmvnorm(rhocur.st, mean = rhocan.st, sigma = diag(lambda^2, p), log = TRUE) -
  #  dmvnorm(rhocan.st, mean = rhocur.st, sigma = diag(lambda^2, p), log = TRUE)
  val <- 0 # SEH: the above is always zero by symmetry
  #Jacobian
  val <- val + p*log(2) - sum(log(1-rhocur) + log(1+rhocur)) - 
    p*log(2) + sum(log(1-rhocan) + log(1+rhocan))
  #likelihood
  val <- val + llike(y, X, rhocan, sigma2) - llike(y, X, rhocur, sigma2)
  #prior
  val <- val + sum(dbeta((rhocan + 1)/2, alpha, beta, log = TRUE)) -
    sum(dbeta((rhocur + 1)/2, alpha, beta, log = TRUE))
  return(min(0, val))
}
#l.acc.prob(y, X, c(0.1,0.2,0.3,0.4), c(0.2,0.3,0.4,0.5), 2, 0.5)
##############################
#Metropolis step
#rhoMetropolis <- function(y, rho, sigma2, lambda) {
rhoMetropolis <- function(y, X, rho, sigma2, lambda) { # SEH: as above
  p <- length(rho) # SEH: as above
  #Random generation
  #rhocan.st <- rmvnorm(1, rho, diag(lambda^2, p))
  rho.st <- log(1 + rho) - log(1 - rho) # SEH
  rhocan.st <- rmvnorm(1, rho.st, diag(lambda^2, p))[1,] # centre at current value in R
  rhocan <- (exp(rhocan.st) - 1 ) / (exp(rhocan.st) + 1)
  #acceptance
  A <- l.acc.prob(y, X, rho, rhocan, sigma2, lambda)
  u <- runif(1)
  val <- rho
  #if(log(u) > A) {
  if(log(u) < A) { # SEH: other way around!
    acc <- TRUE # SEH: keep track of acceptance count for tuning
    val <- rhocan
  } else acc <- FALSE
  #return(val)
  return(list(rho=val, acc=acc))
}

#rhoMetropolis(y, X, c(0.1,0.2,0.3,0.4), 2, 0.3)
#############################
#ALL TOGETHER

MCMCoutput <- function(N, y, prior, tuning) { #prior list (IG, betaalpha, betabeta)
  a <- prior[[1]][1]
  b <- prior[[1]][2]
  alpha <- prior[[2]]
  beta <- prior[[3]]
  p <- length(alpha) # SEH: as above
  # create "design matrix"
  X <- matrix(rep(0, p*(length(y)-p)), nrow = length(y)-p, ncol = p) #i should be pulling from data NOT yT
  for(j in 1:(length(y)-p)){
    xj <- y[j:(j+p-1)]
    X[j,] <- rev(xj)
  } # SEH: moved this bit here
  
  mat <- matrix(rep(0, N*(p+1)), nrow = N, ncol = p+1)
  #sample from prior
  for(j in 1:p){
    mat[1, j] <- rbeta(1, alpha[j], beta[j])
  }
  mat[1, p+1] <- rinvgamma(1, a, b)
  print('Start')
  #UPDATES
  acc <- 0
  for(i in 2:N) {
    #METROPOLIS
    tmp <- rhoMetropolis(y, X, mat[(i-1), 1:p], mat[(i-1),(p+1)], tuning)
    mat[i, 1:p] <- tmp$rho
    acc <- acc + tmp$acc # SEH: monitor acceptances
    #GIBBS
    mat[i, (p+1)] <- sigma2Gibbs(y, X, mat[i,1:p], mat[(i-1),(p+1)], a, b)
    #counter
    if(i %% floor(N/20) == 0){
      print(i)
    }
  }
  colnames(mat) = c(paste("phi[", 1:p, "]", sep=""), "sigma2") # SEH: added column names
  return(list(mcmc=mat, acc=acc/(N-1))) # SEH: return acceptance rate as well as MCMC
}

niters = 1e4
set.seed(66)
Output <- MCMCoutput(niters, y, list(c(a,b), alpha, beta), 0.08) 
Output$acc # SEH: check acceptance rate around 30%-ish
library(durhamSLR) # SEH: from BCM3
diagnostics(Output$mcmc, 3) # SEH: check mixing and convergence

#takes ages and sticks :( 
#SEH: still a bit slow but seems okay now! Remove first half
#of iterations as burn-in before producing these plots
par(mfrow=c(1,1), ask=TRUE) # SEH: so you can tab through
for(i in 1:(p+1)) { #very basic density estimation for comparison
  hist(Output$mcmc[-(1:floor(niters/2)),i], probability = TRUE)
  lines(density(Output$mcmc[-(1:floor(niters/2)),i]))
  abline(v = c(truerho, truesigma2)[i], col = 'red')
}
par(ask=FALSE) # SEH: switch off when you're done

hist(Output$mcmc[-c(1:niters/10), 1],breaks=20,probability=TRUE,
     main = 'Posterior Density for PARCOR 1',xlab=NULL,ylim=c(0,8))
lines(density(Output$mcmc[-c(1:niters/10), 1]),col=4,lwd=3)
abline(v = truerho[1], col = 'red',lwd=2)

hist(Output$mcmc[-c(1:niters/10), 2],breaks=20,probability=TRUE,
     main = 'Posterior Density for PARCOR 2',xlab=NULL)
lines(density(Output$mcmc[-c(1:niters/10), 2]),col=4,lwd=3)
abline(v = truerho[2], col = 'red',lwd=2)

hist(Output$mcmc[-c(1:niters/10), 3],breaks=20,probability=TRUE,
     main = 'Posterior Density for PARCOR 3',xlab=NULL)
lines(density(Output$mcmc[-c(1:niters/10), 3]),col=4,lwd=3)
abline(v = truerho[3], col = 'red',lwd=2)

hist(Output$mcmc[-c(1:niters/10), 4],breaks=20,probability=TRUE,
     main = 'Posterior Density for PARCOR 4',xlab=NULL)
lines(density(Output$mcmc[-c(1:niters/10), 4]),col=4,lwd=3)
abline(v = truerho[4], col = 'red',lwd=2)

hist(Output$mcmc[-c(1:niters/10), 5],breaks=20,probability=TRUE,
     main = 'Posterior Density for Variance',xlab=NULL)
lines(density(Output$mcmc[-c(1:niters/10), 5]),col=4,lwd=3)
abline(v = truesigma2, col = 'red',lwd=2)
 
plot(ts(Output$mcmc[-c(1:niters/10), 1]),main='Trace plot',ylab=NULL)
acf(ts(Output$mcmc[-c(1:niters/10), 1]),main='ACF')

plot(ts(Output$mcmc[-c(1:niters/10), 2]),main='Trace plot',ylab=NULL)
acf(ts(Output$mcmc[-c(1:niters/10), 2]),main='ACF')

plot(ts(Output$mcmc[-c(1:niters/10), 3]),main='Trace plot',ylab=NULL)
acf(ts(Output$mcmc[-c(1:niters/10), 3]),main='ACF')

plot(ts(Output$mcmc[-c(1:niters/10), 4]),main='Trace plot',ylab=NULL)
acf(ts(Output$mcmc[-c(1:niters/10), 4]),main='ACF')

plot(ts(Output$mcmc[-c(1:niters/10), 5]),main='Trace plot',ylab=NULL)
acf(ts(Output$mcmc[-c(1:niters/10), 5]),main='ACF')

#all not looking fine :(
#now, stan validation

#####################
# Stationary region #
#####################
data <- list(T = 1000,
             p = 4,
             y = as.numeric(y),
             a = 2,
             b = 2,
             alpha = c(2,3,3,2),
             beta = c(3,2,2,3))
outputStationary <- stan('stationary.stan',
                         data = data, iter=10000,
                         chains = 4) #not doing anything :D
#need to extract the PARCOR and sigma
outStat <- extract(outputStationary, c('r', 'sigma2'))
hist(out$sigma2,probability=TRUE) #very very bad
hist(out$r[,1],probability = TRUE);abline(v=truerho[1],col='red')
hist(out$r[,2],probability = TRUE);abline(v=truerho[1],col='red')
hist(out$r[,3],probability = TRUE);abline(v=truerho[1],col='red')
hist(out$r[,4],probability = TRUE);abline(v=truerho[1],col='red')

hist(Output$mcmc[-c(1:1000),1], breaks=20, probability = TRUE, ylim = c(0,8),
     main='Posterior Summary', xlab=NULL, col = alpha('purple', 0.4))
hist(outStat$r[-c(1:1000),1], breaks=20, probability = TRUE,
     xlab=NULL, add=TRUE, col = alpha('orange',0.5))
lines(density(Output$mcmc[-c(1:1000),1]), col = 'purple', lwd=3)
lines(density(outStat$r[-c(1:1000),1]), col = 'red', lwd=3)
abline(v=truerho[1],col='cyan',lwd=3)
legend('topleft', legend = c('MH density estimate','NUTS density estimate','True value'),
       fill = c('purple', 10, 'cyan'))

hist(Output$mcmc[-c(1:1000),2], breaks=20, probability = TRUE,
     main='Posterior Summary', xlab=NULL, col = alpha('purple', 0.4))
hist(outStat$r[-c(1:1000),2], breaks=20, probability = TRUE,
     xlab=NULL, add=TRUE, col = alpha('orange',0.5))
lines(density(Output$mcmc[-c(1:1000),2]), col = 'purple', lwd=3)
lines(density(outStat$r[-c(1:1000),2]), col = 'red', lwd=3)
abline(v=truerho[2],col='cyan',lwd=3)
legend('topright', legend = c('MH density estimate','NUTS density estimate','True value'),
       fill = c('purple', 10, 'cyan'))

hist(Output$mcmc[-c(1:1000),3], breaks=20, probability = TRUE,
     main='Posterior Summary', xlab=NULL, col = alpha('purple', 0.4))
hist(outStat$r[-c(1:1000),3], breaks=20, probability = TRUE,
     xlab=NULL, add=TRUE, col = alpha('orange',0.5))
lines(density(Output$mcmc[-c(1:1000),3]), col = 'purple', lwd=3)
lines(density(outStat$r[-c(1:1000),3]), col = 'red', lwd=3)
abline(v=truerho[3],col='cyan',lwd=3)
legend('topleft', legend = c('MH density estimate','NUTS density estimate','True value'),
       fill = c('purple', 10, 'cyan'))

hist(Output$mcmc[-c(1:1000),4], breaks=20, probability = TRUE,
     main='Posterior Summary', xlab=NULL, col = alpha('purple', 0.4))
hist(outStat$r[-c(1:1000),4], breaks=20, probability = TRUE,
     xlab=NULL, add=TRUE, col = alpha('orange',0.5))
lines(density(Output$mcmc[-c(1:1000),4]), col = 'purple', lwd=3)
lines(density(outStat$r[-c(1:1000),4]), col = 'red', lwd=3)
abline(v=truerho[4],col='cyan',lwd=3)
legend('topright', legend = c('MH density estimate','NUTS density estimate','True value'),
       fill = c('purple', 10, 'cyan'))

hist(Output$mcmc[-c(1:1000),5], breaks=20, probability = TRUE,
     main='Posterior Summary', xlab=NULL, col = alpha('purple', 0.4))
hist(outStat$sigma2[-c(1:1000)], breaks=20, probability = TRUE,
     xlab=NULL, add=TRUE, col = alpha('orange',0.5))
lines(density(Output$mcmc[-c(1:1000),5]), col = 'purple', lwd=3)
lines(density(outStat$sigma2[-c(1:1000)]), col = 'red', lwd=3)
abline(v=truesigma2,col='cyan',lwd=3)
legend('topright', legend = c('MH density estimate','NUTS density estimate','True value'),
       fill = c('purple', 10, 'cyan'))
