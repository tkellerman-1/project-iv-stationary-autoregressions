library(rstan)
#library(durhamSLR)
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)
library(scales)

rtop <- function(r){
  p = length(r)
  mat = diag(r, nrow = p) #rkk = rk
  for(k in 2:p){
    for(i in 1:(k-1)){
      mat[i,k] = mat[i,(k-1)]-r[k]*mat[(k-i),(k-1)]
    }
  }
  return(as.numeric(mat[,p]))
}

#######################################
# Stationary with unconstrained prior #
#######################################



####################################
# random walk in stationary region #
####################################

#true variances
sigw <- 0.1
sigr <- 0.1
p <- 3
T <- 1000
#random walk in Rp - independent atm so can just use rnorm
Y <- c(0.2,0.5,0.3,rep(0,T-p))
rst=matrix(0,1000,3)
set.seed(23) #30 4 41 42 23 54
rst0 <- rnorm(p, 0 , sqrt(sigw))
rst[1,] = rst0
for(t in (p+1):T) {
  for(i in 1:p){
    rst[t,i] = rnorm(1,rst[(t-1),i],sqrt(sigw)) 
  }
  rt = (exp(rst[t,]) - 1) / (exp(rst[t,]) + 1)
  phit = rtop(rt)
  Y[t] = rnorm(1, mean = phit %*% Y[(t-1):(t-p)], sd = sqrt(sigr))
}

plot(ts(Y),main = 'Time-Varying AR(3)',ylab = NULL, col = 4,lwd=0.3)
plot(ts(Y),xlim = c(900,1000),ylab = NULL,
       main='Time-Varying AR(3)',col=4,lwd=2)

data <- list(p = 3,
             T = 1000,
             Y = Y,
             a = 2,
             b = 2,
             alpha = 2,
             beta = 2,
             s = 0.1)
outputSRW <- stan('stationarySRW.stan', iter=10000,
               data = data, chains = 4)
pars <- extract(outputSRW, c("sigmaA", "sigmaR"))
plot(pars[[1]],type = 'l') #not too bad mixing
plot(pars[[2]],type = 'l') #same
#they look like they are autocorrelated rather than white noise but could be worse
hist(pars[[1]], breaks = 20)

density(pars[[1]]) #mean is 0.108!!!! yayyyyyyy
density(pars[[2]])
#this works great icl

############################
# Graphs for stationary RW #
############################
#error variance
hist(pars$sigmaA, breaks = 20, probability = TRUE,
     main = 'Posterior density of sigmaA', xlab = NULL)
lines(density(pars$sigmaA), col = 4, lwd = 3)
abline(v = sigr, col = 'red',lwd=2)

hist(pars$sigmaR, breaks = 20, probability = TRUE,
     main = 'Posterior density of sigmaR', xlab = NULL)
lines(density(pars$sigmaR), col = 4, lwd = 3)
abline(v = sigw, col = 'red',lwd=2)



#now look to validate the last 20 steps
data <- list(p=3,
             T=980,
             Y=Y[1:980],
             a=2,
             b=2,
             alpha=2,
             beta=2,
             s=0.1,
             Tstar = 20)
outputSRWPP <- stan('stationarySRWPP.stan', iter=10000,
                data=data,chains=4)
parspp <- extract(outputSRWPP, c("sigmaA", "sigmaR"))
plot(parspp[[1]],type = 'l', ylab='sigmaA',main='Trace plot') #not too bad mixing
plot(parspp[[2]],type = 'l', ylab='sigmaR',main='Trace plot') #same

parspp <- extract(outputSRWPP)
preds <- apply(parspp$Ystar,2,mean)[p:23]
predCIu <- apply(parspp$Ystar, 2, quantile, p=0.975)[p:23]
predCId <- apply(parspp$Ystar, 2, quantile, p=0.025)[p:23]
str(preds)
plot(ts(Y[1:980]), xlim = c(960,1000), lwd = 2,ylab = NULL,
     main='Time-Varying AR(3)')
lines(980:1000, predCIu, col = 'purple',lwd=1.5)
lines(980:1000, predCId, col = 'purple',lwd=1.5)
polygon(c(980:1000,1000:980),c(predCIu, rev(predCId)),col = alpha('purple',0.4),lty=0)
lines(980:1000, Y[980:1000],col=4,lwd=2)
lines(980:1000, preds, col = 'red', lwd = 2)

legend('bottomleft', legend=c('Mean prediction',
                            '95% credible interval',
                            'True value'),
       fill = c('red','purple',4))

###############################################################################
# Random walk in Rp, non-stationary AR coefficients with predictive densities #
###############################################################################

data <- list(p=3,
             T=980,
             Y=Y[1:980],
             a=2,
             b=2,
             alpha=2,
             beta=2,
             s=0.1,
             Tstar = 20)
Outputnonstationary <- stan('nonstationarySRWpp.stan', data=data,
                            iter=10000, chains=4)
parsppns <- extract(Outputnonstationary)
plot(parsppns[[1]],type = 'l', ylab='sigmaA', main='Trace plot') #not too bad mixing
plot(parsppns[[2]],type = 'l', ylab='sigmaR', main='Trace plot') #same

hist(parsppns$sigmaR,probability=TRUE,breaks=20,
     main = 'Posterior summary', xlab='sigmaA')
lines(density(parsppns$sigmaR),col='blue',lwd=2)
abline(v=0.1,col='red',lwd=1.5)

str(parspp$Ystar)
preds.ns <- apply(parsppns$Ystar,2,mean)[p:23]
predCIu.ns <- apply(parsppns$Ystar, 2, quantile, p=0.975)[p:23]
predCId.ns <- apply(parsppns$Ystar, 2, quantile, p=0.025)[p:23]
str(preds.ns)
plot(ts(Y[1:980]), xlim = c(970,1000), lwd = 2,ylab = NULL,
     main='Time-Varying AR(3)')
lines(980:1000, predCIu.ns, col = 'brown',lwd=1.5)
lines(980:1000, predCId.ns, col = 'brown',lwd=1.5)
polygon(c(980:1000,1000:980),c(predCIu.ns, rev(predCId.ns)),col = alpha('brown',0.2),lty=0)
lines(980:1000, Y[980:1000],col=4,lwd=2)
for(i in 2:23){
  boxplot(parspp$Ystar[,i], add=TRUE,col=alpha('purple',0.4),at=(980+i),range=0)
}
boxplot(parspp$Ystar[,3], add=TRUE,at=981)
boxplot(parspp$Ystar[,6], add=FALSE,col=alpha('purple',0.4))

lines(980:1000, predCIu, col = 'purple',lwd=1.5)
lines(980:1000, predCId, col = 'purple',lwd=1.5)
polygon(c(980:1000,1000:980),c(predCIu, rev(predCId)),col = alpha('purple',0.4),lty=0)
lines(980:1000, Y[980:1000],col=4,lwd=2)

lines(980:1000, preds.ns, col = 'green', lwd = 2) #completely expected if the RW goes out of stationary region, which it has
lines(980:1000, preds, col = 'red', lwd = 2)


legend('bottomleft', legend=c('Stationary mean pred',
                              'Stationary 95% CrI',
                              'Non-stationary mean pred',
                              'Non-stationary 95% CrI',
                              'True value'),
       fill = c('red','purple','green','brown',4))



