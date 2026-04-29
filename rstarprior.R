library(invgamma)
library(mvtnorm)
library(scales)

library(rstan)
#library(durhamSLR)
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

r2phi <- function(r){
  p = length(r)
  mat = diag(as.numeric(r)) #rkk = rk
  for(k in 2:p){
    for(i in 1:(k-1)){
      mat[i,k] = mat[i,(k-1)]-r[k]*mat[(k-i),(k-1)]
    }
  }
  return(as.numeric(mat[,p]))
}


#########
# model #
#########

p=4
truer=c(0.3,-0.2,-0.7,0.4)
truerst= log(1+truer)-log(1-truer)
truesigma2=4
y=arima.sim(n=1000,list(ar=r2phi(c(0.3,-0.2,-0.7,0.4))),sd=2)
plot(y)

########
# stan #
########

data=list(T = 1000,
          p = 4,
          y = as.numeric(y),
          a = 2,
          b = 2,
          s=1)
outputStationaryp <- stan('stationarymodifiedprior.stan',
                         data = data, iter=10000,
                         chains = 4)
outStat <- extract(outputStationaryp, c('rst', 'sigma2'))

##########
# Graphs #
##########

#the acfs are trivial and not included in the
report.
plot(ts(outStat$rst[,1]),main='Trace plot',ylab=NULL)
acf(ts(outStat$sigma2),main='ACF')
plot(ts(outStat$rst[,2]),main='Trace plot',ylab=NULL)
acf(ts(outStat$rst[,2]),main='ACF')
plot(ts(outStat$rst[,3]),main='Trace plot',ylab=NULL)
acf(ts(outStat$rst[,3]),main='ACF')
plot(ts(outStat$rst[,4]),main='Trace plot',ylab=NULL)
acf(ts(outStat$rst[,4]),main='ACF')
plot(ts(as.numeric(outStat$sigma2)),main='Trace plot',ylab=NULL)
acf(ts(outStat$sigma2),main='ACF')

hist(outStat$sigma2,breaks=20,probability=TRUE,
     main='Posterior Density Estimate',xlab='sigma2',col=alpha('orange',0.4))
lines(density(outStat$sigma2),col='red',lwd=2)
abline(v=truesigma2,col=2,lwd=2)
hist(outStat$rst[,1],breaks=20,probability=TRUE,
     main='Posterior Density Estimate',xlab='rst1',col=alpha('orange',0.4))
lines(density(outStat$rst[,1]),col='red',lwd=2)
abline(v=truerst[1],col=2,lwd=2)
hist(outStat$rst[,2],breaks=20,probability=TRUE,
     main='Posterior Density Estimate',xlab='rst2',col=alpha('orange',0.4))
lines(density(outStat$rst[,2]),col='red',lwd=2)
abline(v=truerst[2],col=2,lwd=2)
hist(outStat$rst[,3],breaks=20,probability=TRUE,
     main='Posterior Density Estimate',xlab='rst3',col=alpha('orange',0.4))
lines(density(outStat$rst[,3]),col='red',lwd=2)
abline(v=truerst[3],col=2,lwd=2)
hist(outStat$rst[,4],breaks=20,probability=TRUE,
     main='Posterior Density Estimate',xlab='rst4',col=alpha('orange',0.4))
lines(density(outStat$rst[,4]),col='red',lwd=2)
abline(v=truerst[4],col=2,lwd=2)
