##################################
### Graphs for AR(p) processes ###
##################################

#arima.sim or other methods don't generate non-stationary time series so this is done manually.

T=1000

########################
# Non-stationary AR(1) #
########################

y1=rep(0,T)
set.seed(44)
for(t in 2:T){
  y1[t] = rnorm(1,1.1*y1[t-1], 1)
}

plot(ts(y1),col=3,lwd=2, main = 'AR(1) Model 1',ylab=NULL,xlab='Time')

########################
# Stationary AR(1)     #
########################

y2=rep(0,T)
set.seed(44)
for(t in 2:T){
  y2[t] = rnorm(1,0.9*y2[t-1], 1)
}

plot(ts(y2),col=4,lwd=1.5, main = 'AR(1) Model 2',ylab=NULL,xlab='Time')

########################
# Random Walk AR(1)    #
########################

y3=rep(0,T)
set.seed(44)
for(t in 2:T){
  y3[t] = rnorm(1,y3[t-1], 1)
}

plot(ts(y3),col=2,lwd=1.5, main = 'AR(1) Model 3',ylab=NULL,xlab='Time')


#########################
# Stationary region     #
#########################

#Transformation [-1,1]^p -> C_p
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

#sample
samp=runif(10000,-1,1)
sampmat=matrix(samp,nrow=2000,ncol=5)
statsamp=matrix(0,2000,5)
for(i in 1:nrow(statsamp)){
  statsamp[i,]=r2phi(sampmat[i,])
}
pairs(statsamp,col='2',pch=20, labels=c(1:5))
