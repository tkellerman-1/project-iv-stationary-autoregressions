library(tvReg)
library(rstan)
#library(durhamSLR)
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)
library(scales)

########
# Data #
########
?tvReg::RV

data(RV)
#concentrate on the RV data set
head(RV)

plot(ts(RV$RV),main='Daily Realised Variance for S&P 500, 1990-2006',
     col = 'magenta', ylab=NULL)
plot(ts(RV$RQ_lag_sqrt))
RVlog=log(RV$RV)
RVlog=RVlog-mean(RVlog)
plot(ts(RVlog), main = 'Zero-centred RV for S&P 500, 1990-2006',
     col='magenta',ylab=NULL) #centres on 0


#########################################
# Model order via predictive capability #
#########################################

#run a much smaller model for not many iterations, up to 7 times to see if there's any influence from weeks

thinRV=rep(0,1000)
for(t in 1:1000){
  thinRV[t]=RVlog[4*t]
}
plot(ts(thinRV)) #if the stan takes too long in the first instance
thinRV=as.numeric(thinRV)[1:980]

#This is a very basic prediction validation scheme.
#To actually select a model, we will use mean squared error.
standata1 <- list(p=1,       
                 T=980,
                 Y=thinRV,
                 a=2,
                 b=3,
                 alpha=3,
                 beta=2,
                 s=0.05,
                 kappap=1,
                 Tstar=20)

standata2=standata1
standata3=standata1
standata4=standata1
standata5=standata1
standata6=standata1
standata7=standata1

standata2[[1]]=2
standata3[[1]]=3
standata4[[1]]=4
standata5[[1]]=5
standata6[[1]]=6
standata7[[1]]=7

applicationstan1 <- stan('stationarySRWPP.stan', data=standata1,
                        iter=1000, chains=1) #Rhats dont work
applicationstan2 <- stan('stationarySRWPP.stan', data=standata2,
                         iter=1000, chains=1) #chain mixes better ish
applicationstan3 <- stan('stationarySRWPP.stan', data=standata3,
                         iter=1000, chains=1) #
applicationstan4 <- stan('stationarySRWPP.stan', data=standata4,
                         iter=1000, chains=1)
applicationstan5 <- stan('stationarySRWPP.stan', data=standata5,
                         iter=1000, chains=1)
applicationstan6 <- stan('stationarySRWPP.stan', data=standata6,
                         iter=1000, chains=1)
applicationstan7 <- stan('stationarySRWPP.stan', data=standata7,
                         iter=1000, chains=1)



#PLOT: DO NOT REDRAW
thinRV=rep(0,1000)
for(t in 1:1000){
  thinRV[t]=RVlog[4*t]
}
plot(ts(thinRV[1:980]), xlim=c(975,1000),ylab=NULL, main='Predictive densities')
lines(980:1000,thinRV[980:1000],col='red',lwd=2)
legend('topleft', legend = c(1,2,3,4,5,6,7), 
       lwd=rep(1.5,7),col = c('red','orange','yellow','green','blue','purple','violet'))

##Predictions - manually changed p and applicationstan*p* and the colour of the credible intervals
p=1
model=applicationstan1
useful=extract(model, c("Ystar"))
preds=apply(useful$Ystar,2,mean)
predCIu=apply(useful$Ystar, 2, quantile, p=0.975)
predCId=apply(useful$Ystar, 2, quantile, p=0.025)
lines(980:1000, predCIu[p:(20+p)], col = alpha('violet',alpha=0.8),lwd=1.5)
lines(980:1000, predCId[p:(20+p)], col = alpha('violet',alpha=0.8),lwd=1.5)


#increasing p by 1 adds about 16MB of data for 1000 iters so expect approx 640MB per 10000 iters * 4 chains
#time taken from using whole chain will be absolutely vast.


#prediction intervals get smaller up to about p=4 then it kind of levels off

#mean square prediction error - should all basically be the same
#again, manually changed p and computed again
p=7
model=applicationstan7
useful=extract(model, c("Ystar"))
preds=apply(useful$Ystar,2,mean)
mse7=mean((thinRV[(1000-(20+p-1)):1000]-preds)^2)

########################
# Full model inference #
########################

#total chain takes roughly 20GB to load in memory all at once,
#I only have 16GB thus only used thinned data set.

#full model inference takes about 90-100 hours and the first attempt had 20000 divergent transitions
#so were not doing that again unless we have something satisfactory

#the two options are to run a thinRV for the full 20000 samples
#or run the full T=4000 time series for not very long
#we'll use the first in the full report because that is a more robust estimate, even if its not strictly
#useful for predictions using this data set

thinRV=rep(0, floor(4264/4))
for(i in 1:floor(4264/4)){
  thinRV[i]=RVlog[4*i]
}
validation=thinRV[1046:1066]
thinRV=thinRV[1:1046] #withhold last 20
plot(ts(thinRV),col='magenta',main='Thinned log zero-centred daily RV',ylab=NULL) 

fullstandata <- list(p=7,       
                T=1046,
                Y=thinRV,
                a=2,
                b=3,
                alpha=3,
                beta=2,
                s=0.05,
                kappap=1,
                Tstar=20)
applicationstan <- stan('stationarySRWPP.stan',data=fullstandata,
                        iter=10000,chains=4)
hmcout=extract(applicationstan, c("sigmaR", "sigmaA", "kappa"))
predsout=extract(applicationstan, "Ystar")
str(as.numeric(hmcout$sigmaR))
plot(ts(as.numeric(hmcout$sigmaR)),main='Trace plot',ylab='sigmaR') #odd
plot(ts(as.numeric(hmcout$sigmaA)),main='Trace plot',ylab='sigmaA') #fine
plot(ts(as.numeric(hmcout$kappa)),main='Trace plot',ylab='kappa') #odd

str(density(as.numeric(hmcout$sigmaR)))

optimiser=function(x,y){
  xmax=x[1]
  ymax=x[2]
  N=length(y)
  for(i in 1:N){
      if(y[i]>ymax){
        xmax=x[i]
        ymax=y[i]
      }
  }
  return(c(xmax,ymax))
}
dsigmaR=density(as.numeric(hmcout$sigmaR))
dsigmaA=density(as.numeric(hmcout$sigmaA))
dkappa=density(as.numeric(hmcout$kappa))

sigmaRMAP=optimiser(dsigmaR$x,dsigmaR$y)
sigmaAMAP=optimiser(dsigmaA$x,dsigmaA$y)
kappaMAP=optimiser(dkappa$x,dkappa$y)

#quantile applied to underlying data produces same quantiles
#as the quantile.density function. Could have also just done this to get the MAP anyway
quant=function(x){
  return(c(quantile(x,p=0.025),quantile(x,p=0.975)))
}
sigmaRcrI=quant(as.numeric(hmcout$sigmaR))
sigmaAcrI=quant(as.numeric(hmcout$sigmaA))
kappaCrI=quant(as.numeric(hmcout$kappa))

hist(as.numeric(hmcout$sigmaR),probability=TRUE,breaks=20,
     main = 'Posterior Density Estimate',xlab='sigmaR')
lines(density(as.numeric(hmcout$sigmaR)),lwd=2,col='red')
abline(v=sigmaRMAP[1],col='blue',lwd=2)
points(sigmaRMAP[1],sigmaRMAP[2],pch=16,col='blue')
legend('topright', legend=c('Density estimate', 'MAP estimate'), lty=c(1,0),pch=c(1,16),col=c('red','blue'))
#abline(v=quantile(as.numeric(hmcout$sigmaR),p=0.025),col='purple')
#abline(v=quantile(as.numeric(hmcout$sigmaR),p=0.975),col='purple')
#polygon(c(sigmaRcrI[1],sigmaRcrI[1],sigmaRcrI[2],sigmaRcrI[2]), c(0,30,30,0),col=alpha('pink',0.3))
lines(sigmaRcrI, rep(0,2), col='purple',lwd=4)
legend('topright', legend=c('Density estimate', 'MAP estimate', '95% Credible interval'),
       lty=c(1,0,1),pch=c(NULL,16,NULL),col=c('red','blue','purple'),lwd=c(2,0,4))


hist(as.numeric(hmcout$sigmaA),probability=TRUE,breaks=20,
     main = 'Posterior Density Estimate',xlab='sigmaA')
lines(density(as.numeric(hmcout$sigmaA)),lwd=2,col='red')
abline(v=sigmaAMAP[1],col='blue',lwd=2)
points(sigmaAMAP[1],sigmaAMAP[2],pch=16,col='blue')
lines(sigmaAcrI, rep(0,2), col='purple',lwd=4)
legend('topleft', legend=c('Density estimate', 'MAP estimate', '95% Credible interval'),
       lty=c(1,0,1),pch=c(NULL,16,NULL),col=c('red','blue','purple'),lwd=c(2,0,4))



hist(as.numeric(hmcout$kappa),probability=TRUE,breaks=20,
     main = 'Posterior Density Estimate',xlab='kappa')
lines(density(as.numeric(hmcout$kappa)),lwd=2,col='red')
abline(v=kappaMAP[1],col='blue',lwd=2)
points(kappaMAP[1],kappaMAP[2],pch=16,col='blue')
lines(kappaCrI, rep(0,2), col='purple',lwd=4)
legend('topleft', legend=c('Density estimate', 'MAP estimate', '95% Credible interval'),
       lty=c(1,0,1),pch=c(NULL,16,NULL),col=c('red','blue','purple'),lwd=c(2,0,4))



preds=apply(as.matrix(predsout$Ystar,ncol=27), 2, mean)[7:27]
predsCIu=apply(as.matrix(predsout$Ystar,ncol=27),2, quantile,p=0.975)[7:27]
predsCId=apply(as.matrix(predsout$Ystar,ncol=27),2, quantile,p=0.025)[7:27]
length(predsCIu) 
preds
plot(thinRV,xlim=c(1000,1066),type='l')
lines(1046:1066, predsCIu, lwd=1.5,col='purple')
lines(1046:1066,predsCId, lwd=1.5,col='purple')
polygon(c(1046:1066,1066:1046),c(predsCIu, rev(predsCId)),col = alpha('purple',0.4),lty=0)
lines(1046:1066,preds,lwd=2,col='red')
lines(1046:1066, validation,col='blue',lwd=1.5)
legend('topleft', legend=c('Mean Prediction','95% Credible Interval', 'True Value'),
       lty=c(1,1,1), lwd=c(2,1.5,1.5), col=c('red','purple','blue'))

truethinRV=rep(0, floor(4264/4))
for(i in 1:floor(4264/4)){
  truethinRV[i]=RV$RV[4*i]
}
truevalidation=truethinRV[1046:1066]
truepreds=exp(preds+mean(log(RV$RV)))
plot(ts(log(RV$RV)-mean(log(RV$RV))))
truepredsCIu=exp(predsCIu+mean(log(RV$RV)))
truepredsCId=exp(predsCId+mean(log(RV$RV)))

plot(ts(truethinRV),xlim=c(1026,1066),ylim=c(0,2e-4),ylab='DRV',main='Predictions for daily realised variance')
lines(1046:1066, truepredsCIu, lwd=1.5,col='purple')
lines(1046:1066, truepredsCId, lwd=1.5,col='purple')
polygon(c(1046:1066,1066:1046),c(truepredsCIu, rev(truepredsCId)),col = alpha('purple',0.4),lty=0)
lines(1046:1066,truevalidation,lwd=1.5,col='blue')
lines(1046:1066, truepreds,lwd=2,col='red')
legend('topleft', legend=c('Mean Prediction','95% Credible Interval', 'True Value'),
       lty=c(1,1,1), lwd=c(2,1.5,1.5), col=c('red','purple','blue'))

