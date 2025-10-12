'''
AR(1) process:
y_t = phi*y+t-1 + epsilon_t
we will take epsilon_t ~ Normal(0, v) in the first instance and go from there

'''

t <- seq(0, 10, 0.01)
y <- rep(0, length(t))

phi = 0.3
v=1
set.seed(1)
for(i in 2:length(t)){
  y[i] = rnorm(1, phi*y[i-1], v)
}
plot(t,y, type = 'l')

#make this a function, want to vary phi and v
series_reg <- function(phi, v) {
  y <- rep(0, length(t))
  for(i in 2:length(t)){
    y[i] = rnorm(1, phi*y[i-1], v)
  }
  return(y)
}
y.data <- data.frame("t" = t,
                     "ts1" = series_reg(0.1, 1),
                     "ts2" = series_reg(0.4, 1),
                     "ts3" = series_reg(0.7, 1),
                     "ts4" = series_reg(-0.2, 1))


ggplot2::ggplot(data = y.data) +
  ggplot2::geom_line(mapping = ggplot2::aes(x = t, y = ts1, color = "red")) +
  ggplot2::geom_line(mapping = ggplot2::aes(x = t, y = ts2, color = "blue")) +
  ggplot2::geom_line(mapping = ggplot2::aes(x = t, y = ts3, color = "green")) +
  ggplot2::geom_line(mapping = ggplot2::aes(x = t, y = ts4, color = "purple"))
  

#want a general AR(p) process to model, requires effort

#this is going to be one iteration of the time series i.e. y_t = f(y_t-1), 
#then another function will apply that to data

arp_process <- function(phi, xt1){
  #want to use state space representation
  p <- length(phi)
  
  g <- diag(x=1, nrow = p, ncol = p)[-p,]
  G <- rbind(phi, g)
  f <- c(1, rep(0, p-1))
  
  omegat <- c(rnorm(1, 0, v), rep(0,p-1))
  
  xt <- G%*%xt1 + omegat
  yt <- t(f)%*%xt
  
  rm(list = c('p', 'g', 'G', 'f', 'xt', 'omegat'))
  
  return(yt)
}
#need some fake data
arp_process(rep(0.1, 4), c(5,4,3,6))
#good thing: this makes a functional style operation possible which means we can use apply()

#we want fake data now - would have to specify the first p values of the series then let it go

yt <- c(2,3,2,3,4,5, rep(0,995))

for(i in 7:length(yt)){
  yt[i] <- arp_process(rep(0.1,6), yt[(i-1):(i-6)])
}
plot(ts(yt))
acf(yt)


#we have made some headway into the actual generation of the model. Now we turn to inference.