
# INPUTS
n=25  #number of observations
N=2500  #particle sample size
true.mu = 0
sigma = 1
mu.hyper = 0
sigma.hyper = 10
data=rnorm(n,true.mu,sigma)
epsilon = 1
time.steps = 20
weights = matrix(1/N,time.steps,N)
mu=matrix(NA,time.steps,N)
d=matrix(NA,time.steps,N)
rho=function(y,x) abs(sum(y)-sum(x))/n

for(t in 1:time.steps){
  if(t==1){
    for(i in 1:N){
      d[t,i]= epsilon +1
      while(d[t,i]>epsilon) {
        proposed.mu=rnorm(1,0,sigma.hyper) #<--prior draw
        x=rnorm(n, proposed.mu, sigma)
        d[t,i]=rho(data,x)
        }
      mu[t,i]= proposed.mu
    }
    }
  else{
      epsilon = c(epsilon,quantile(d[t-1,],.75))
      mean.prev <- sum(mu[t-1,]*weights[t-1,])
      var.prev <- sum((mu[t-1,] - mean.prev)^2*weights[t-1,])
      for(i in 1:N){
        d[t,i]= epsilon[t]+1
      while(d[t,i]>epsilon[t]) {
        sample.particle <- sample(N, 1, prob = weights[t-1,])
        proposed.mu0 <- mu[t-1, sample.particle]
        proposed.mu <- rnorm(1, proposed.mu0, sqrt(2*var.prev))
        x <- matrix(rnorm(n,proposed.mu, sigma),n,1)
        d[t,i]=rho(data,x) 
        }
      mu[t,i]= proposed.mu
      mu.weights.denominator<-
        sum(weights[t-1,]*dnorm(proposed.mu,mu[t-1,],sqrt(2*var.prev)))
      mu.weights.numerator<-dnorm(proposed.mu,0,sigma.hyper)
      weights[t,i] <- mu.weights.numerator/mu.weights.denominator
      }
      }
  weights[t,] <- weights[t,]/sum(weights[t,])
  }