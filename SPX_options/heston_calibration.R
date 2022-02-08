### Packages ###

library(tidyverse)

### Functions ###

charFunc <- function(u, spot, r=0, div=0, ttm, v, sigma, kappa, theta, rho){
  d <- sqrt( (rho*sigma*u*1i-kappa)^2 + sigma^2*(u*1i+u^2) )
  d <- -d
  g <- (kappa - rho*sigma*u*1i + d) / (kappa - rho*sigma*u*1i - d)
  tempM <- (kappa-rho*sigma*u*1i+d)*ttm-2*log((1-g*exp(d*ttm))/(1-g))
  M <- (r-div)*u*1i*ttm+(kappa*theta)/(sigma^2)*tempM
  N <- (kappa-rho*sigma*u*1i+d)/(sigma^2)*((1-exp(d*ttm))/(1-g*exp(d*ttm)))
  res <- exp(M+N*v+1i*u*log(spot))
  return(res)
}

HestonCall <- function(spot, strike, r=0, div=0, ttm, v, sigma, kappa, theta, rho){
  integrand1 <- function(u){ 
    num1 <- charFunc(u-1i, spot, r, div, ttm, v, sigma, kappa, theta, rho)
    den1 <- charFunc(-1i, spot, r, div, ttm, v, sigma, kappa, theta, rho)
    dummy1 <- exp(-1i*u*log(strike))*num1/(1i*u*den1)
    integrand1 <- Re(dummy1)
  }
  integrand2 <- function(u){ 
    dummy2 <- exp(-1i*u*log(strike))*charFunc(u, spot, r, div, ttm, v, sigma, kappa, theta, rho)/(1i*u)
    integtand2 <- Re(dummy2)
  }
  Pi1 <- 0.5 + 1/pi * integrate(integrand1,0,Inf,stop.on.error = FALSE)$value
  Pi2 <- 0.5 + 1/pi * integrate(integrand2,0,Inf,stop.on.error = FALSE)$value
  res <- spot*exp(-div*ttm)*Pi1 - strike*exp(-r*ttm)*Pi2
  return(res)
}

simulateHeston <- function(S0,nPaths,nDates,r,ttm,v,sigma,kappa,theta,rho){
  
  nSim <- nPaths/2
  dt <- ttm/nDates
  Smat <- matrix(nrow = nPaths, ncol = nDates)
  Vmat <- matrix(nrow = nPaths, ncol = nDates)
  X <- log(rep(S0, nPaths))
  V <- rep(v, nPaths)
  
  for(i in 1:nDates){
    
    z1 <- rnorm(nSim); z2 <- rnorm(nSim)
    dW1 <- sqrt(dt)*z1; dW2 <- rho*dW1+sqrt(1-rho^2)*sqrt(dt)*z2
    dW1_hat <- sqrt(dt)*(-z1); dW2_hat <- rho*dW1_hat+sqrt(1-rho^2)*sqrt(dt)*(-z2)
    
    X[1:nSim] <- X[1:nSim] + (r-0.5*pmax(V[1:nSim],0)) * dt + sqrt(pmax(V[1:nSim],0)) * dW1
    X[(nSim+1):(nSim*2)] <- X[(nSim+1):(nSim*2)] + (r-0.5*pmax(V[(nSim+1):(nSim*2)],0)) * dt + sqrt(pmax(V[(nSim+1):(nSim*2)],0)) * dW1_hat 
    
    V[1:nSim]<-V[1:nSim]+kappa*(theta-pmax(V[1:nSim],0))*dt+sigma*sqrt(pmax(V[1:nSim],0))*dW2
    V[(nSim+1):(2*nSim)]<-V[(nSim+1):(2*nSim)]+kappa*(theta-pmax(V[(nSim+1):(2*nSim)],0))*dt+sigma*sqrt(pmax(V[(nSim+1):(2*nSim)],0))*dW2_hat
    
    Smat[,i] <- exp(X)
    Vmat[,i] <- V
  }
  
  result <- list(Smat, Vmat)
  
}

BSM <- function(S,K,r,sigma,TTM,call=T){
  d1 <- (log(S/K) + (r + sigma^2/2) * TTM)/(sigma * sqrt(TTM))
  d2 <- d1 - sigma * sqrt(TTM)
  if(call) price <- pnorm(d1) * S - exp(-r * TTM) * K * pnorm(d2)
  else price <- exp(-r * TTM) * K * pnorm(-d2) - pnorm(-d1) * S 
  return(price)
}

BSMDelta <- function(S,K,r,sigma,TTM){
  d1 <- (log(S/K) + (r + sigma^2/2) * TTM)/(sigma * sqrt(TTM))
  return(pnorm(d1))
}

BSMVega <- function(S,K,r,sigma,TTM){
  d1 <- (log(S/K) + (r + sigma^2/2) * TTM)/(sigma * sqrt(TTM))
  return(S * dnorm(d1) * sqrt(TTM))
}

impVol <- function(price,S,K,r,TTM){
  f <- function(sigma) price - BSM(S,K,r,sigma,TTM)
  impVol <- uniroot(f,c(0,2))
  return(impVol$root)
}

### Options on Google, provided for assignment in Financial Engineering @ CBS ###

data <- read_delim('data/GOOGLData.csv', delim = ";")

l <- length(data$Price)

# Spot is(/was) equal to 2786, assume risk free rate of zero

S0 <- 2786; r <- 0

### Calculating IV ###

for(i in 1:l){
  data[i,'impliedVolatility'] <- impVol(data[[i,1]],S0,data[[i,2]],r,data[[i,3]]/365)
}

### Plotting IV ###

jpeg('google_IV.jpg', width = 1920, height = 1080, res = 300)

impVolPlot <- ggplot(data = data , mapping = aes(log(Strike/S0), impliedVolatility)) + 
  geom_line(aes(colour = factor(Expiry))) + xlab('log(K/S)') + ylab('Implied volatility') +
  labs(colour = "Expiry (in days)"); impVolPlot

dev.off()

### Delta and Vega ###

for(i in 1:l){
  data[i, 'delta'] <- BSMDelta(S0,data[[i,2]],r,data[[i,4]],data[[i,3]]/365)
  data[i, 'vega'] <- BSMVega(S0,data[[i,2]],r,data[[i,4]],data[[i,3]]/365) / 100
}

deltaPlot <- ggplot(data = data , mapping = aes(log(Strike/S0), delta)) + 
  geom_line(aes(colour = factor(Expiry))) + xlab('log(K/S)') + ylab('Delta') + 
  labs(colour = "Expiry (in days)"); deltaPlot

vegaPlot <- ggplot(data = data , mapping = aes(log(Strike/S0), vega)) + 
  geom_line(aes(colour = factor(Expiry))) + xlab('log(K/S)') + ylab('Vega') +
  labs(colour = "Expiry (in days)"); vegaPlot 

### Calibrating Heston to observed prices ###

# We calibrate by minimizing squared sum of errors
lossFunction <- function(parms){
  sum <- 0
  for(i in 1:l){
    sum <- sum + ( (data[[i,1]] - HestonCall(S0,data[[i,2]],0,0,data[[i,3]]/365,parms[1],parms[2],parms[3],parms[4],parms[5])) / (data[[i,6]] * 100) )^2
  }
  return(sum)
}

# Providing logically sound bounds on paramters speeds up the process (a lot)
par <- optim(par = c(0.2^2, 1, 2, 0.2^2, -0.25), lossFunction, method = "L-BFGS-B",
             lower = c(0.01,0.01,0.01,0.01,-1), upper = c(0.99,10,10,0.99,1))


# Rounded estimates
v <- 0.065
sigma <- 1
kappa <- 2
theta <- 0.09
rho <- -0.3

### Calculating and plotting implied volatility of prices calculated with calibrated Heston model ###

for(i in 1:l){
  HestonPrice <- HestonCall(S0,data[[i,2]],0,0,data[[i,3]]/365,v,sigma,kappa,theta,rho)
  data[i,'HestonImp'] <- impVol(HestonPrice,S0,data[[i,2]],r,data[[i,3]]/365)
}

impVolPlot2 <- ggplot(data = data , mapping = aes(log(Strike/S0), impliedVolatility)) + 
  geom_line(aes(x = log(Strike/S0), y = impliedVolatility, colour = factor(Expiry))) + 
  geom_line(aes(x = log(Strike/S0), y = HestonImp, colour = factor(Expiry)), linetype="dashed") +
  xlab('log(K/S)') + 
  ylab('Implied volatility') + 
  labs(colour = "Expiry"); impVolPlot2 


