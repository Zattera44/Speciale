charFunc <- function(u, spot, ttm, v, sigma, kappa, theta, rho, lambda, mu, delta){
  d <- sqrt( (rho*sigma*u*1i-kappa)^2 + sigma^2*(u*1i+u^2) )
  d <- -d
  g <- (kappa - rho*sigma*u*1i + d) / (kappa - rho*sigma*u*1i - d)
  C <- (kappa*theta)/(sigma^2)*((kappa-rho*sigma*u*1i+d)*ttm-2*log((1-g*exp(d*ttm))/(1-g)))
  D <- (kappa-rho*sigma*u*1i+d)/(sigma^2)*((1-exp(d*ttm))/(1-g*exp(d*ttm)))
  M <- lambda*ttm*((exp(1i*u*mu-0.5*u^2*delta^2)-1)-1i*u*(exp(mu+0.5*delta^2)-1))
  res <- exp(M+C+D*v+1i*u*log(spot))
  return(res)
}

BatesCall <- function(spot, strike, ttm, v, sigma, kappa, theta, rho, lambda, mu, delta){
  integrand1 <- function(u){ 
    num1 <- charFunc(u-1i, spot, ttm, v, sigma, kappa, theta, rho, lambda, mu, delta)
    den1 <- charFunc(-1i, spot, ttm, v, sigma, kappa, theta, rho, lambda, mu, delta)
    dummy1 <- exp(-1i*u*log(strike))*num1/(1i*u*den1)
    integrand1 <- Re(dummy1)
  }
  integrand2 <- function(u){ 
    dummy2 <- exp(-1i*u*log(strike))*charFunc(u, spot, ttm, v, sigma, kappa, theta, rho, lambda, mu, delta)/(1i*u)
    integtand2 <- Re(dummy2)
  }
  Pi1 <- 0.5 + 1/pi * integrate(integrand1,0,Inf,stop.on.error = FALSE)$value
  Pi2 <- 0.5 + 1/pi * integrate(integrand2,0,Inf,stop.on.error = FALSE)$value
  res <- spot*Pi1 - strike*Pi2
  return(res)
}

lossFunction <- function(parms){
  sum <- 0
  for(i in 1:l){
    sum <- sum + ((data$Price[i] - BatesCall(S0,data$Strike[i],data$Expiry[i]/365,
                                             parms[1],parms[2],parms[3],
                                             parms[4],parms[5],parms[6],
                                             parms[7],parms[8])) / (data$Vega[i] * 100) )^2
  }
  return(sum)
}

BSM <- function(S,K,sigma,TTM,call=T){
  d1 <- (log(S/K) + (sigma^2/2) * TTM)/(sigma * sqrt(TTM))
  d2 <- d1 - sigma * sqrt(TTM)
  if(call) price <- pnorm(d1) * S - K * pnorm(d2)
  else price <- K * pnorm(-d2) - pnorm(-d1) * S 
  return(price)
}

BSMVega <- function(S,K,sigma,TTM){
  d1 <- (log(S/K) + (sigma^2/2) * TTM)/(sigma * sqrt(TTM))
  return(S * dnorm(d1) * sqrt(TTM))
}

impVol <- function(price,S,K,TTM){
  f <- function(sigma) price - BSM(S,K,sigma,TTM)
  impVol <- uniroot(f,c(-1,1))
  return(impVol$root)
}

library(tidyverse)

data <- read_csv2('/Users/tk/Documents/GitHub/Speciale/data/calibration_data.csv')

data <- data[,c("Strike", "Price", "Expiry", "impVol")]

#maturities <- c(10,17,30,66)

#data <- data[data$Expiry %in% maturities,]

l <- length(data$Price); S0 <- 4240

for(i in 1:l){
  data[i,'Vega'] <- BSMVega(S0,data$Strike[i],data$impVol[i],data$Expiry[i]/365) / 100
}

#par <- optim(par = c(0.05, 0.3, 0.5, 0.05, -0.7, 0.1, -0.1, 0.1), lossFunction, method = "L-BFGS-B", lower = c(0.01, 0.01, 0.01, 0.01, -0.99, 0.01, -0.99, 0.01), upper = c(10, 10, 10, 10, 0.99, 10, 0.99, 10))

v <- 0.06809358
sigma <- 0.65239094
kappa <- 0.40659087
theta <- 0.09848509
rho <- -0.84379886
lambda <- 0.14531402
mu <- -0.21302818
delta <- 0.24271881

for(i in 1:l){
  BatesPrice <- BatesCall(S0,data$Strike[i],data$Expiry[i]/365,v,sigma,kappa,theta,rho,lambda,mu,delta)
  if(BatesPrice<0){ 
    next
  }
  else{
  data[i,'BatesImp'] <- impVol(BatesPrice,S0,data$Strike[i],data$Expiry[i]/365)
  }
}

data <- data[!(is.na(data$BatesImp)),]

impVolPlot <- ggplot(data = data , mapping = aes(log(Strike/S0), impVol)) + 
  geom_line(aes(x = log(Strike/S0), y = impVol, colour = factor(Expiry))) + 
  geom_line(aes(x = log(Strike/S0), y = BatesImp, colour = factor(Expiry)), linetype="dashed") +
  xlab('log(K/S)') + 
  ylab('Implied volatility') + 
  ggtitle('Calibration result in the Bates model') +
  labs(colour = "Expiry"); impVolPlot

data$error <- abs(data$impVol-data$BatesImp)
mean(data$error)
median(data$error)
sd(data$error)
hist(data$error,col='blue', xlab='Absolute error', main='Absolute errors in the Bates calibration')

