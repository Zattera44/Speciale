charFunc <- function(u, spot, ttm, v, sigma, kappa, theta, rho){
  d <- sqrt( (rho*sigma*u*1i-kappa)^2 + sigma^2*(u*1i+u^2) )
  d <- -d
  g <- (kappa - rho*sigma*u*1i + d) / (kappa - rho*sigma*u*1i - d)
  tempM <- (kappa-rho*sigma*u*1i+d)*ttm-2*log((1-g*exp(d*ttm))/(1-g))
  M <- (kappa*theta)/(sigma^2)*tempM
  N <- (kappa-rho*sigma*u*1i+d)/(sigma^2)*((1-exp(d*ttm))/(1-g*exp(d*ttm)))
  res <- exp(M+N*v+1i*u*log(spot))
  return(res)
}

HestonCall <- function(spot, strike, ttm, v, sigma, kappa, theta, rho){
  integrand1 <- function(u){ 
    num1 <- charFunc(u-1i, spot, ttm, v, sigma, kappa, theta, rho)
    den1 <- charFunc(-1i, spot, ttm, v, sigma, kappa, theta, rho)
    dummy1 <- exp(-1i*u*log(strike))*num1/(1i*u*den1)
    integrand1 <- Re(dummy1)
  }
  integrand2 <- function(u){ 
    dummy2 <- exp(-1i*u*log(strike))*charFunc(u, spot, ttm, v, sigma, kappa, theta, rho)/(1i*u)
    integtand2 <- Re(dummy2)
  }
  Pi1 <- 0.5 + 1/pi * integrate(integrand1,0,Inf,stop.on.error = FALSE)$value
  Pi2 <- 0.5 + 1/pi * integrate(integrand2,0,Inf,stop.on.error = FALSE)$value
  res <- spot*Pi1 - strike*Pi2
  return(res)
}

BSMVega <- function(S,K,sigma,TTM){
  d1 <- (log(S/K) + (sigma^2/2) * TTM)/(sigma * sqrt(TTM))
  return(S * dnorm(d1) * sqrt(TTM))
}

lossFunction <- function(parms){
  sum <- 0
  for(i in 1:l){
    sum <- sum + ( (data$Price[i] - HestonCall(S0,data$Strike[i],data$Expiry[i]/365,parms[1],parms[2],parms[3],parms[4],parms[5])) / (data$Vega[i] * 100) )^2
  }
  return(sum)
}

BSM <- function(S,K,sigma,TTM){
  d1 <- (log(S/K) + (sigma^2/2) * TTM)/(sigma * sqrt(TTM))
  d2 <- d1 - sigma * sqrt(TTM)
  price <- pnorm(d1) * S - K * pnorm(d2)
  return(price)
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
  data[i, 'Vega'] <- BSMVega(S0,data$Strike[i],data$impVol[i],data$Expiry[i]/365) / 100
}

#par <- optim(par = c(0.1^2, 0.8, 2.01, 0.1^2, -0.75), lossFunction, method = "L-BFGS-B", lower = c(0.01,0.01,0.1,0.01,-0.99), upper = c(0.99,2,5,0.99,0.99))

v <- 0.08118792
sigma <- 0.78055835
kappa <- 2.12469802
theta <- 0.05161678
rho <- -0.74892488

for(i in 1:l){
  HestonPrice <- HestonCall(S0,data$Strike[i],data$Expiry[i]/365,v,sigma,kappa,theta,rho)
  if(HestonPrice<0){ 
    next
  }
  else{ 
    data[i,'HestonImp'] <- impVol(HestonPrice,S0,data$Strike[i],data$Expiry[i]/365)
  }
}

data <- data[!(is.na(data$HestonImp)),]

impVolPlot <- ggplot(data = data , mapping = aes(log(Strike/S0), impVol)) + 
  geom_line(aes(x = log(Strike/S0), y = impVol, colour = factor(Expiry))) + 
  geom_line(aes(x = log(Strike/S0), y = HestonImp, colour = factor(Expiry)), linetype="dashed") +
  xlab('log(K/S)') + 
  ylab('Implied volatility') + 
  ggtitle('Calibration result in the Heston model') +
  labs(colour = "Expiry"); impVolPlot

data$error <- abs(data$impVol - data$HestonImp)
mean(data$error)
median(data$error)
sd(data$error)
hist(data$error, col='blue', xlab='Absolute error', main='Absolute errors in the Heston calibration')



