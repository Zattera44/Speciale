library(tidyverse)

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

BSMVega <- function(S,K,r,sigma,TTM){
  d1 <- (log(S/K) + (r + sigma^2/2) * TTM)/(sigma * sqrt(TTM))
  return(S * dnorm(d1) * sqrt(TTM))
}

lossFunction <- function(parms){
  sum <- 0
  for(i in 1:l){
    sum <- sum + ( (data[[i,2]] - HestonCall(S0,data$Strike[i],0,0,data$Expiry[i],parms[1],parms[2],parms[3],parms[4],parms[5])) / (data$Vega[i] * 100) )^2
  }
  return(sum)
}

BSM <- function(S,K,r,sigma,TTM){
  d1 <- (log(S/K) + (r + sigma^2/2) * TTM)/(sigma * sqrt(TTM))
  d2 <- d1 - sigma * sqrt(TTM)
  price <- pnorm(d1) * S - exp(-r * TTM) * K * pnorm(d2)
  return(price)
}

impVol <- function(price,K,TTM,S,r,increment=0.00001,limit=1000){
  limit <- price/limit
  continue <- TRUE
  sigma <- increment
  while(continue){
    eps <- price - BSM(S,K,r,sigma,TTM)
    if(abs(eps) < limit){
      continue <- FALSE
    }
    sigma <- sigma + increment
    if(sigma > 2){
      continue <- FALSE
      sigma <- NA
    }
  }
  return(sigma)
}

for(i in 1:l){
  data[i, 'vega'] <- BSMVega(S0,data$Strike[i],r,data$impVol[i],data$Expiry[i]/365) / 100
}

par <- optim(par = c(0.2^2, 1, 2, 0.2^2, -0.25), lossFunction, method = "L-BFGS-B",
             lower = c(0.01,0.01,0.01,0.01,-1), upper = c(0.99,10,10,0.99,1))


v <- 0.08525995
sigma <- 2.99525485
kappa <- 2.19443040
theta <- 0.11726362
rho <- -0.84045985


for(i in 1:l){
  HestonPrice <- HestonCall(S0,data[[i,2]],0,0,data[[i,3]],v,sigma,kappa,theta,rho)
  data[i,'HestonImp'] <- impVol(HestonPrice,S0,data[[i,2]],r,data[[i,3]])
}

impVolPlot <- ggplot(data = data , mapping = aes(log(Strike/S0), impliedVolatility)) + 
  geom_line(aes(x = log(Strike/S0), y = impliedVolatility, colour = factor(Expiry))) + 
  geom_line(aes(x = log(Strike/S0), y = HestonImp, colour = factor(Expiry)), linetype="dashed") +
  xlab('log(K/S)') + 
  ylab('Implied volatility') + 
  labs(colour = "Expiry"); impVolPlot
