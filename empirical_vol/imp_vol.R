library(quantmod)
library(tidyverse)

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

dates <- c('2022-03-21',
           '2023-03-17')

SPX <- getOptionChain('^SPX', NULL)

expiry_dates <- as.Date(names(SPX), format = '%b.%d.%Y')
expiry_days <- expiry_dates - Sys.Date()

i = 1
for(df in SPX){
    #if(as.character(as.Date(names(SPX)[i], format = '%b.%d.%Y')) %in% dates){
      if(i == 1){
        merged <- df$call
        merged$Expiry <- as.numeric(expiry_days[i])
      }
      else{
        if(is.null(df$calls)) next
        current <- df$calls
        current$Expiry <- as.numeric(expiry_days[i])
        merged <- rbind(merged, current)
      }
    #}
  i = i + 1
}

data <- merged

S0 <- getQuote('^SPX')$Last; r <- 0

data$Price <- (data$Bid+data$Ask)/2
data$Moneyness <- log(data$Strike/S0)
data <- data[abs(data$Moneyness) < 0.7,]
data <- data[data$IV > mean(data$IV)-sd(data$IV),]
data <- data[data$IV<mean(data$IV)+sd(data$IV),]
data <- data[c('Strike','Price', 'Expiry', 'IV', 'Moneyness')]
data <- data[!(is.na(data$Price)),]
l <- length(data$Price)

for(i in 1:l){
  data[i, 'impVol'] <- impVol(data$Price[i], data$Strike[i], data$Expiry[i]/365, S0, r)
  print(i)
}  

data <- data[!(is.na(data$impVol)),]

nits <- 0
errors = TRUE
while(errors){
  print(nits)
  len <- 0
  for(expiry in unique(data$Expiry)){
    k = 1
    remove <- vector()

    for(i in 2:(length(data[data$Expiry==expiry, 'impVol'])-1)){

      if((data[data$Expiry==expiry, 'impVol'][i-1] < data[data$Expiry==expiry, 'impVol'][i]) && (data[data$Expiry==expiry, 'impVol'][i+1] < data[data$Expiry==expiry, 'impVol'][i])){
        remove[k] <- i
        k <- k + 1
      }
      len <- len + length(remove)
    }
    data[data$Expiry==expiry, 'impVol'][remove] <- NA
  }
  if(len == 0){
    errors <- FALSE
  }
  data <- data[!(is.na(data$impVol)),]
  nits <- nits + 1
}

maturities <- unique(data$Expiry)

for(mat in maturities){
plot <- ggplot(data = data[data$Expiry==mat,] , mapping = aes(Moneyness, impVol)) + 
  geom_line(aes(colour = factor(Expiry))) + xlab('log(K/S)') + ylab('Implied volatility') +
  labs(colour = "Expiry")
print(plot)
}

filtered <- c(6,10,16,17,22,24,30,38,52,66,94,136,157,169,199,220,248,311,339,367,458)

for(mat in filtered){
  plot <- ggplot(data = data[data$Expiry==mat,] , mapping = aes(Moneyness, impVol)) + 
    geom_line(aes(colour = factor(Expiry))) + xlab('log(K/S)') + ylab('Implied volatility') +
    labs(colour = "Expiry")
  print(plot)
}

copy <- copy[copy$Expiry %in% filtered,]

plot <- ggplot(data = copy , mapping = aes(Moneyness, impVol)) + 
  geom_line(aes(colour = factor(Expiry))) + xlab('log(K/S)') + ylab('Implied volatility') +
  labs(colour = "Expiry"); plot

test <- c(10,17,30,66,94,136,220,339)

copy2 <- copy[copy$Expiry %in% test,]

plot <- ggplot(data = copy2 , mapping = aes(Moneyness, impVol)) + 
  geom_line(aes(colour = factor(Expiry))) + xlab('log(K/S)') + ylab('Implied volatility') +
  labs(colour = "Expiry"); plot





