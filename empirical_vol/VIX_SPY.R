library(moments)

spy <- read.csv2('data/SPY.csv', sep = ',')
longspy <- read.csv2('data/longSPY.csv', sep = ',')
vix <- read.csv2('data/VIX.csv', sep = ',')
longvix <- read.csv2('data/longVIX.csv', sep = ',')

spy <- spy[c('Date','Adj.Close')]
longspy <- longspy[c('Date','Adj.Close')]
vix <- vix[c('Date','Adj.Close')]
longvix <- longvix[c('Date','Adj.Close')]


spy$Adj.Close <- as.numeric(spy$Adj.Close)
longspy$Adj.Close <- as.numeric(longspy$Adj.Close)
vix$Adj.Close <- as.numeric(vix$Adj.Close)
longvix$Adj.Close <- as.numeric(longvix$Adj.Close)

cor(spy$Adj.Close, vix$Adj.Close)

change <- diff(vix$Adj.Close)

cor(change[-length(change)], change[-1])

#jpeg('spy_vix.jpg', width = 1920, height = 1080, res = 300)
plot(as.Date(vix$Date), vix$Adj.Close/max(vix$Adj.Close), type='l', col ='green',
     xlab='', ylab='')
lines(as.Date(vix$Date), spy$Adj.Close/max(spy$Adj.Close), col = 'blue')
#dev.off()

logSpy <- diff(log(longspy$Adj.Close))
sim <- rnorm(n = length(logSpy), mean = mean(logSpy), sd = sd(logSpy))

#jpeg('spy_hist.jpg', width = 1920, height = 1080, res = 300)
hist(logSpy, breaks = 100, col = 'blue', freq = FALSE, main = '', xlab = 'Log returns')
lines(density(sim), col = 'red', lwd = 2)
#dev.off()

skewness(logSpy); skewness(sim)
kurtosis(logSpy); kurtosis(sim)

longdiff <- diff(log(longvix$Adj.Close))
sim2 <- rnorm(n = length(longdiff), mean = mean(longdiff), sd = sd(longdiff))

jpeg('vix_hist.jpg', width = 1920, height = 1080, res = 300)
hist(longdiff, breaks = 100, freq = FALSE, col = 'green',
     xlab = 'ln(VIX_t) - ln(VIX_t-1)', main = '')
lines(density(sim2), col = 'red', lwd = 2)
dev.off()

#jpeg('vix_diff.jpg', width = 1920, height = 1080, res = 300)
plot(longdiff~as.Date(longvix$Date)[-1], type = 'l', col = 'green',
     xlab = '', ylab='ln(VIX_t) - ln(VIX_t-1)')
#dev.off()

skewness(longdiff); skewness(sim2)
kurtosis(longdiff); kurtosis(sim2)
