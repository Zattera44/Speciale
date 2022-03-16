simulateVol <- function(nPaths,nDates,ttm,v,sigma,kappa,theta,rho){
  
  dt <- ttm/nDates
  Vmat <- matrix(nrow = nPaths, ncol = nDates)
  V <- rep(v, nPaths)
  
  for(i in 1:nDates){
    
    z <- rnorm(nPaths); dW <- sqrt(dt)*rnorm(nPaths)
    dZ <- rho*dW+sqrt(1-rho^2)*sqrt(dt)*z

    V[1:nPaths]<-V[1:nPaths]+kappa*(theta-pmax(V[1:nPaths],0))*dt+sigma*sqrt(pmax(V[1:nPaths],0))*dZ

    Vmat[,i] <- V
  }
  
  result <- Vmat
  
}

nPaths <- 1000; nDates <- 365; ttm <- 1; sigma <- 0.1; rho <- -0.1
v <- 0.01; kappa <- 1; theta <- 0.1

sim1 <- simulateVol(nPaths,nDates,ttm,v,sigma,kappa,theta,rho)

v <- 0.01; kappa <- 2

sim2 <- simulateVol(nPaths,nDates,ttm,v,sigma,kappa,theta,rho)

v <- 0.2; kappa <- 1

sim3 <- simulateVol(nPaths,nDates,ttm,v,sigma,kappa,theta,rho)


v <- 0.2; kappa <- 2

sim4 <- simulateVol(nPaths,nDates,ttm,v,sigma,kappa,theta,rho)

kLvL <- rowMeans(t(sim1))
kHvL <- rowMeans(t(sim2))
kLvH <- rowMeans(t(sim3))
kHvH <- rowMeans(t(sim4))

jpeg('kappa_theta.jpg', width = 1920, height = 1080, res = 300)
plot(kLvL, ylim=c(0,0.21), 
     type='l', lty=2, lwd=3, col='blue', 
     ylab = 'Pathwise average volatility', xlab='Day',
     main = 'Mean reversion in the Heston model')
lines(kHvL, col='blue', lwd=3)
lines(kLvH, lty=2, lwd=3, col='pink')
lines(kHvH, col='pink', lwd=3)
abline(0.1, 0, col = 'grey')
dev.off()



