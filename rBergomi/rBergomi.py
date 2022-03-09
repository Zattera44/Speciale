import numpy as np   
import scipy as sp 

### This class is inspired by the repo found @ https://github.com/ryanmccrickerd/rough_bergomi,
###                as well as the paper found @ https://arxiv.org/abs/1708.02563.

### The method 'sim_Volterra()' utilises results derived in https://arxiv.org/abs/1507.03004.

class rBergomi():
    def __init__(self, npaths: int, nsteps: int, T: float, \
                    alpha: float, eta: float, rho: float):
        self.npaths = npaths
        self.nsteps = nsteps
        self.T = T
        self.dt = 1.0/nsteps
        self.s = int(nsteps * T)
        self.t = np.linspace(0, T, 1 + self.s)[np.newaxis,:]
        self.alpha = alpha
        self.eta = eta
        self.rho = rho
        self.dW = self.__sim_dW()
        self.dZ = self.__sim_dZ()

    def __cov(self):
        cov = np.array([[0.,0.],[0.,0.]])
        cov[0,0] = 1/self.nsteps
        cov[0,1] =  cov[1,0] = 1./((1.*self.alpha+1) * self.nsteps**(1.*self.alpha+1))
        cov[1,1] = 1./((2.*self.alpha+1) * self.nsteps**(2.*self.alpha+1))
        return cov

    def __sim_dW(self):
        return np.random.multivariate_normal([0,0], self.__cov(), (self.npaths, self.s))     

    def __sim_dZ(self):
        return np.random.randn(self.npaths, self.s) * np.sqrt(self.dt)

    def __sim_dB(self):
        return self.rho * self.dW[:,:,0] + np.sqrt(1-self.rho**2) * self.dZ

    def __g(self, x):
        return x**self.alpha

    def __b(self, k):
        return ((k**(self.alpha+1)-(k-1)**(self.alpha+1))/(self.alpha+1))**(1/self.alpha)

    def sim_volterra_hybrid(self):
        dW = self.__sim_dW()
        X_caret = np.zeros((self.npaths, self.s+1))

        for i in np.arange(1, 1 + self.s, 1):
            X_caret[:,i] = dW[:,i-1,1]

        convolve_term = np.zeros(self.s+1)
        for k in np.arange(2,self.s+1,1):
            x = self.__b(k)/self.nsteps
            convolve_term[k] = self.__g(x)

        X_hat = np.zeros((self.npaths, len(convolve_term)+len(dW[0,:,0])-1))
        for j in np.arange(0,self.npaths,1):
            X_hat[j,:] = sp.signal.convolve(convolve_term,dW[j,:,0],mode='full',method='auto')
        X_hat = X_hat[:,:self.s+1]

        return np.sqrt(2*self.alpha+1) * (X_caret + X_hat)

    def sim_Gatheral(self):
        pass

    def sim_V(self, xi = 1.0, eta = 1.0, hybrid = True):
        if hybrid:
            volterra = self.sim_Volterra() 
            return xi*np.exp(eta*volterra-0.5*eta**2*self.t**(2*self.alpha+1))
        else:
            dW = self.dW[:,:,0]
            increments = eta * np.sqrt(2*self.alpha+1) * (self.t - u)**self.alpha * dW

    def sim_S(self, S0 = 100): #WIP
        dB = self.sim_dB()
        V = self.sim_V()
        increments = np.sqrt(V[:,:-1]) * dB - 0.5 * V[:,:-1] * self.dt
        integral = np.cumsum(increments, axis = 1)

    def xi(self, heston = True, theta =  0.035, kappa = 1, v = 0.04):
        if heston:
            xi = theta + (v - theta) * np.exp(-kappa*self.t)
        return xi