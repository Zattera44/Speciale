import numpy as np    

class rBergomi():
    def __init__(self, npaths: int, nsteps: int, T: float, \
                    H: float, eta: float, rho: float):

        self.npaths = npaths
        self.nsteps = nsteps
        self.T = T
        self.H = H
        self.eta = eta
        self.rho = rho
        self.dt = 1.0/nsteps
        self.s = int(nsteps * T)
        self.t = np.linspace(0, T, 1 + self.s)[np.newaxis,:]
        self.alpha = 0.5*(2.0*H-1.0)

    def xi(heston = True):
        pass

    def __cov(self):
        cov = np.array([0,0],[0,0])
        cov[0,0] = 1/self.nsteps
        cov[0,1] = cov[1,0] = 1./((1.*self.a+1) * self.nsteps**(1.*self.a+1))
        cov[1,1] = 1./((2.*self.a+1) * self.nsteps**(2.*self.a+1))
        return cov

    def dW(self):
        return self.__cov()
        #return np.random.multivariate_normal([0,0], self.__cov, (self.npaths, self.s))

    def dZ(self):
        return np.random.randn(self.npaths, self.s) * np.sqrt(self.dt)

    def dB(self):
        pass