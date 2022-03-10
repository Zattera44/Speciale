import numpy as np   
import scipy.special
import scipy.signal
from itertools import product

### This class is inspired by the repo found @ https://github.com/ryanmccrickerd/rough_bergomi,
###                as well as the paper found @ https://arxiv.org/abs/1708.02563.

### The method 'sim_Volterra()' utilises results derived in https://arxiv.org/abs/1507.03004.

class rBergomi():
    def __init__(self, npaths = 10000, nsteps = 52, T = 1, \
                    H = -0.1, eta = 1.9, rho = -0.3, alpha = -0.43, \
                    hybrid = True):
        self.npaths = npaths
        self.nsteps = nsteps
        self.T = T
        self.dt = 1.0/nsteps
        self.s = int(nsteps * T)
        self.t_hybrid = np.linspace(0, T, 1 + self.s)[np.newaxis,:]
        self.t = np.linspace(0, T, self.s)
        self.H = H
        #self.alpha = H - 0.5
        self.alpha = alpha
        self.eta = eta
        self.rho = rho
        self.hybrid = hybrid

        if hybrid:
            self.dW = self.__sim_dW()
            self.dZ = self.__sim_dZ()
        else:
            self.normal = self.__normal()

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

    def __G(self, x):
        gamma = 0.5 - self.H
        F = scipy.special.hyp2f1(1.0,gamma,2.0-gamma,1.0/x)
        return (1.0-2.0*gamma)/(1.0-gamma)*x**(-gamma)*F

    def __check_symmetric(self, a):
        return np.allclose(a, a.T)

    def WZ(self):
        idx = self.t
        v, u = np.zeros(pow(self.s,2)), np.zeros(pow(self.s,2))

        k = 0
        for i in product(idx, idx):
            v[k], u[k] = i
            k = k + 1

        D = np.sqrt(2.0*self.H)/(self.H+0.5)
        res = self.rho * D * (v**(self.H+0.5)-(v-np.minimum(u,v))**(self.H+0.5))
        return np.reshape(res, (self.s,self.s))

    def WW(self):
        cov = np.zeros((self.s,self.s))

        i = 0
        for v in self.t:
            j = 0
            for u in self.t:
                if i == j:
                    cov[i,j] = u**(2.0*self.H)
                if i > j:
                    x = v/u
                    cov[i,j] = u**(2.0*self.H)*self.__G(x)
                if i < j:
                    cov[i,j] = 0
                j += 1
            i += 1

        return cov + cov.T - np.diag(cov.diagonal())

    def ZZ(self):
        cov = np.zeros((self.s,self.s))

        i = 0
        for v in self.t:
            j = 0
            for u in self.t:
                cov[i,j] = np.min([v,u])
                j += 1
            i += 1

        return cov

    def joint_cov(self):
        res = np.block([[self.WW(), self.WZ()],[self.WZ().T, self.ZZ()]])
        if self.__check_symmetric(res):
            return res
        else:
            print('Covariance matrix is not symmetric...')

    def paths(self):
        normal = np.random.rand(self.s*2,self.npaths)
        cov = self.joint_cov()
        sigma = scipy.linalg.cholesky(cov,lower=True)
        control = np.all(np.matmul(sigma,sigma.T) - cov < 1e-06)
        if control:
            print("Good news, everyone! I think I perfected a scheme that will simulate all paths!")
            return np.matmul(sigma, normal)
        else: 
            print("Bad news, everyone! I don't think the simulation is going to make it...")

    def sim_Gatheral(self):
        pass

    def sim_volterra(self):
        if self.hybrid:
            dW = self.__sim_dW()
            X_caret = np.zeros((self.npaths,self.s+1))
            for i in np.arange(1,1+self.s,1):
                X_caret[:,i] = dW[:,i-1,1]
            convolve_term = np.zeros(self.s+1)
            for k in np.arange(2,self.s+1,1):
                x = self.__b(k)/self.nsteps
                convolve_term[k] = self.__g(x)

            X_hat = np.zeros((self.npaths,len(convolve_term)+len(dW[0,:,0])-1))
            for j in np.arange(0,self.npaths,1):
                X_hat[j,:] = scipy.signal.convolve(convolve_term,dW[j,:,0],mode='full',method='auto')
                #X_hat[j,:] = np.convolve(convolve_term,dW[j,:,0])
            X_hat = X_hat[:,:self.s+1]

            return np.sqrt(2 * self.alpha + 1) * (X_caret + X_hat)

    def sim_V(self, xi = 0.235**2):
        volterra = self.sim_volterra() 
        if self.hybrid:
            return xi * np.exp(self.eta * volterra-0.5 * self.eta**2 * self.t_hybrid**(2 * self.alpha + 1))

    def xi(self, heston = True, theta =  0.035, kappa = 1, v = 0.04):
        if heston:
            xi = theta + (v - theta) * np.exp(-kappa*self.t_hybrid)
        return xi