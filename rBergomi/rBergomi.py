import numpy as np   
import scipy.special
import scipy.signal
from itertools import product

### This class is inspired by the repo found @ https://github.com/ryanmccrickerd/rough_bergomi,
###                as well as the paper found @ https://arxiv.org/abs/1708.02563.

### The method 'sim_volterra()' utilises results derived in https://arxiv.org/abs/1507.03004.
### The simulation scheme when hybrid = False is the one found in ...

class rBergomi():
    def __init__(self, npaths = 10000, nsteps = 52, T = 1, \
                    H = 0.1, eta = 1.9, rho = -0.3, alpha = None, \
                    hybrid = True):

        self.npaths = npaths
        self.nsteps = nsteps
        self.T = T

        self.dt = 1.0/nsteps
        self.s = int(nsteps * T)

        if hybrid:
            self.t = np.linspace(0, T, 1 + self.s)#[np.newaxis,:]
        else:
            self.t = np.linspace(T/self.s,T,self.s)

        self.H = H

        if alpha is not None:
            self.alpha = alpha
        else:
            self.alpha = H - 0.5

        self.eta = eta
        self.rho = rho
        self.hybrid = hybrid

        if hybrid:
            self.dW = self.__sim_dW()
            self.dZ = self.__sim_dZ()

        self.volterra = self.__sim_volterra() 

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

    def __sim_paths(self):

        v, u = np.zeros(pow(self.s,2)), np.zeros(pow(self.s,2))

        k = 0
        for i in product(self.t, self.t):
            v[k], u[k] = i
            k = k + 1

        index_matrix = np.reshape(np.minimum(v,u),(self.s,self.s))

        def G(H, x):
            F = scipy.special.hyp2f1(1.0,0.5-H,1.5+H,1.0/x)
            return (2*H*(x**(H-0.5))*scipy.special.gamma(0.5+H)*F)/scipy.special.gamma(H+1.5)

        WW = np.zeros((self.s,self.s))
        WZ = np.zeros((self.s,self.s))

        ZZ = index_matrix

        D = np.sqrt(2.0*self.H)/(self.H+0.5)

        i = 0
        for v in self.t:
            j = 0
            for u in self.t:
                if i == j:
                    WW[i,j] = u**(2.0*self.H)
                if i > j:
                    x = v/u
                    WW[i,j] = u**(2.0*self.H)*G(self.H,x)
                if i < j:
                    WW[i,j] = 0
                WZ[i,j] = D * (v**(self.H+0.5)-(v-np.minimum(u,v))**(self.H+0.5))
                j += 1
            i += 1

        WW = WW + WW.T - np.diag(WW.diagonal())

        cov = np.block([[WW, WZ],[WZ.T, ZZ]])

        mean = np.zeros(self.s*2)

        #A = scipy.linalg.cholesky(cov, lower=True)
        #normal = np.random.randn(steps*2,npaths)
        #paths = np.matmul(A, normal)

        paths = np.random.multivariate_normal(mean=mean,cov=cov, size=self.npaths)

        return paths

    def __sim_volterra(self):
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
                #X_hat[j,:] = scipy.signal.convolve(convolve_term,dW[j,:,0],mode='full',method='auto')
                X_hat[j,:] = np.convolve(convolve_term,dW[j,:,0])
            X_hat = X_hat[:,:self.s+1]

            return np.sqrt(2 * self.alpha + 1) * (X_caret + X_hat)
        
        else:
            paths = np.zeros((self.npaths,self.s+1))
            paths[:,0] = 0
            paths[:,1:] = self.__sim_paths()[:,:self.s]
            return paths

    def sim_V(self, xi = 0.235**2):
        t = np.linspace(0,self.T,self.s+1)
        return xi * np.exp(self.eta * self.volterra-0.5 * self.eta**2 * t**(2 * self.alpha + 1))

    def xi(self, heston = True, theta =  0.035, kappa = 1, v = 0.04):
        if heston:
            xi = theta + (v - theta) * np.exp(-kappa*self.t)
        return xi