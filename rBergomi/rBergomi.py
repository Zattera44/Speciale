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

        # Discretize
        self.nsteps = nsteps # Number of steps per year
        self.T = T # Time to expiry
        self.dt = 1.0/nsteps # Stepsize
        self.s = int(nsteps * T) # Number of steps in total

        # The grid changes slightly depending on which simulation scheme is used
        if hybrid:
            self.t = np.linspace(0, T, 1 + self.s)
        else:
            self.t = np.linspace(T/self.s,T,self.s)

       
        self.H = H # Hurst index

        # To be consistent with notation in paper
        if alpha is not None:
            self.alpha = alpha
        else:
            self.alpha = H - 0.5

        self.eta = eta # "Vol of vol"
        self.rho = rho # Correlation between vol and asset
        self.hybrid = hybrid # Boolean

        # Get paths used in hybrid scheme
        if hybrid:
            self.dW = self.__sim_dW()
            self.dZ = self.__sim_dZ()

        # Simulate the Volterra paths used to simulate V   
        self.paths = self.__sim_paths() 
        self.volterra = self.__sim_volterra() 

    # Covariance matrix used in hybrid scheme
    def __cov(self):
        cov = np.array([[0.,0.],[0.,0.]])
        cov[0,0] = 1/self.nsteps
        cov[0,1] =  cov[1,0] = 1./((1.*self.alpha+1) * self.nsteps**(1.*self.alpha+1))
        cov[1,1] = 1./((2.*self.alpha+1) * self.nsteps**(2.*self.alpha+1))
        return cov

    # Paths related to V in hybrid scheme
    def __sim_dW(self):
        return np.random.multivariate_normal([0,0], self.__cov(), (self.npaths, self.s))     

    # Paths related to S in hybrid scheme
    def __sim_dZ(self):
        return np.random.randn(self.npaths, self.s) * np.sqrt(self.dt)

    # Based on above path, create paths with correct correlation
    def __sim_dB(self):
        return self.rho * self.dW[:,:,0] + np.sqrt(1-self.rho**2) * self.dZ

    # Used to simulate Volterra paths in hybrid scheme
    def __g(self, x):
        return x**self.alpha

    # Used to simulate Volterra paths in hybrid scheme
    def __b(self, k):
        return ((k**(self.alpha+1)-(k-1)**(self.alpha+1))/(self.alpha+1))**(1/self.alpha)

    # Simulate paths through the method proposed by Gatheral
    # Only used if hybrid = False
    def __sim_paths(self):

        v, u = np.zeros(pow(self.s,2)), np.zeros(pow(self.s,2))

        k = 0
        for i in product(self.t, self.t):
            v[k], u[k] = i
            k = k + 1

        # Matrix keeping track of which minimum time index in each entry
        index_matrix = np.reshape(np.minimum(v,u),(self.s,self.s))


        # Used to find E[WW], see Gatheral
        def G(H, x): 
            F = scipy.special.hyp2f1(1.0,0.5-H,1.5+H,1.0/x)
            return (2*H*(x**(H-0.5))*scipy.special.gamma(0.5+H)*F)/scipy.special.gamma(H+1.5)

        # Initialize E[WW] and E[WZ]
        WW = np.zeros((self.s,self.s))
        WZ = np.zeros((self.s,self.s))

        # Correlation of Browninan motion with itself
        ZZ = index_matrix

        # Used to find E[WZ], see Gatheral
        D = np.sqrt(2.0*self.H)/(self.H+0.5)

        # Construction E[WW] and E[WZ]
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

        # Utilizing that E[WW] is symmetric 
        WW = WW + WW.T - np.diag(WW.diagonal())

        # Creating covariance matrix for (W,Z)~N(0,Sigma)
        cov = np.block([[WW, WZ],[WZ.T, ZZ]]) # = Sigma

        mean = np.zeros(self.s*2)

        # Simulating the specified number of paths by drawing 
        # from the joint distribution of W and Z
        paths = np.random.multivariate_normal(mean=mean,cov=cov, size=self.npaths)

        ### Slower, but maybe more illustrating way of achieving the above ###

        #A = scipy.linalg.cholesky(cov, lower=True)
        #normal = np.random.randn(steps*2,npaths)
        #paths = np.matmul(A, normal)

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
            paths[:,1:] = self.paths[:,:self.s]
            return paths

    def sim_V(self, xi = 0.235**2):
        t = np.linspace(0,self.T,self.s+1)
        return xi * np.exp(self.eta * self.volterra-0.5 * self.eta**2 * t**(2 * self.alpha + 1))

    def sim_S(self, S0 = 100):

        V = self.sim_V()
        if self.hybrid:
            dB = self.__sim_dB()
        else:
            #dB  = np.zeros((self.npaths,self.s+1))
            #dB[:,0] = 0
            #dB[:,1:] = self.paths[:,self.s:]
            dB = self.paths[:,self.s:]

        increments = np.sqrt(V[:,:-1]) * dB - 0.5 * V[:,:-1] * self.dt

        integral = np.cumsum(increments, axis = 1)

        S = np.zeros_like(V)
        S[:,0] = S0
        S[:,1:] = S0 * np.exp(integral)
        return integral

    def xi(self, heston = True, theta =  0.035, kappa = 1, v = 0.04):
        if heston:
            xi = theta + (v - theta) * np.exp(-kappa*self.t)
        return xi