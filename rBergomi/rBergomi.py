import numpy as np
import scipy.special

import jax.numpy as jnp
import jax.scipy as jsp
from jax import grad, value_and_grad

from itertools import product

class roughBergomi:
    def __init__(self, n_paths = 100000, n_steps = 52,
                    T = 1, H = 0.1, eta = 1, rho = -.3,
                    xi = 0.235**2):

        # number of paths
        self.n_paths = n_paths
        # number of steps per year
        self.n_steps = n_steps
        # time horizon, in years
        self.T = T
        # size of time step
        self.dt = 1.0/n_steps
        # number of total steps 
        self.s = int(n_steps * T)
        # time grid 
        self.t = np.linspace(T/self.s,T,self.s)

        # parameters of rough Bergomi model
        self.H = H
        self.eta = eta 
        self.rho = rho
        self.xi = xi

    def G(self, x): 
        F = scipy.special.hyp2f1(1.0,0.5-self.H,1.5+self.H,1.0/x)
        return (2*self.H*(x**(self.H-0.5))*scipy.special.gamma(0.5+self.H)*F)/scipy.special.gamma(self.H+1.5)

    def idx_matrix(self):
        v, u = np.zeros(pow(self.s, 2)), np.zeros(pow(self.s, 2))
        k = 0
        for i in product(self.t, self.t):
            v[k], u[k] = i
            k = k + 1
        return np.reshape(np.minimum(v, u), (self.s, self.s))

    def covariance_matrix(self):
        ZZ = self.idx_matrix
        WW = np.zeros((self.s, self.s))
        WZ = np.zeros((self.s, self.s))
        D = np.sqrt(2.0*self.H)/(self.H+0.5)
        i = 0
        for v in self.t:
            j = 0
            for u in self.t:
                if i == j:
                    WW[i,j] = u**(2.0*self.H)
                if i > j:
                    x = v/u
                    WW[i,j] = u**(2.0*self.H)*self.G(x)
                if i < j:
                    WW[i,j] = 0
                WZ[i,j] = D * (v**(self.H+0.5)-(v-np.minimum(u,v))**(self.H+0.5))
                j += 1
            i += 1
        WW = WW + WW.T - np.diag(WW.diagonal())
        return np.block([[WW, WZ],[WZ.T, ZZ]])

    def simulate_paths(self):
        cov = self.covariance_matrix
        mean = np.zeros(self.s*2)
        paths = np.random.multivariate_normal(mean = mean, cov = cov, size = self.n_paths)
        volterra, W = np.zeros((self.n_paths , self.s+1)), np.zeros((self.n_paths, s+1)) 
        volterra[:,0] = 0; volterra[:,1:] = paths[:,:self.s]
        W[:,0] = 0; W[:,1:] = paths[:,self.s:]
        dW = np.diff(W, axis = 1)
        return volterra, W, dW

    def simulate_V(self, volterra):
        t = np.linspace(0, self.T, self.s+1)
        return self.xi * np.exp(self.eta * volterra - 0.5 * self.eta**2 * t**(2 * (self.H - 0.5) + 1))

    def simulate_S(self, spot, V, dW):
        t = np.linspace(0, self.T, self.s + 1)
        increments = np.sqrt(V[:,:-1]) * dW - 0.5 * V[:,:-1] * self.dt
        integral = np.cumsum(increments, axis = 1)
        S = np.zeros_like(V)
        S[:,0] = spot
        S[:,1:] = spot * np.exp(integral)

    def call_price(self, strike, S):
        last = S[:,self.s + 1]
        itm = last[last-strike>0]
        return np.sum(itm)/self.n_paths

    def price_delta(self, spot, strike, V, dW):
        def price(spot):
            S = jnp.empty(shape=(self.n_paths, self.s+1))
            S = S.at[:,0].set(spot)

            for i in range(self.s):
                inc = jnp.exp(jnp.sqrt(V[:,i]) * dW[:,i] - 0.5 * V[:,i] * self.dt)
                S = S.at[:,i+1].set(S[:,i]*inc)

            return jnp.sum(jnp.where(S[:,i+1]-strike<0,0,S[:,i+1]-strike))/self.n_paths

        return value_and_grad(price, argnums=0)(spot)


