import numpy as np
import scipy.special
from scipy.stats import norm
from scipy.optimize import brentq
import jax.numpy as jnp
from jax import value_and_grad
from itertools import product
from tqdm import tqdm_notebook

class roughBergomi:
    def __init__(self
                ,n_paths = 100000
                ,n_steps = 52
                ,T = 1):

        self.n_paths = n_paths
        self.n_steps = n_steps
        self.T = T

        self.dt = 1.0/n_steps
        self.s = int(n_steps * T)
        self.t = np.linspace(T/self.s,T,self.s)
   
    def G(self
            ,x
            ,H): 

        F = scipy.special.hyp2f1(1.0,0.5-H,1.5+H,1.0/x)
        return (2*self.H*(x**(H-0.5))*scipy.special.gamma(0.5+H)*F)/scipy.special.gamma(H+1.5)

    def idx_matrix(self):

        v, u = np.zeros(pow(self.s, 2)), np.zeros(pow(self.s, 2))

        k = 0
        for i in product(self.t, self.t):
            v[k], u[k] = i
            k = k + 1

        return np.reshape(np.minimum(v, u), (self.s, self.s))

    def covariance_matrix(self
            ,H
            ,rho):

        ZZ = self.idx_matrix()
        WW = np.zeros((self.s, self.s))
        WZ = np.zeros((self.s, self.s))
        D = np.sqrt(2.0*H)/(H+0.5)

        i = 0
        for v in self.t:

            j = 0
            for u in self.t:

                if i == j:
                    WW[i,j] = u**(2.0*H)

                if i > j:
                    x = v/u
                    WW[i,j] = u**(2.0*H)*self.G(x, H)

                if i < j:
                    WW[i,j] = 0

                WZ[i,j] = rho * D * (v**(H+0.5)-(v-np.minimum(u,v))**(H+0.5))

                j += 1
            i += 1

        WW = WW + WW.T - np.diag(WW.diagonal())
        return np.block([[WW, WZ],[WZ.T, ZZ]])

    def simPaths(self, H = 0.1, rho = -0.3):

        self.H = H
        self.rho = rho

        cov = self.covariance_matrix(H, rho)
        mean = np.zeros(self.s*2)
        paths = np.random.multivariate_normal(mean = mean, cov = cov, size = self.n_paths)

        volterra, W = np.zeros((self.n_paths , self.s+1)), np.zeros((self.n_paths, self.s+1)) 

        volterra[:,0] = 0
        volterra[:,1:] = paths[:,:self.s]

        W[:,0] = 0
        W[:,1:] = paths[:,self.s:]

        dW = np.diff(W, axis = 1)

        return volterra, W, dW

    def simV(self
            ,volterra
            ,xi = 0.235**2
            ,eta = 1):

        t = np.linspace(0, self.T, self.s+1)
        return xi * np.exp(eta * volterra - 0.5 * eta**2 * t**(2 * self.H))

    def simS(self
            ,V
            ,dW
            ,spot=1.0):

        increments = np.sqrt(V[:,:-1]) * dW - 0.5 * V[:,:-1] * self.dt
        integral = np.cumsum(increments, axis = 1)

        S = np.zeros_like(V)
        S[:,0] = spot

        if isinstance(spot, np.ndarray):
            spot = spot[:,None]

        S[:,1:] = spot * np.exp(integral)

        return S

    def callPrice(self, S, strike = 1.0):
        last = S[:,self.s]
        itm = last[last-strike>0]
        return np.sum(itm-strike)/self.n_paths

    def priceDeriv(self
            ,dW
            ,volterra
            ,spot = 1.0
            ,strike = 1.0
            ,xi = 0.235**2
            ,eta = 1
            ,payoff = 'European'
            ,deriv = 0):

        def price(spot, xi):

            t = np.linspace(0, self.T, self.s+1)
            V = xi * jnp.exp(eta * volterra - 0.5 * eta**2 * t**(2 * self.H))

            if payoff != 'European':

                increments = jnp.sqrt(V[:,:-1]) * dW - 0.5 * V[:,:-1] * self.dt
                integral = jnp.cumsum(increments, axis = 1)

                S = jnp.empty(shape=(self.n_paths, self.s+1))
                S = S.at[:,0].set(spot)
                S = S.at[:,1:].set(spot * jnp.exp(integral))

                if payoff == 'Asian':
                    return jnp.mean(jnp.maximum(jnp.mean(S, axis=1) - strike, 0.0))

                if payoff == 'Cliquet':
                    return jnp.mean(jnp.sum(jnp.maximum((S[:,1:] / S[:,:-1]) - 1, 0), axis=1))

            else:

                increments = jnp.sqrt(V[:,:-1]) * dW - 0.5 * V[:,:-1] * self.dt
                integral = jnp.sum(increments, axis = 1)
                S = spot * jnp.exp(integral)

                return jnp.mean(jnp.maximum(S - strike, 0.0))

        return value_and_grad(price, argnums=deriv)(spot, xi)

    def initialize_spots(self
            ,vol=0.3
            ,spot=1.0):

        lognorm = np.random.lognormal
        return lognorm(mean=np.log(spot), sigma=vol, size=self.n_paths)
        
    def trainingSet(self
            ,dW
            ,volterra
            ,spot = 1.0
            ,strike = 1.0
            ,xi = 0.235**2
            ,eta = 1
            ,payoff = 'European'
            ,deriv = 0):

        def pathwise_payoff(spot, xi, volterra, dW, payoff):

            t = np.linspace(0, self.T, self.s+1)
            V = xi * jnp.exp(eta * volterra - 0.5 * eta**2 * t**(2 * self.H))

            if payoff != 'European':

                increments = jnp.sqrt(V[:-1]) * dW - 0.5 * V[:-1] * self.dt
                integral = jnp.cumsum(increments)

                S = jnp.empty(shape=(1, self.s+1))
                S = S.at[:,0].set(spot)
                S = S.at[:,1:].set(spot * jnp.exp(integral))

                if payoff == 'Asian':
                    return jnp.mean(jnp.maximum(jnp.mean(S, axis=1) - strike, 0.0))

                if payoff == 'Cliquet':
                    return jnp.mean(jnp.sum(jnp.maximum((S[:,1:] / S[:,:-1]) - 1, 0)))

            else:

                increments = jnp.sqrt(V[:-1]) * dW - 0.5 * V[:-1] * self.dt
                integral = jnp.sum(increments)
                S = spot * jnp.exp(integral)

                return jnp.mean(jnp.maximum(S - strike, 0.0))
        
        if deriv == 0:
            spots = self.initialize_spots(spot=spot)
            xis = np.repeat(xi, self.n_paths)
            x = spots
        if deriv == 1: 
            spots = np.repeat(spot, self.n_paths)
            xis = np.random.uniform(low=0.01, high=0.15, size=self.n_paths)
            x = xis
       
        payoffs, derivs = [], []
        for path in tqdm_notebook(range(self.n_paths), desc='Simulating training set'):
            p, d = value_and_grad(pathwise_payoff,argnums=deriv)(spots[path], xis[path], volterra[path, :], dW[path,:], payoff)
            payoffs.append(p)
            derivs.append(d)

        return x.reshape([-1,1]), np.array(payoffs).reshape([-1,1]), np.array(derivs).reshape([-1,1])

    def testSet(self
            ,dW
            ,volterra
            ,lower = 0.35
            ,upper = 1.65
            ,n = 100
            ,spot = 1.0
            ,strike = 1.0
            ,xi = 0.235**2
            ,eta = 1
            ,payoff = 'European'
            ,deriv = 0):

        if payoff != 'Cliquet':
            spots = np.linspace(lower, upper, n).reshape((-1, 1))
            xTest = spots
            xis = np.repeat(xi, n)
        else:
            xis = np.linspace(lower, upper, n).reshape((-1, 1))
            xTest = xis
            spots = np.repeat(spot, n)

        yTest, dydxTest = [], []
        for i in tqdm_notebook(range(n), desc='Simulating test set'):
            y, dydx = self.priceDeriv(dW, volterra
                                    ,spot = spots[i], strike = strike
                                    ,xi = xis[i], eta = eta
                                    ,payoff = payoff, deriv = deriv)
            yTest.append(y)
            dydxTest.append(dydx)
        return xTest, np.array(yTest).reshape((-1, 1)), np.array(dydxTest).reshape((-1, 1))    

def bsPrice(spot, strike, vol, T):
    d1 = (np.log(spot/strike) + 0.5 * vol * vol * T) / vol / np.sqrt(T)
    d2 = d1 - vol * np.sqrt(T)
    return spot * norm.cdf(d1) - strike * norm.cdf(d2)

def bsImpVol(P, spot, strike, T):
    P = np.maximum(P, np.maximum(spot - strike, 0))
    def error(s):
        return bsPrice(spot, strike, s, T) - P
    return brentq(error, 1e-9, 1e+9)