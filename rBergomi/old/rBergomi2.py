import numpy as np
import scipy.special
import jax.numpy as jnp
from jax import grad, value_and_grad
import jax
from itertools import product
from tqdm import tqdm_notebook
import random

class roughBergomi:
    def __init__(self, n_paths = 100000, n_steps = 52,
                    T = 1, H = 0.1, eta = 1, rho = -.3,
                    xi = 0.235**2, spot = 1.0, strike = 1.0):

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

        self.spot = spot
        self.strike = strike
        self.spots = None

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
        ZZ = self.idx_matrix()
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
                WZ[i,j] = self.rho * D * (v**(self.H+0.5)-(v-np.minimum(u,v))**(self.H+0.5))
                j += 1
            i += 1
        WW = WW + WW.T - np.diag(WW.diagonal())
        return np.block([[WW, WZ],[WZ.T, ZZ]])

    def simulate_paths(self):
        cov = self.covariance_matrix()
        mean = np.zeros(self.s*2)
        paths = np.random.multivariate_normal(mean = mean, cov = cov, size = self.n_paths)
        volterra, W = np.zeros((self.n_paths , self.s+1)), np.zeros((self.n_paths, self.s+1)) 
        volterra[:,0] = 0; volterra[:,1:] = paths[:,:self.s]
        W[:,0] = 0; W[:,1:] = paths[:,self.s:]
        dW = np.diff(W, axis = 1)
        return volterra, W, dW

    def simulate_V(self, volterra):
        t = np.linspace(0, self.T, self.s+1)
        return self.xi * np.exp(self.eta * volterra - 0.5 * self.eta**2 * t**(2 * self.H))

    def simulate_S(self, V, dW, spot=None):
        if spot is None:
            spot = self.spot
        increments = np.sqrt(V[:,:-1]) * dW - 0.5 * V[:,:-1] * self.dt
        integral = np.cumsum(increments, axis = 1)
        S = np.zeros_like(V)
        S[:,0] = spot
        if isinstance(spot, np.ndarray):
            spot = spot[:,None]
        S[:,1:] = spot * np.exp(integral)
        return S

    def call_price(self, S, strike = None):
        if strike is None:
            strike = self.strike
        last = S[:,self.s]
        itm = last[last-strike>0]
        return np.sum(itm-strike)/self.n_paths

    #def price_delta(self, V, dW, spot = None, strike = None):
    #    if spot is None:
    #        spot = self.spot
    #    if strike is None:
    #        strike = self.strike
    #    def price(spot):
    #        S = jnp.empty(shape=(self.n_paths, self.s+1))
    #        S = S.at[:,0].set(spot)
    #        for i in range(self.s):
    #            inc = jnp.exp(jnp.sqrt(V[:,i]) * dW[:,i] - 0.5 * V[:,i] * self.dt)
    #            S = S.at[:,i+1].set(S[:,i]*inc)
    #        return jnp.sum(jnp.where(S[:,i+1]-strike<0,0,S[:,i+1]-strike))/self.n_paths
    #        #return jnp.mean(jax.lax.max(S[:,i+1]-strike,0.0))
    #    return value_and_grad(price, argnums=0)(spot)

    def price_delta(self, V, dW, spot = None, strike = None):
        if spot is None:
            spot = self.spot
        if strike is None:
            strike = self.strike
        def price(spot):
            increments = jnp.sqrt(V[:,:-1]) * dW - 0.5 * V[:,:-1] * self.dt
            integral = jnp.sum(increments, axis = 1)
            S = spot * jnp.exp(integral)
            return jnp.sum(jax.lax.max(S-strike,0.0))/self.n_paths
        return value_and_grad(price, argnums=0)(spot)

    def asian(self, V, dW, spot = None, strike = None):
        if spot is None:
            spot = self.spot
        if strike is None:
            strike = self.strike
        def price(spot):
            increments = jnp.sqrt(V[:,:-1]) * dW - 0.5 * V[:,:-1] * self.dt
            integral = jnp.cumsum(increments, axis = 1)
            S = jnp.empty(shape=(self.n_paths, self.s+1))
            S = S.at[:,0].set(spot)
            S = S.at[:,1:].set(spot * jnp.exp(integral))
            return jnp.sum(jax.lax.max(jnp.mean(S, axis=1)-strike,0.0))/self.n_paths
        return value_and_grad(price, argnums=0)(spot)

    def cliquet(self, V, dW, floor=-0.1, cap=0.1, min = 0, spot = 1.0):
        if spot is None:
            spot = self.spot
        def price(min):
            increments = jnp.sqrt(V[:,:-1]) * dW - 0.5 * V[:,:-1] * self.dt
            integral = jnp.cumsum(increments, axis = 1)
            S = jnp.empty(shape=(self.n_paths, self.s+1))
            S = S.at[:,0].set(spot)
            S = S.at[:,1:].set(spot * jnp.exp(integral))
            return jnp.mean(jnp.maximum(jnp.sum(jnp.minimum(jnp.maximum((S[:,1:] / S[:,:-1]) - 1, floor), cap), axis=1), min))
        return value_and_grad(price, argnums=0)(min)

    def initialize_spots(self, vol=0.3, spot=None):
        if spot is None:
            spot = self.spot
        lognorm = np.random.lognormal
        spots = lognorm(mean=np.log(self.spot), sigma=vol, size=self.n_paths)
        self.spots = spots
        return spots

    def payoff_delta(self, V, dW, strike = None):
        if strike is None:
            strike = self.strike
        def payoff(spot, V, dW):
            increments = jnp.sqrt(V) * dW - 0.5 * V * self.dt
            integral = jnp.sum(increments)
            S = spot * jnp.exp(integral)
            return jnp.max(jnp.array([S - strike, 0.0]))
        if self.spots is None:
            spots = self.initialize_spots()
        else:
            spots = self.spots
        payoffs, deltas = [], []
        for path in tqdm_notebook(range(self.n_paths), desc='Simulating training set'):
            spot = spots[path]
            v = V[path, :-1]
            dw = dW[path,:]
            p, d = value_and_grad(payoff,argnums=0)(spot, v, dw)
            payoffs.append(p)
            deltas.append(d)
        return spots.reshape([-1,1]), np.array(payoffs).reshape([-1,1]), np.array(deltas).reshape([-1,1])

    def asian_payoff(self, V, dW, strike = None):
        if strike is None:
            strike = self.strike
        def payoff(spot, V, dW):
            increments = jnp.sqrt(V) * dW - 0.5 * V * self.dt
            integral = jnp.cumsum(increments)
            S = jnp.empty(shape=(1, self.s+1))
            S = S.at[:,0].set(spot)
            S = S.at[:,1:].set(spot * jnp.exp(integral))
            return jnp.max(jnp.array([jnp.mean(S) - strike, 0.0]))
        if self.spots is None:
            spots = self.initialize_spots()
        else:
            spots = self.spots
        payoffs, deltas = [], []
        for path in tqdm_notebook(range(self.n_paths), desc='Simulating training set'):
            spot = spots[path]
            v = V[path, :-1]
            dw = dW[path,:]
            p, d = value_and_grad(payoff,argnums=0)(spot, v, dw)
            payoffs.append(p)
            deltas.append(d)
        return spots.reshape([-1,1]), np.array(payoffs).reshape([-1,1]), np.array(deltas).reshape([-1,1])

    def cliquet_payoff(self, V, dW, floor=-0.1, cap=0.1, minmin=0, minmax=0.15, spot=1.0):
        def payoff(min, V, dW):
            increments = jnp.sqrt(V) * dW - 0.5 * V * self.dt
            integral = jnp.cumsum(increments)
            S = jnp.empty(shape=(1, self.s+1))
            S = S.at[:,0].set(spot)
            S = S.at[:,1:].set(spot * jnp.exp(integral))
            return jnp.mean(jnp.maximum(jnp.sum(jnp.minimum(jnp.maximum((S[:,1:] / S[:,:-1]) - 1, floor), cap), axis=1), min))
        mins = np.random.uniform(minmin, minmax, size=self.n_paths)
        payoffs, deltas = [], []
        for path in tqdm_notebook(range(self.n_paths), desc='Simulating training set'):
            min = mins[path]
            v = V[path, :-1]
            dw = dW[path,:]
            p, d = value_and_grad(payoff,argnums=0)(min, v, dw)
            payoffs.append(p)
            deltas.append(d)
        return mins.reshape([-1,1]), np.array(payoffs).reshape([-1,1]), np.array(deltas).reshape([-1,1])

    def cliquet_eta(self, dW, volterra, spot = 1.0, eta=1.0):
        if spot is None:
            spot = self.spot
        def price(eta):
            t = np.linspace(0, self.T, self.s+1)
            V = self.xi * jnp.exp(eta * volterra - 0.5 * eta**2 * t**(2 * self.H))
            increments = jnp.sqrt(V[:,:-1]) * dW - 0.5 * V[:,:-1] * self.dt
            integral = jnp.cumsum(increments, axis = 1)
            S = jnp.empty(shape=(self.n_paths, self.s+1))
            S = S.at[:,0].set(spot)
            S = S.at[:,1:].set(spot * jnp.exp(integral))
            return jnp.mean(jnp.sum(jnp.maximum((S[:,1:] / S[:,:-1]) - 1, 0), axis=1))
        return value_and_grad(price, argnums=0)(eta)

    def test(self, V, dW, lower=0.35, upper=1.65, n=100, n_paths = None):
        if n_paths is not None:
            n_paths_orginial = self.n_paths
            self.n_paths = n_paths
            volterra, W, dW = self.simulate_paths()
            V = self.simulate_V(volterra)
        xTest = np.linspace(lower, upper, n).reshape((-1, 1))
        yTest, dydxTest = [], []
        for spot in tqdm_notebook(xTest, desc='Simulating test set'):
            y, dydx = self.price_delta(V, dW, spot=spot)
            yTest.append(y)
            dydxTest.append(dydx)
        if n_paths is not None:
            self.n_paths = n_paths_orginial
        return xTest, np.array(yTest).reshape((-1, 1)), np.array(dydxTest).reshape((-1, 1))

    def asian_test(self, V, dW, lower=0.35, upper=1.65, n=100, n_paths = None):
        if n_paths is not None:
            n_paths_orginial = self.n_paths
            self.n_paths = n_paths
            volterra, W, dW = self.simulate_paths()
            V = self.simulate_V(volterra)
        xTest = np.linspace(lower, upper, n).reshape((-1, 1))
        yTest, dydxTest = [], []
        for spot in tqdm_notebook(xTest, desc='Simulating test set'):
            y, dydx = self.asian(V, dW, spot=spot)
            yTest.append(y)
            dydxTest.append(dydx)
        if n_paths is not None:
            self.n_paths = n_paths_orginial
        return xTest, np.array(yTest).reshape((-1, 1)), np.array(dydxTest).reshape((-1, 1))

    def cliquet_test(self, V, dW, floor = -0.15, cap = 0.15, lower=0.0, upper=0.1, n=100, n_paths = None):
        if n_paths is not None:
            n_paths_orginial = self.n_paths
            self.n_paths = n_paths
            volterra, W, dW = self.simulate_paths()
            V = self.simulate_V(volterra)
        xTest = np.linspace(lower, upper, n).reshape((-1, 1))
        yTest, dydxTest = [], []
        for min in tqdm_notebook(xTest, desc='Simulating test set'):
            y, dydx = self.cliquet(V, dW, floor=floor, cap=cap, min = min, spot=self.spot)
            yTest.append(y)
            dydxTest.append(dydx)
        if n_paths is not None:
            self.n_paths = n_paths_orginial
        return xTest, np.array(yTest).reshape((-1, 1)), np.array(dydxTest).reshape((-1, 1))

    def cliquet_eta_test(self, dW, volterra, lower=0.5, upper=2, n=100, n_paths = None):
        if n_paths is not None:
            n_paths_orginial = self.n_paths
            self.n_paths = n_paths
            volterra, W, dW = self.simulate_paths()
            V = self.simulate_V(volterra)
        xTest = np.linspace(lower, upper, n).reshape((-1, 1))
        yTest, dydxTest = [], []
        for e in tqdm_notebook(xTest, desc='Simulating test set'):
            y, dydx = self.cliquet_eta(dW, volterra, spot=self.spot, eta=e)
            yTest.append(y)
            dydxTest.append(dydx)
        if n_paths is not None:
            self.n_paths = n_paths_orginial
        return xTest, np.array(yTest).reshape((-1, 1)), np.array(dydxTest).reshape((-1, 1))