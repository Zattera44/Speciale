import numpy as np
from scipy.stats import norm
from scipy.optimize import brentq

### Most of the code is from ###
# https://github.com/differential-machine-learning/notebooks/blob/master/DifferentialMLTF2.ipynb
# Adjustments have been made in order to augment the training set with sigma

def bsPrice(spot, strike, vol, T):
    d1 = (np.log(spot/strike) + 0.5 * vol * vol * T) / vol / np.sqrt(T)
    d2 = d1 - vol * np.sqrt(T)
    return spot * norm.cdf(d1) - strike * norm.cdf(d2)

def bsImpVol(P, spot, strike, T):
    P = np.maximum(P, np.maximum(spot - strike, 0))
    def error(s):
        return bsPrice(spot, strike, s, T) - P
    return brentq(error, 1e-9, 1e+9)
    
def bsDelta(spot, strike, vol, T):
    d1 = (np.log(spot/strike) + 0.5 * vol * vol * T) / vol / np.sqrt(T)
    return norm.cdf(d1)

def bsVega(spot, strike, vol, T):
    d1 = (np.log(spot/strike) + 0.5 * vol * vol * T) / vol / np.sqrt(T)
    return spot * np.sqrt(T) * norm.pdf(d1)
    
class BlackScholes:
    
    def __init__(self, 
                 vol=0.2,
                 T1=1, 
                 T2=2, 
                 K=1,
                 volMult=1.5):
        
        self.spot = 1
        self.vol = vol
        self.T1 = T1
        self.T2 = T2
        self.K = K
        self.volMult = volMult
                        
    def trainingSet(self, m):
            
        returns = np.random.normal(size=[m, 2])

        vol0 = self.vol * self.volMult
        R1 = np.exp(-0.5*vol0*vol0*self.T1 + vol0*np.sqrt(self.T1)*returns[:,0])
        sigma = abs(np.random.normal(scale=0.15,size=m)) #np.random.uniform(0.05, 0.5, m)
        R2 = np.exp(-0.5*sigma**2*(self.T2-self.T1) \
                    + sigma*np.sqrt(self.T2-self.T1)*returns[:,1])
        S1 = self.spot * R1
        S2 = S1 * R2 

        pay = np.maximum(0, S2 - self.K)
        R2a = np.exp(-0.5*sigma**2*(self.T2-self.T1) \
                    - sigma*np.sqrt(self.T2-self.T1)*returns[:,1])
        S2a = S1 * R2a             
        paya = np.maximum(0, S2a - self.K)
        
        X = S1
        Y = 0.5 * (pay + paya)

        Z1 =  np.where(S2 > self.K, R2, 0.0).reshape((-1,1)) 
        Z2 =  np.where(S2a > self.K, R2a, 0.0).reshape((-1,1)) 
        Z = 0.5 * (Z1 + Z2)
        V1 = np.where(S2 > self.K, S2*(-sigma*(self.T2-self.T1)+np.sqrt(self.T2-self.T1)*returns[:,1]), 0.0).reshape((-1,1))
        V2 = np.where(S2a > self.K, S2a*(-sigma*(self.T2-self.T1)+np.sqrt(self.T2-self.T1)*(-returns[:,1])), 0.0).reshape((-1,1))
        V = 0.5 * (V1 + V2)

        return X.reshape([-1,1]), sigma.reshape([-1,1]), Y.reshape([-1,1]), Z.reshape([-1,1]), V.reshape([-1,1])
    
    def testSet(self, lower=0.35, upper=1.65, num=100):
        
        spots = np.linspace(lower, upper, num).reshape((-1, 1))
        prices = bsPrice(spots, self.K, self.vol, self.T2 - self.T1).reshape((-1, 1))
        deltas = bsDelta(spots, self.K, self.vol, self.T2 - self.T1).reshape((-1, 1))
        return spots, prices, deltas   