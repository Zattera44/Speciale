#%%

import math
import pandas as pd
import numpy as np
from scipy import stats
from scipy.optimize import fsolve, minimize
from scipy import integrate
import matplotlib.pyplot as plt
import seaborn as sns

#%%

def BSM(S, K, r, sigma, tau):
    d1 = (np.log(S/K) + (r + np.power(sigma, 2) /2) * tau)/(sigma * np.sqrt(tau))
    d2 = d1 - sigma * np.sqrt(tau)
    price = stats.norm.cdf(d1) * S - stats.norm.cdf(d2) * K * np.exp(-r * tau)
    return price

#%%

def ImpliedVol(S, K, r, price, tau):
    f = lambda sigma: price - BSM(S, K, r, sigma, tau)
    return fsolve(f, [0.1])[0]

#%%

def delta(S, r, data):
    d1 = (np.log(S/data['Strike']) + (r + np.power(data['Implied Vol'], 2)/2)*data['Expiry'] / 365)/(data['Implied Vol']*np.sqrt(data['Expiry'] / 365))
    delta = stats.norm.cdf(d1)
    return delta

#%%

def vega(S, r, data):
    d1 = (np.log(S/data['Strike']) + (r + np.power(data['Implied Vol'], 2)/2)*data['Expiry'] / 365)/(data['Implied Vol']*np.sqrt(data['Expiry'] / 365))
    vega = S * np.sqrt(data['Expiry'] / 365) * stats.norm.pdf(d1) / 100
    return vega

#%%

data = pd.read_csv('GOOGLData.csv', sep=';')

#%%

data.head(10)

#%%

#Defining parameters from the assignment
S = 2786
r = 0

#%%

length = data.shape[0]
Est_Vol = np.zeros((length,1))
for i in range(length):
    Est_Vol[i] = ImpliedVol(S, data.loc[i][1], r, data.loc[i][0], data.loc[i][2] / 365)

#%%

Est_Vol = pd.DataFrame(Est_Vol, columns = ['Implied Vol'])
Est_Vol.head(10)

#%%

data = pd.concat([data,Est_Vol], axis = 1)

#%%

data.head(9)

#%%

data['X_Plot'] = np.log(data['Strike']/S)

#%%

hues = sns.color_palette('mako', 5)

#%%

sns.set_style("darkgrid", {"axes.facecolor": ".9"})

sns.lineplot(x='X_Plot', y='Implied Vol', data=data, hue='Expiry', palette = hues)

plt.show()

#%%

#Calculation of Delta and Vega
data['Delta'] = delta(S, r, data)
data['Vega'] = vega(S, r, data)

#%%

sns.lineplot(x='X_Plot', y='Delta', data=data, hue='Expiry', palette = hues)

#%%

sns.lineplot(x='X_Plot', y='Vega', data=data, hue='Expiry', palette = hues)

#%% md


## Opgave 2.1

#%%

def charFunction(u, S, tau, v, sigma, kappa, theta, rho, r = 0, div = 0):
    d = np.sqrt((rho * sigma * u * 1j - kappa) ** 2 + sigma ** 2 * (u*1j + u ** 2))
    d = -d
    g = (kappa - rho * sigma * u * 1j + d) / (kappa - rho * sigma * u * 1j - d)
    M = (r - div) * u*1j*tau + (kappa * theta)/(sigma **2) * ((kappa - rho*sigma*u*1j + d)*tau - 2 * np.log((1-g*np.exp(d*tau))/(1-g)))
    N = (kappa - rho*sigma*u*1j + d)/(sigma ** 2) *((1-np.exp(d*tau))/(1-g*np.exp(d*tau)))
    result = np.exp(M + N*v + 1j * u * np.log(S))
    return result

#%%

def HestonCall(S, K, tau, v, sigma, kappa, theta, rho, r=0, div=0):
    def integral1(u):
        char1 = charFunction(u - 1j, S, tau, v, sigma, kappa, theta, rho)
        char2 = charFunction(-1j, S, tau, v, sigma, kappa, theta, rho)
        z = (np.exp(-1j * u * np.log(K)) * char1 / char2) / (1j * u)
        return z.real
    def integral2(u):
        char = charFunction(u, S, tau, v, sigma, kappa, theta, rho)
        z2 = (np.exp(-1j * u * np.log(K)) * char) / (1j * u)
        return z2.real
    Pi1 = 0.5 + 1/math.pi * integrate.quad(integral1, 0, np.inf, limit = 100)[0]
    Pi2 = 0.5 + 1/math.pi * integrate.quad(integral2, 0, np.inf, limit = 100)[0]
    result = S * np.exp(-div*tau)*Pi1 - K * np.exp(-r * tau)*Pi2
    #print(integrate.quad(integral2, 0, np.inf, limit=50, epsabs=0.01)[1])
    #print(tau, v, sigma, kappa, theta, rho)
    return result

#%%



def loss(x0):
    err = 0
    for i in range(len(data)):
        err = err + ((data.loc[i][0] - HestonCall(S, data.loc[i][1], data.loc[i][2] /365, x0[0], x0[1], x0[2], x0[3], x0[4])) / (data.loc[i][3] * 100) ) ** 2
    return err

#%%



x0 = np.array([0.2 , 1., 2., 0.2 ** 2, -0.25])
bnds = np.array([[0.01, 0.99], [0.01, 10], [0.01, 10], [0.01, 0.99], [-1, 1]])


#%%



x0 = np.array([0.2, 1., 2., 0.2, -0.2])
res = minimize(loss, x0, bounds=bnds, options={'maxiter':200, 'disp':True}) #tol=1e-2
print(res)

print(HestonCall(S, 2500, 1, 0.02, 10.0, 0.01, 0.01, -1.0))


