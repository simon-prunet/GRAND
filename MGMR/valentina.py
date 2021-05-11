import numpy as np
import matplotlib.pyplot as plt

Ep = 1e17 # eV
Ne = 6*(Ep/1e10)
C = -np.log(0.63)/4000.
X0 = 36.7 # g/cm2
Xmx = 840. + 70.*np.log10(Ep/1e20) # g/cm2
J=-0.04 * 1.44e-9 * Ne

def ctret(ct,d):
    return (-d**2/(2.*ct))

def XX(z):
    return 1000.*np.exp(-C*z)

def ss(X):
    return 3.*(X/X0)/(X/X0 + 2*Xmx/X0)

def dssdXX(X):
    return 6.*(Xmx/X0**2) / (X/X0 + 2*Xmx/X0)**2

def ft(ct):
    z = -ct
    X = XX(z)
    s = ss(X) 
    res = np.exp( (X-Xmx-1.5*X*np.log(s))/X0 )
    return res

def dftdct(ct):
    z = -ct
    X = XX(z)
    s = ss(X)
    dsdX = dssdXX(X)
    res = ft(ct)/X0 * (1 - 1.5*np.log(s) - 1.5*X/s*dsdX) * C * X
    return (res)

def E(ct,d):
    ctr = ctret(ct,d)
    fac = J*4.*ctr**2/d**4
    res = fac * (ctr * dftdct(ctr) + ft(ctr))
    return(res)

d=300.0
ct = np.linspace(1,50,1001)
# d=700.
# ct = np.linspace(1,150,1001)
ctr = ctret(ct,d)
plt.plot(ct,E(ct,d)*1e6)
plt.xlabel(r'$ct(m)$')
plt.ylabel(r'$E_x(\mu V)$')
plt.show()