import numpy as np
import mpmath as mp
import scipy.integrate as si
from scipy.special import erfc

R_earth = 6371.0 # km
Rs=315.0
k=0.1218 #km-1

def norm(x):
    '''
    Compute Euclidian norm of coordinate vector
    '''
    return np.sqrt(x[0]**2+x[1]**2+x[2]**2)

def compute_cos_zenith(xa,xb):
    '''
    Compute cosine of zenith angle at point A
    xa and xb are the cartesian coordinates of points A and B
    Only assumes origin of coordinates is at earth center, so that 
    xa is radial vector going upwards at point A.
    Light goes from B to A here.
    '''
    xab = xb-xa
    nxab = norm(xab)
    nxa  = norm(xa)
    coszen = (xa[0]*xa[0]+xab[1]*xa[1]+xab[2]*xa[2])/(nxab*nxa)

    return coszen

def s(h,R,coszen):
    '''
    Compute distance to point A along AB segment, as a fonction of local altitude
    '''
    return (-R*coszen + np.sqrt(h**2 + 2.0*h*R + R**2*coszen**2))

def dsdh(h,R,coszen):
    '''
    Return derivative of distance to point A along AB segment, w.r.t. local altitude
    '''
    return((R+h)/np.sqrt(h**2 + 2.0*h*R + R**2*coszen**2))

def approx_dsdh(h,R,coszen):
    '''
    Return second order expansion in h/R of derivative of distance to point A along AB segment, 
    w.r.t. local altitude
    '''
    prf = 1.0/np.abs(coszen)
    hoR = h/R
    cos2zen = coszen**2
    tan2zen = (1.0-cos2zen)/cos2zen
    return (prf*( 1.0 - tan2zen*hoR +3.0/2.0*tan2zen/cos2zen*hoR**2))


def compute_geometry(xa,xb):
    '''
    Origin of coordinates at earth center
    z axis goes up at point A
    '''

    coszen = compute_cos_zenith(xa,xb)
    R = norm(xa) # Effective radius (earth plus altitude of point A from ground)
    hB = norm(xb) - R # Altitude of point B w.r.t. point A

    return (coszen, R, hB)

def refractivity(h):

    ''' 
    Return refractivity as a function of height
    '''
    return (Rs * np.exp(-k*h))

def exact_refractivity_integral(xa,xb):
    '''
    Will compute the exact integral of the refractivity on the AB segment.
    The integral will be computed numerically, using the exact formula for ds/dh
    '''

    # First, compute zenith angle cosine at point A, R=||OA||, and 
    coszen, R, hB = compute_geometry(xa,xb)
    dh = R - R_earth # Altitude of point A
    
    # The following will multiply the refractivity integral, as the refractivity height dependence is computed from ground level,
    # while the height in the integral is w.r.t. point A 
    prefactor = np.exp(-k*dh) 

    def integrand(h):
        return (prefactor * refractivity(h) * dsdh(h,R,coszen))

    return (si.quad(integrand,0,hB)[0])

def approximate_refractivity_integral(xa,xb):
    '''
    Same as before, but dsdh exact formula is Taylor expanded to second order in h/R
    where R is the effective radius (||OA||) and h is the altitude w.r.t. point A
    '''
    coszen, R, hB = compute_geometry(xa,xb)
    dh = R - R_earth
    prefactor = np.exp(-k*dh)

    
    prefactor *= Rs/(k*np.fabs(coszen))
    cos2zen = coszen**2
    tan2zen = (1.-cos2zen)/cos2zen
    ee =  np.exp(-k*hB)
    return (prefactor*( 1 - ee -tan2zen*(1-(1+k*hB)*ee)/(k*R) +tan2zen/cos2zen*(6-(6+6*k*hB+3*k**2*hB**2)*ee)/(k*R)**2 )) 


###############################################

def s(h,zen):
    s = -R_earth*np.cos(zen) + np.sqrt(h**2 + 2*h*R_earth + np.cos(zen)**2 * R_earth**2)
    return s

def exact_integral(zendeg,hx):

    zen = np.deg2rad(zendeg)
    sx = s(hx,zen)
    def integrand(s):
        h = -R_earth + np.sqrt(R_earth**2+s**2+2*R_earth*s*np.cos(zen))
        return np.exp(-k*h)

    res = Rs * si.quad(integrand,0,sx)[0]
    return(res)

def inclined_approx(zendeg,hx):

    ''' 
    For use with 5<zendeg<89.99 degrees
    '''
    zen = np.deg2rad(zendeg)
    sx = s(hx,zen)

    def beta0(w):
        res = np.cosh(w)/np.sinh(w)-0.5/np.sinh(w/2.)
        return res

    def beta1(w):
        res = 1./np.sinh(w)-np.cosh(w)**2/np.sinh(w)**3 +1./(8.*np.sinh(w/2.)**3) -3./(16.*np.sinh(w/2.))
        return res

    alpha0 = 1.0
    alpha1 = 3./8.

    def IK1(u,z,w):
        
        smallw = np.sqrt(np.pi/(2.*z))*(alpha0 + alpha1/z)*erfc(np.sqrt(2.*z)*np.sinh(w/2.))*np.exp(-z+u)
        largew = (beta0(w)/z + beta1(w)/z**2)*np.exp(-z*np.cosh(w)+u)
        return smallw+largew
    
    u = k*R_earth
    z = u*np.sin(zen)
    w1 = np.arcsinh(np.cos(zen)/np.sin(zen))
    w2 = np.arcsinh( (sx + R_earth*np.cos(zen))/(R_earth*np.sin(zen)) )
    return Rs *(R_earth*np.sin(zen)) *( IK1(u,z,w1) - IK1(u,z,w2))





