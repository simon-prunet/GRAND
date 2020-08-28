import aotools
import numpy as np
from solver import *

def zernike_noll(j,X,Y):

    (n,m) = aotools.zernIndex(j)
    return zernike_nm(n,m,X,Y)

def zernike_radial_n(n,m,r):
    return (np.sqrt(2.0*(n+1)) * aotools.zernikeRadialFunc(n,m,r))

def zernike_nm(n,m,X,Y):

    R = np.sqrt(X**2 + Y**2)

    if ((R>1).any()):
        print ("Coordinates should be normalized so that they pertain to the unit disc.")
        return

    phi = np.arctan2(Y,X)

    if (m==0):
        return (np.sqrt(n+1.)*aotools.zernikeRadialFunc(n,0,R))
    else:

        if (m>0):
            return (np.sqrt(2.0*(n+1))*aotools.zernikeRadialFunc(n,m,R)*np.cos(m*phi))
        if (m<0):
            m=np.abs(m)
            return (np.sqrt(2.0*(n+1))*aotools.zernikeRadialFunc(n,m,R)*np.sin(m*phi))



def zernike_array_noll(coeffs,X,Y,js):
    '''
    Computes the direct Zernike transform \sum_j c_j Zernike_j(X,Y),
    where c_j are the coefficients and j \in js.
    '''

    # Make sure coeffs and js are of same length
    if (len(js)!=len(coeffs)):
        print ("list of js and coeffs should be the same size")
        return

    res = np.zeros_like(X)
    i=0
    for j in js:
        res += coeffs[i] * zernike_noll(j,X,Y)
        i+=1

    return (res)

def zernike_array_noll_transpose(vec,X,Y,js):
    '''
    Computes the transpose of the Zernike transform
    c_j = \sum_i vec_i Zernike_j(X_i,Y_i) = <vec_i,Zernike_j(X,Y)>
    '''

    res = np.zeros(len(js))
    i=0
    for j in js:
        res[i] = np.dot(vec,zernike_noll(j,X,Y))
        i+=1
    return (res)

def gram_matrix(coeffs,X,Y,js,weights=1.0,regul=0):

    '''
    Computes A^T.A.c where A is the direct Zernike transform
    '''

    # Make sure coeffs and js are of same length
    if (len(js)!=len(coeffs)):
        print ("list of js and coeffs should be the same size")
        return
    vec = zernike_array_noll(coeffs,X,Y,js)
    res = zernike_array_noll_transpose(weights*vec,X,Y,js)
    res += regul*coeffs
    return (res)

def compute_zernike_coeffs(vec,X,Y,js,weights=None,regul=0):

    '''
    Compute Zernike coefficients from vector using normal equations
    Possibility to give input weights per antenna location.
    '''
    if (weights is None):
        weights = np.ones_like(vec)
    rhs = zernike_array_noll_transpose(weights*vec,X,Y,js)
    res = cg_solve(rhs,gram_matrix,np.zeros(len(js)),args=(X,Y,js,weights,regul))

    return (res)

def compute_j_list(nmax,mmax_in):
    '''
    Construct a list of Noll indices with n<=nmax, abs(m)<=mmax
    '''
    l=[]
    mmax = np.min((mmax_in,nmax))
    (n,m)=(0,0)
    j=1
    while (n<=nmax):
        (n,m) = aotools.zernIndex(j)
        if (abs(m)<=mmax and n<=nmax):
            l.append(j)
        j+=1

    return (l)

def compute_j_list_n(n,mmax_in):
    '''
    Construct a list of Noll indices of a given n, for all valid (abs(m)<=n)
    '''
    mmax = np.min((mmax_in,n))
    l=[]
    if ((n-mmax)%2):
        m=-mmax+1
    else:
        m=-mmax

    while (m<=mmax):
        j=(n*(n+1))/2+abs(m)
        tt = n%4 < 2
        if (m>=0 and not tt):
            j += 1
        if (m<=0 and tt):
            j += 1
        l.append(j)
        m += 2

    return l


### The two following routines define respectively a quadrature set (points and weights) for the Zernike polynomials on the unit
### disc, as well as a Zernike transform that relies on the Zernike orthogonality on said quadrature set.

def zernike_quad(nmax,mmax):
    '''
    Construct a list of quadrature points and associated weights. These must be adapted to compute
    scalar products involved in the analysis (Zernike transform), hence must have >= 2mmax+1 radial branches,
    and >= nmax+1 radial points (which allows integration of polynomials up to degree 2nmax+1). Note that
    the radial measure is proportional to radius for Zernike polynomials orthogonality, hence the shifted Jacobi polynomial
    P_{nmax+1}^{0,1}(2r-1) is chosen for the radial quadrature.
    '''

    # First define radial points and weights
    from scipy.special import roots_jacobi
    
    xx,w8 = roots_jacobi(nmax+1,0,1) # xx in [-1,1]
    r = (xx+1.)/2.
    theta = np.arange(2*mmax+1)/(2.*mmax+1.)*2.*np.pi
    ctheta = np.cos(theta)
    stheta = np.sin(theta)

    X = np.outer(r,ctheta).flatten()
    Y = np.outer(r,stheta).flatten()
    # Normalize weights: factor of 4 for shifted polynomial variable change, factor of (2*mmax+1)/2. for angular integration
    w8 /= 2.*(2*mmax+1.)
    W = np.outer(w8,np.ones_like(ctheta)).flatten()

    return X,Y,W

def zernike_array_noll_inverse(vec,X,Y,js,W):
    ''' 
    Computes the Zernike coefficients from sampled 2D map using quadrature points and weights
    '''
    res = np.zeros(len(js))
    
    for i,j in enumerate(js):
        res[i] = np.dot(vec,zernike_noll(j,X,Y)*W)

    return (res)





