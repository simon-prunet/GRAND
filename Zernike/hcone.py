import numpy as np
import zernike as zk
import scipy.special as ss
import aotools

def jacobi_quad(nmax):

    # Here 2 = d-1 = 3-1 in 3D
    xx,ww = ss.roots_jacobi(nmax+1,0,2)
    t = (xx+1.)/2.

    return (t, ww)


def cone_jacobi(N,n,t):
    '''
    returns the orthogonal polynomials on the vertical direction 
    for the cone interior. Here t in [0,1]
    '''

    jac = ss.jacobi(N-n,0,2+2*n)
    res = jac(2*t-1) * t**n / np.sqrt(8./(2*N+3.))

    return res



def cone_polynomial(N,n,m,X,Y,t):

    '''
    Cone polynomial, orthogonal on cone interior
    Based on Zernike polynomials on the disc and 
    Jacobi polynomials for the cone axis
    Here X,Y are assumed to be in the unit disc, and t in [0,1],
    i.e. these are normalized coordinates
    '''

    cj = cone_jacobi(N,n,t)
    zz = zk.zernike_nm(n,m,X,Y)

    return np.outer(cj,zz).flatten()

def cone_polynomial_noll(N,j,X,Y,t):
    '''
    same as before, but Noll index instead of (n,m)
    '''
    n,m = aotools.zernIndex(j)
    cj = cone_jacobi(N,n,t)
    zz = zk.zernike_noll(j,X,Y)

    return np.outer(cj,zz).flatten()


def cone_p_noll(N,j,XX,YY,ZZ):
    '''
    Same as before, but input coordinates are not factorized
    '''
    n,m = aotools.zernIndex(j)
    cj = cone_jacobi(N,n,ZZ)
    zz = zk.zernike_noll(j,XX,YY)

    return cj*zz

def compute_j_N_list(nmax,mmax):

    '''
    Computes the list of quantum numbers N,n,m such that
    abs(m)<=min(n,mmax)
    0<=n<=N<=nmax
    '''
    ll=[]
    for n in range(nmax+1):

        jlistn = zk.compute_j_list_n(n,mmax)
        # For each element of jlistn corresponding to a different m, we have 
        # nmax+1-n cone jacobi terms
        for j in jlistn:
            n,m = aotools.zernIndex(j)
            for N in range(n,nmax+1):
                ll.append((j,N))

    return ll


def cone_quad(nmax,mmax,normalized=False,return_XYt=True):

    ''' 
    Computes cone samples and quadrature weights for orthonormal polynomials
    on the cone interior
    Contiguous samples are within each disc
    If normalized is set to True, each disc is normalized to the unit disc, these are the
    coordinates that are relevant for Zernike polynomials. Only use unnormalized coordinates
    for display purposes.
    If return_XYt is set to true, return independently X,Y coordinates of unit disc and t coordinate of cone axis.
    These are the coordinates used to compute the cone polynomials.
    '''
    X,Y,W = zk.zernike_quad(nmax,mmax)
    t,wt = jacobi_quad(nmax)
    W8 = np.outer(wt,W).flatten()
    if (normalized):
        XX = np.outer(np.ones_like(t),X).flatten()
        YY = np.outer(np.ones_like(t),Y).flatten()
    else:
        XX = np.outer(t,X).flatten()
        YY = np.outer(t,Y).flatten()
    ZZ = np.outer(t,np.ones_like(X)).flatten()

    if (return_XYt):
        return (XX,YY,ZZ,W8,X,Y,t)
    else:
        return (XX,YY,ZZ,W8)

def cone_array_noll_inverse(vec,X,Y,t,jNlist,W8):

    '''
    compute inverse cone polynomial transform using orthonormal polynomials on quadrature set
    '''
    res = np.zeros(len(jNlist))
    for i,(j,N) in enumerate(jNlist):
        res[i] = np.dot(vec,cone_polynomial_noll(N,j,X,Y,t)*W8)

    return (res)

def cone_array_noll(coeffs,X,Y,t,jNlist):
    '''
    Computes the direct cone polynomial transform \sum_Nj c_Nj Cone_Nj(X,Y,t),
    where c_Nj are the coefficients and N,j \in jNlist.
    '''

    # Make sure coeffs and js are of same length
    if (len(jNlist)!=len(coeffs)):
        print ("list of js and coeffs should be the same size")
        return

    res = np.zeros(X.size*t.size)
    for i,(j,N) in enumerate(jNlist):
        res += coeffs[i] * cone_polynomial_noll(N,j,X,Y,t)

    return (res)

def cone_array_p_noll(coeffs,XX,YY,ZZ,jNlist):
    '''
    Same as before, but coordinates are not factorized
    '''
    # Make sure coeffs and js are of same length
    if (len(jNlist)!=len(coeffs)):
        print ("list of js and coeffs should be the same size")
        return

    res = np.zeros(XX.size)
    for i,(j,N) in enumerate(jNlist):
        res += coeffs[i] * cone_p_noll(N,j,XX,YY,ZZ)

    return (res)


def random_sample(n):
    '''
    Samples uniformly in a cone, returning samples in normalized coordinates
    '''
    x = np.random.uniform(size=n)
    r = np.sqrt(x)
    phi = np.random.uniform(high=2.0*np.pi,size=n)
    z = np.random.uniform(size=n)

    X = r*np.cos(phi)
    Y = r*np.sin(phi)
    Z = np.power(z,1./3.)

    return X,Y,Z

