import aotools
import numpy as np
from solver import *

def zernike_noll(j,X,Y):

    (n,m) = aotools.zernIndex(j)
    return zernike_nm(n,m,X,Y)

def zernike_nm(n,m,X,Y):

    R = np.sqrt(X**2 + Y**2)

    if ((R>1).any()):
        print ("Coordinates should be normalized so that they pertain to the unit disc.")
        return

    phi = np.arctan2(Y,X)

    if (m==0):
        return (aotools.zernikeRadialFunc(n,0,R))
    else:

        if (m>0):
            return (np.sqrt(2.0*(n+1))*aotools.zernikeRadialFunc(n,m,R)*np.cos(m*phi))
        if (m<0):
            m=np.abs(m)
            return (np.sqrt(2.0*(n+1))*aotools.zernikeRadialFunc(n,m,R)*np.sin(m*phi))



def zernike_array_noll(coeffs,X,Y):
    '''
    Computes the direct Zernike transform \sum_j c_j Zernike_j(X,Y),
    where c_j are the coefficients 
    '''

    jmax = len(coeffs)

    res = np.zeros_like(X)
    # Noll j index starts at 1
    for j in np.arange(1,jmax+1):
        res += coeffs[j-1] * zernike_noll(j,X,Y)

    return (res)

def zernike_array_noll_transpose(vec,jmax,X,Y):
    '''
    Computes the transpose of the Zernike transform
    c_j = \sum_i vec_i Zernike_j(X_i,Y_i) = <vec_i,Zernike_j(X,Y)>
    '''
    res = np.zeros(jmax)

    for j in np.arange(1,jmax+1):
        res[j-1] = np.dot(vec,zernike_noll(j,X,Y))
    return (res)

def gram_matrix(coeffs,X,Y):

    '''
    Computes A^T.A.c where A is the direct Zernike transform
    '''

    jmax = len(coeffs)
    vec = zernike_array_noll(coeffs,X,Y)
    res = zernike_array_noll_transpose(vec,jmax,X,Y)

    return (res)

def compute_zernike_coeffs(vec,jmax,X,Y):

    '''
    Compute Zernike coefficients from vector using normal equations
    '''
    rhs = zernike_array_noll_transpose(vec,jmax,X,Y)
    res = cg_solve(rhs,gram_matrix,np.zeros(jmax),args=(X,Y))

    return (res)


