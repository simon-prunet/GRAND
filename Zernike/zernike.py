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

def gram_matrix(coeffs,X,Y,js):

    '''
    Computes A^T.A.c where A is the direct Zernike transform
    '''

    # Make sure coeffs and js are of same length
    if (len(js)!=len(coeffs)):
        print ("list of js and coeffs should be the same size")
        return
    vec = zernike_array_noll(coeffs,X,Y,js)
    res = zernike_array_noll_transpose(vec,X,Y,js)

    return (res)

def compute_zernike_coeffs(vec,X,Y,js):

    '''
    Compute Zernike coefficients from vector using normal equations
    '''
    rhs = zernike_array_noll_transpose(vec,X,Y,js)
    res = cg_solve(rhs,gram_matrix,np.zeros(len(js)),args=(X,Y,js))

    return (res)

def compute_j_list(nmax,mmax):
    '''
    Construct a list of Noll indices with n<=nmax, abs(m)<=mmax
    '''
    l=[]
    (n,m)=(0,0)
    j=1
    while (n<=nmax):
        (n,m) = aotools.zernIndex(j)
        if (abs(m)<=mmax and n<=nmax):
            l.append(j)
        j+=1

    return (l)



