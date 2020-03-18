import numpy as np
from alegendre import *
from scipy.optimize import brentq
from math import ceil


class basis:

    def __init__(self,kmax=20,mmax=4,theta0=3.0):
        
        # These values will be used to determine sampling, and max values 
        self.kmax=kmax
        self.mmax=mmax
        # Cone aperture
        self.theta0 = theta0
        self.orders = np.ndarray((kmax+1,mmax+1))
        self.orders.fill(np.nan)
        # Initialize computations
        self.dsize = np.zeros(1)
        alegendreeval.alegendre_eval_init(dsize)
        # Compute n_k's
        for m in range(mmax+1):
            self.orders[:,m] = self.compute_orders(np.deg2rad(self.theta0),m,self.kmax)
            
    def phase_function(self,dnu,dmu,theta):

        mydnu = np.array(dnu)
        mydmu = np.array(dmu)
        mytheta = np.array(theta)
        alpha = np.zeros(1)
        alphader = np.zeros(1)
        vallogp = np.zeros(1)
        vallogq = np.zeros(1)
        valp = np.zeros(1)
        valq = np.zeros(1)
        dsize = np.zeros(1)

        alegendreeval.alegendre_eval(mydnu,mydmu,mytheta,alpha,alphader,vallogp,vallogq,valp,valq)
        return (alpha)

    def phase_counter(self,dnu,dmu,j,theta):

        phase = self.phase_function(dnu,dmu,theta)
        return (phase-2.5*np.pi - np.pi*(j-1))

    def compute_orders(self,theta,dmu,kmax=20):

        orders = np.zeros(kmax+1)
        numin = np.max((10.0,dmu))
        for j in range(kmax+1):
            numax=2*numin
            while(self.phase_counter(numax,dmu,j+1,theta)<0):
                numax += numin
            nuk = brentq(self.phase_counter,numin,numax,args=(dmu,j+1,theta))
            orders[j] = nuk
            print ('nuk = %f'%nuk)

        return (orders)

    def compute_roots(self,dnu,dmu,theta):
        
        nproots = np.zeros(1,dtype=int)
        nqroots = np.zeros(1,dtype=int)
        alegendreeval.alegendre_nroots(dnu,dmu,nproots,nqroots)
        roots=np.zeros(nproots)
        Cone.alegendre_array_proots(dnu,dmu,roots)
        ou = np.where(roots<theta)
        return roots[ou]

    def compute_pnumu(self,dnu,dmu,thetas):

        res = np.zeros_like(thetas)
        Cone.alegendre_array_p_eval(dnu,dmu,thetas,res)
        return res

    def compute_jacobi_quad(self,dmu,n):

        nn=int(ceil(n/2.))
        t_array = np.zeros(nn)
        w_array = np.zeros(nn)
        Cone.alegendre_array_jacobi(dmu,n,t_array,w_array)

        return t_array,w_array

    def compute_scaled_jacobi_quad(self,n,theta):

        # Scales the roots and weights for integration with t in [0,theta]
        t,w = self.compute_jacobi_quad(0.0,n)
        if (n%2): 
            #odd. pi/2 is part of the roots, appears only once
            t = np.hstack((t,np.pi-t[:-1][::-1]))
            w = np.hstack((w,w[:-1][::-1]))
        else:
            t = np.hstack((t,np.pi-t[::-1]))
            w = np.hstack((w,w[::-1]))

        ctheta = np.cos(theta)
        st = np.arccos((1-ctheta)/2. * np.cos(t) + (1+ctheta)/2.)
        sw = (1-ctheta)/2. * w
        return st,sw

    def grid_coordinates(self,phis,thetas):

        # Construct a 2D set of coordinates in phi, theta, tensor product of the inputs
        # Fast varying coordinate is phi
        ph,th = np.meshgrid(phis,thetas)
        return ph, th

    def generate_map(self,k,m_in,amplitude):
        # Simple single coefficient map
        if (m_in<=0):
            cos=True
            m=-m_in
        else:
            cos=False
            m=m_in
        if (m>self.mmax):
            print ('Required value of m=%d is larger than mmax=%d'%(m,self.mmax))
            return
        if (k>self.kmax):
            print ('Required value of k=%d is larger than kmax=%d'%(k,self.kmax))
            return

        thetas, weights = self.compute_scaled_jacobi_quad(self.kmax,self.theta0)
        phis = np.arange(2*self.mmax)*np.pi/float(self.mmax)
        if (cos):
            azfunc = np.cos(m*phis)
        else:
            azfunc = np.sin(m*phis)
        thfunc = self.compute_pnumu(self.orders[k,m],m,thetas)
        return np.outer(thfunc,azfunc)

    def show_map(self,map,theta=3.0,up=50):

        dims = map.shape
        # Upscale by a factor 10
        import cv2
        from matplotlib import pyplot as plt
        fmap = np.fft.fft(map,axis=-1)
        tmp = np.zeros((dims[0],up*dims[1]),dtype=np.complex)
        tmp[:,:dims[1]/2+1] = fmap[:,:dims[1]/2+1]
        tmp[:,-dims[1]/2+1:] = fmap[:,-dims[1]/2+1:]
        im2 = np.fft.ifft(tmp)*up
        im3 = cv2.resize(im2.real,(up*dims[1],up*dims[0]))
        #return(im3)

        plt.subplot(projection='polar')
        th = np.linspace(0,theta,up*dims[0])
        ph = np.linspace(0,2*np.pi,up*dims[1])
        T,P = np.meshgrid(th,ph)
        plt.pcolormesh(P,T,im3.T)
        plt.colorbar()
        plt.show()
        return(im3)
#def analyse(map):




