import hcone
import numpy as np


## These routines compute cartesian coordinates in (North,West,Up) frame from normalized cone coordinated X,Y,t
def cone_to_geo(X,Y,t,orig_alt=2700.,zenith_angle=120.0,azimuth=0.0,scaling=4000.,th_cone=3.0):

    '''
    Transform cone normalised coordinates to (North,West,Up) coordinates.
    X,Y are unit disc coordinates (hence the multiplication by t)
    t is a cone axis coordinate (zero at orig position), in [0,1]
    scaling is the physical length of the cone.
    '''
    disc_factor = np.tan(np.deg2rad(th_cone)/2.)
    xx=X*t*disc_factor
    yy=Y*t*disc_factor
    zz=t.copy()
    xx *= scaling; yy *= scaling; zz *= scaling
    # Compute rotation matrix. All rotations are expressed with respect to North and Up axes
    # First rotation is North(zenith_angle)
    # Second rotation is Up(azimuth)
    th = np.deg2rad(zenith_angle)
    az = np.deg2rad(azimuth)

    rot1 = np.array([[1,0,0],[0,np.cos(th),-np.sin(th)],[0,np.sin(th),np.cos(th)]])
    rot2 = np.array([[np.cos(az),-np.sin(az),0],[np.sin(az),np.cos(az),0],[0,0,1]])

    rot = np.dot(rot2,rot1)

    return np.dot(rot,np.array([xx,yy,zz])) + np.array([0,0,orig_alt])[:,None]


def geo_to_cone(xx,yy,zz,orig_alt=2700.,zenith_angle=120.0,azimuth=0.0,scaling=4000.,th_cone=3.0):
    '''
    Inverse function of the one above
    '''
    disc_factor = np.tan(np.deg2rad(th_cone)/2.)

    arr = np.array([xx,yy,zz])
    arr -= np.array([0,0,orig_alt])

    th = np.deg2rad(zenith_angle)
    az = np.deg2rad(azimuth)

    rot1 = np.array([[1,0,0],[0,np.cos(th),-np.sin(th)],[0,np.sin(th),np.cos(th)]])
    rot2 = np.array([[np.cos(az),-np.sin(az),0],[np.sin(az),np.cos(az),0],[0,0,1]])

    rot = np.dot(rot2,rot1).transpose()

    arr = np.dot(rot,arr)

    # Final scalings
    arr /= scaling
    arr /= np.array([arr[2]*disc_factor,arr[2]*disc_factor,1.])

    return arr

def refractive_index(zz,rs=325.,kr=0.0001218):
    '''
    Return air refractive index as a function of altitude
    n-1 = 1e-6*rs*exp(-kr*zz)
    '''
    return 1.0+1e-6*rs*np.exp(-kr*zz)


def example(nsamples=1000,orig_alt=2700.,scaling=4000.,zenith_angle=135.,nmax=10,mmax=4):

    from matplotlib import pyplot as plt
    from mpl_toolkits.mplot3d import Axes3D

    XX,YY,ZZ,W8,X,Y,t = hcone.cone_quad(nmax,mmax,return_XYt=True,normalized=True)
    jNlist=hcone.compute_j_N_list(nmax,mmax)

    n=refractive_index(cone_to_geo(XX,YY,ZZ,zenith_angle=zenith_angle,orig_alt=orig_alt,scaling=scaling)[2,:])
    coeffs=hcone.cone_array_noll_inverse(n,X,Y,t,jNlist,W8)
    vec=hcone.cone_array_p_noll(coeffs,XX,YY,ZZ,jNlist)

    fig=plt.figure()
    ax=fig.add_subplot(111,projection='3d')
    pl=ax.scatter(XX*ZZ,YY*ZZ,ZZ,c=n-vec)
    plt.colorbar(pl,ax=ax)
    plt.title('Residuals at quadrature points')
    plt.show()

    Xr,Yr,Zr=hcone.random_sample(nsamples)
    ref=refractive_index(cone_to_geo(Xr,Yr,Zr,zenith_angle=zenith_angle,orig_alt=orig_alt,scaling=scaling)[2,:])
    val=hcone.cone_array_p_noll(coeffs,Xr,Yr,Zr,jNlist)

    fig=plt.figure()
    ax=fig.add_subplot(111,projection='3d')
    pl=ax.scatter(Xr*Zr,Yr*Zr,Zr,c=ref-val)
    plt.colorbar(pl,ax=ax)
    plt.title('Residuals at random positions')
    plt.show()

