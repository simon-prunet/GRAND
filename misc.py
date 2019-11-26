import numpy as np
from scipy.signal import hilbert
import os
import sys
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
c = 299792458. # Speed of light in m/s
const_dt = 26630.

class AntennaPlane (object):

    def __init__(self,plane_directory,position_file='antpos.dat',Xmax_file='Xmax.dat',RASPASS_outfile='RASPASSdata.out'):

        '''
        Initializes Antenna Plane parameters (distance to injection point, antenna parameters, etc.)
        :param plane_directory: Name of directory containing all relevant files to a given plane of antennas.
        '''
        
        self.has_arrival_times = False
        self.plane_directory = plane_directory
        self.position_file = os.path.join(self.plane_directory,position_file)
        if (not os.path.exists(self.position_file)):
            print("Cannot read antenna position %s\n Aborting."%self.position_file)
        self.Xmax_file = os.path.join(self.plane_directory,Xmax_file)
        if (not os.path.exists(self.Xmax_file)):
            print("Cannot read Xmax file %s\n Aborting."%self.Xmax_file)
        self.RASPASS_outfile = os.path.join(self.plane_directory,RASPASS_outfile)
        if (not os.path.exists(self.RASPASS_outfile)):
            print("Cannot read RASPASS out file %s\n Aborting."%self.RASPASS_outfile)

        self.dist_to_Xmax = self.get_dist_to_xmax(self.Xmax_file)
        
        self.decay, self.axis = self.get_decay_axis(self.RASPASS_outfile) # Coordinates of injection point, and unit shower vector
        self.Xmax_pos = self.decay + self.dist_to_Xmax*self.axis # Coordinates of Xmax point

        # Exponential index law parameters
        self.Rs = 325.
        self.Kr = 0.1218 * 1e-3 # m^-1


        self.number_of_antennas = self.get_antenna_number(self.position_file)
        self.antenna_coordinates = self.get_antenna_coordinates(self.position_file)

    def compute_ref_arrival_time(self,antenna_pos):
        '''
        Given the antenna coordinates, and coordinates of Xmax position, computes reference 
        arrival time at antenna of a fiducial ray emitted from Xmax position.
        '''
        # We start by computing the sine of the zenithal angle of the fiducial light ray emanating
        # from the Xmax position and reaching later tha antenna
        ray_vec = (self.Xmax_pos - antenna_pos)
        n_ray_vec = np.linalg.norm(ray_vec)
        if (n_ray_vec < 1e-15):
            return 0.
        ray_vec /= n_ray_vec
        sinaz = np.dot(ray_vec,np.array([0.,0.,1.]))
        hX = self.Xmax_pos[2]
        hA = antenna_pos[2]
        # Note that delta_t_XA should be always positive, as hX-hA and sinaz have the same sign
        delta_t_XA = ( (hX-hA) - self.Rs/self.Kr*(np.exp(-self.Kr*hX)-np.exp(-self.Kr*hA))*1e-6 )/(c*sinaz)
        return delta_t_XA*1e9 #in ns

    def compute_mean_arrival_time(self,tracefile,use_hilbert=False):

        '''
        Computes mean arrival time over the three E field components from antenna trace file.
        File format is t(ns) | Ex (uV) | Ey (uV) | Ez (uV)
        Uses max of Hilbert enveloppe for each components
        :param tracefile: Input ZHAireS antenna trace file
        returns arrival time (ns)
        '''

        try:
            arr=np.loadtxt(tracefile)
        except:
            print("Problem loading file %s\n Aborting."%tracefile)
            sys.exit(2)

        arrival_time = arr[:,0]
        Ex = arr[:,1]
        Ey = arr[:,2]
        Ez = arr[:,3]

        power = np.sqrt(Ex**2+Ey**2+Ez**2)
        if (power.sum() < 1):
            # Not enough signal to safely determine arrival time
            return np.nan
        if (use_hilbert):
            t1 = arrival_time[np.argmax(np.abs(hilbert(Ex)))]
            t2 = arrival_time[np.argmax(np.abs(hilbert(Ey)))]
            t3 = arrival_time[np.argmax(np.abs(hilbert(Ez)))]
            return (t1+t2+t3)/3.

        else: # Use power max
            t = arrival_time[np.argmax(power)]
            return t

    def get_arrival_times(self,relative=True,use_hilbert=False):
        
        '''
        Computes all mean arrival times of pulses for all antennas in the plane.
        '''

        self.arrival_times = np.zeros(self.number_of_antennas)
        for i in range(self.number_of_antennas):
            tracefile = os.path.join(self.plane_directory,'a%d.trace'%i)
            self.arrival_times[i] = self.compute_mean_arrival_time(tracefile,use_hilbert=use_hilbert)
            if (relative):
                # Subtract arrival time of fiducial ray emitted from Xmax position to antenna, taking into account ref. index
                self.arrival_times[i] -= self.compute_ref_arrival_time(self.antenna_coordinates[i,:]) + const_dt

        self.has_arrival_times=True


    def get_antenna_number(self,position_file):

        '''
        Returns the number of antennas from a file giving the antennas coordinates (antpos.dat)
        :param position_file: input file giving the antennas coordinates in m.
        '''

        try:
            f=open(position_file,'r')
        except:
            print('Cannot read from file %s\n Aborting.'%position_file)
            sys.exit(2)

        na = len(f.readlines())
        return na

    def get_antenna_coordinates(self,position_file):

        '''
        Returns a numpy array of shape (#antennas,3) containing the antenna coordinates in m, where the origin
        is defined in x,y relative to the injection (decay) position, while the z origin is at sea level.
        This is the reference frame for all Radio Morphing calculations.
        :param position_file: input file giving the antennas coordinates in m.
        '''
        na = self.get_antenna_number(position_file)
        coordinates = np.loadtxt(position_file)
        if (coordinates.shape[0] != na):
            print("Something went wrong when reading antenna coordinates from file %s\n Aborting."%position_file)
            sys.exit(2)

        return coordinates

    def get_dist_to_xmax(self,Xmax_file):
        '''
        Returns distance from Xmax point in meters
        '''
        f = open(Xmax_file,'r')
        line = f.readline()
        f.close
        xmax = np.float(line.split('=')[-1].split('\n')[0])
        return xmax

    def get_decay_axis(self,RASPASS_outfile):
        '''
        Return coordinates of injection point (decay), and coordinates of (unit) shower axis vector
        '''
        f = open(RASPASS_outfile,'r')
        decay = np.zeros(3)
        axis = np.zeros(3)
        for line in f:
            spl = line.split('=')
            key=spl[0]
            val=np.float(spl[-1].split('\n')[0])
            if (key=='INJECTIONX'):
                decay[0]=val
            if (key=='INJECTIONY'):
                decay[1]=val
            if (key=='INJECTIONZ'):
                decay[2]=val
            if (key=='INJECTIONAXISX'):
                axis[0]=val
            if (key=='INJECTIONAXISY'):
                axis[1]=val
            if (key=='INJECTIONAXISZ'):
                axis[2]=val

        return decay, axis


    def visualize_arrival_times(self):
        
        '''
        Returns a 3D scatter plot of the antenna positions, with color as a function of arrival times
        '''
        if (not self.has_arrival_times):
            self.get_arrival_times()
        fig = plt.figure()
        ax = fig.add_subplot(111,projection='3d')
        sc = ax.scatter(self.antenna_coordinates[:,0],self.antenna_coordinates[:,1],self.antenna_coordinates[:,2],c=self.arrival_times,cmap=plt.cm.get_cmap('RdYlBu'))
        plt.colorbar(sc)
        plt.show()


class AllAntennas (object):

    def __init__(self,simu_dir,index_file='MasterIndex',position_file='antpos.dat',Xmax_file='Xmax.dat',RASPASS_outfile='RASPASSdata.out'):

        self.simu_dir = simu_dir
        self.index_file = os.path.join(self.simu_dir,index_file)
        # First, get number of antenna planes and name of plane directories
        if (not os.path.exists(self.index_file)):
            print ('Index file %s does not exist \n Aborting.'%self.index_file)
            sys.exit(2)
        f = open(self.index_file,'r')
        self.nplanes = len(f.readlines()) - 1 # First line is header
        f.close()
        self.plane_dirs = np.empty((self.nplanes,),dtype=object)
        self.planes = np.empty((self.nplanes,),dtype=AntennaPlane)
        f = open(self.index_file,'r')
        f.readline() # Get rid of header
        for i in range(self.nplanes):
            plane_dir = f.readline().split()[0]
            self.plane_dirs[i] = os.path.join(self.simu_dir,plane_dir)
            self.planes[i] = AntennaPlane(self.plane_dirs[i])
        #
        self.has_arrival_times = False

    def visualize_arrival_times(self,use_hilbert=False,zoom=4.,elevation=6.,azimuth=-134.,x_fig_inch=18.3,y_fig_inch=11.8):
        
        '''
        Returns a 3D scatter plot of the antenna positions, with color as a function of arrival times
        '''
        if (not self.has_arrival_times):
            for i in range(self.nplanes):
                self.planes[i].get_arrival_times(use_hilbert=use_hilbert)
            self.has_arrival_times = True
        fig = plt.figure(figsize=(x_fig_inch,y_fig_inch))
        ax = fig.add_subplot(111,projection='3d')
        if (azimuth is not None and elevation is not None):
            ax.view_init(elev=elevation,azim=azimuth)

        for i in range(self.nplanes):
            pl = self.planes[i]
            sc = ax.scatter(pl.antenna_coordinates[:,0],pl.antenna_coordinates[:,1],pl.antenna_coordinates[:,2],c=pl.arrival_times,cmap=plt.cm.get_cmap('RdYlBu'))
        xlim = ax.get_xlim3d(); ylim=ax.get_ylim3d(); zlim=ax.get_zlim3d()
        xc = (xlim[0]+xlim[1])/2.; yc = (ylim[0]-ylim[1])/2.; zc = (zlim[0]+zlim[1])/2.
        dx = (xlim[1]-xlim[0])/2.; dy = (ylim[1]-ylim[0])/2.; dz = (zlim[1]-zlim[0])/2.
        dxyz = np.max((dx,dy,dz))/zoom
        ax.auto_scale_xyz([xc-dxyz,xc+dxyz],[yc-dxyz,yc+dxyz],[zc-dxyz,zc+dxyz])
        plt.colorbar(sc)
        plt.show()
        return fig



