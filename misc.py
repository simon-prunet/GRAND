import numpy as np
from scipy.signal import hilbert
import os
import sys
import matplotlib.pyplot as plt


class AntennaPlane (object):

    def __init__(self,plane_directory):

        '''
        Initializes Antenna Plane parameters (distance to injection point, antenna parameters, etc.)
        :param plane_directory: Name of directory containing all relevant files to a given plane of antennas.
        '''
        
        self.has_arrival_times = False
        self.plane_directory = plane_directory
        self.position_file = os.path.join(self.plane_directory,'antpos.dat')
        
        self.number_of_antennas = self.get_antenna_number(self.position_file)
        self.antenna_coordinates = self.get_antenna_coordinates(self.position_file)

    def compute_mean_arrival_time(self,tracefile):

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

        t1 = arrival_time[np.argmax(np.abs(hilbert(Ex)))]
        t2 = arrival_time[np.argmax(np.abs(hilbert(Ey)))]
        t3 = arrival_time[np.argmax(np.abs(hilbert(Ez)))]

        return (t1+t2+t3)/3.

    def get_arrival_times(self):
        
        '''
        Computes all mean arrival times of pulses for all antennas in the plane.
        '''

        self.arrival_times = np.zeros(self.number_of_antennas)
        for i in range(self.number_of_antennas):
            tracefile = os.path.join(self.plane_directory,'a%d.trace'%i)
            self.arrival_times[i] = self.compute_mean_arrival_time(tracefile)

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

