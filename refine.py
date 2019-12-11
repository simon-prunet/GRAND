import numpy as np
import discretize

import os, sys
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
c = 299792458. # Speed of light in m/s
const_dt = 26630.

class Simulation (object):

    def __init__(self,ZHAireS_inputfile,Xmax_file='Xmax.dat'):

        # Get injection altitude, direction of shower, direction of magnetic field from input file
        try:
            f = open(ZHAireS_inputfile,"r")
        except:
            print ("ZHAires_inputfile %s is not available.\n Aborting."%ZHAireS_inputfile)
            exit(2)

        self.ZHAireS_inputfile = ZHAireS_inputfile
        
        for line in f:
            s = line.split()
            if ('InjectionAltitude' in s):
                self.InjectionAltitude = np.float(s[1])
            if ('PrimaryZenAngle' in s):
                self.PrimaryZenAngle = np.float(s[1])
            if ('PrimaryAzimAngle' in s):
                self.PrimaryAzimAngle = np.float(s[1])
            if ('RandomSeed' in s):
                self.RandomSeed = np.float(s[1])
            if ('GeomagneticField' in s and len(s)==7):
                self.BfieldInclination = np.float(s[3])

            
        f.close()
        try:
            f = open(Xmax_file,"r")
        except:
            print ("Xmax_file %s is not available.\n Aborting."%Xmax_file)
            exit(2)
        line = f.readline()
        if ('XMAX' not in line):
            print ("Xmax_file %s is not in proper format.\n Aborting."%Xmax_file)
            exit(2)
        s = line.split("=")
        self.Xmax = np.float(s[1])
        f.close()
        

        # Compute unit vector of shower propagation. az_rad and zen_rad are seen from source perspective.
        az_rad  = np.deg2rad(180.0-self.PrimaryAzimAngle)
        zen_rad = np.deg2rad(180.0+self.PrimaryZenAngle)
        self.ShowerAxis = np.array([np.cos(az_rad)*np.sin(zen_rad),np.sin(az_rad)*np.sin(zen_rad),np.cos(zen_rad)])

        # Compute cartesian position of Xmax point
        self.Xmax_position = np.array([0.0,0.0,self.InjectionAltitude])+self.Xmax * self.ShowerAxis 

        # Compute B field unit vector
        az_bfield = 0.0 # North is magnetic in this coordinate system
        zen_bfield = np.deg2rad(90.0+self.BfieldInclination)
        self.BfieldAxis = np.array([np.cos(az_bfield)*np.sin(zen_bfield),np.sin(az_bfield)*np.sin(zen_bfield),np.cos(zen_bfield)])
        

    def in_cone(self,xyz,max_ang=3.0, from_Xmax=True, above_ground=True):

        '''
        This routine will take the shower parameters and a maximum opening angle to compute if a a given cartesian position is
        within a cone centered on the shower axis and whose origin is at the Xmax position (from_Xmax=True) or from injection position.
        This will be used for primary refinement.
        '''

        if (not from_Xmax):
            xyz_start = np.array([0.0,0.0,self.InjectionAltitude])
        else:
            xyz_start = self.Xmax_position

        dxyz = xyz - xyz_start
        if (np.linalg.norm(dxyz) < 1e-10):
            return True
        elif (np.dot(dxyz,self.ShowerAxis) > np.cos(np.deg2rad(max_ang)) * np.linalg.norm(dxyz)*np.linalg.norm(self.ShowerAxis)):
            if (above_ground):
                if (xyz[2] > 0.0): #We are above the ground
                    return True
                else:
                    return False
            else:
                return True
        else:
            return False


class Grid (Simulation):

    def __init__(self,ZHAireS_inputfile,Xmax_file='Xmax.dat',grid_width=100.0):

        super().__init__(ZHAireS_inputfile,Xmax_file)
        self.grid_width = grid_width * 1000 # in meters
        self.axis_u = np.cross(self.ShowerAxis,self.BfieldAxis)
        self.axis_u /= np.linalg.norm(self.axis_u)
        self.axis_v = np.cross(self.ShowerAxis,self.axis_u)
        self.axis_w = self.ShowerAxis

        self.x0 = self.Xmax_position - self.grid_width/2 * (self.axis_u+self.axis_v)

        self.TreeMesh = discretize.TreeMesh(h=[(self.grid_width,1)]*3,axis_u=self.axis_u,axis_v=self.axis_v,axis_w=self.axis_w,x0=self.x0)


