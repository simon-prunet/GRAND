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
                self.BfieldInclination = s[3]
            
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

        # Compute B field unit vector
        az_bfield = 0.0 # North is magnetic in this coordinate system
        zen_bfield = np.deg2rad(90.0+self.Bfield_Inclination)
        self.BfieldAxis = np.array([np.cos(az_bfield)*np.sin(zen_bfield),np.sin(az_bfield)*np.sin(zen_bfield),np.cos(zen_bfield)])
        
