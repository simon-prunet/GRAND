################################
#### by A. Zilles
################################
#!/usr/bin/env python

import sys

from sys import argv
import os
import glob

import time

import numpy as np
from numpy import *

#### Astropy logger
#from astropy import log
#log.setLevel("INFO")
#log.info("...")

import logging
logging.basicConfig(filename="example_usingclass.log", level=logging.INFO)
logger = logging.getLogger('Main')

import tqdm

import matplotlib.pyplot as plt
from matplotlib.pyplot import cm 
#from mpl_toolkits.mplot3d import Axes3D

#from astropy import units as u
#muV_m = u.u * u.V / u.m

from astropy.table import Table
from astropy.table import hstack
import h5py

from os.path import split, join, realpath
root_dir = realpath(join(split(__file__)[0], "..")) # = $PROJECT
sys.path.append(join(root_dir, "lib", "python"))
#import radio_simus 
from radio_simus.signal_treatment import p2p
from radio_simus.shower import *
from radio_simus.detector import detector, create_from_file, get_array, get_slopes, find_antennaposition, find_antennaslope
from radio_simus.__init__ import arrayfile 



###############################
#### WORK in PROGRESS: define new functions and test

##def find_antennaposition(positions, antIDs, ID):
    ##index = np.where(antIDs == int(ID))[0][0]
    ###print(index)
    ##return positions[index]

#def find_antennaID(positions, antIDs, pos):
    #pos=np.around(pos)
    #print(type(positions), type(pos))
    #print(positions[0], pos)
    #index = np.where((positions[:]== pos).all())
    #print('indexx', index, positions[index,0])
    #return antIDs[index]
### not yet working
##pos_ID = find_antennaID(positions, antIDs, np.asarray(f.meta["position"]))
##print("----", f.meta["position"], pos_ID, "original", ID)
#################################

if __name__ == '__main__':
    
    if ( len(sys.argv)<2 ):
        print("""
        Example on how to use classes
        -- Analysis trigger for events, create a list of events (class objects) and  trigger 1/0 to class attributes
        -- create a png with statistic for triggering
        
        Use: python3 example_usingclass.py <folder event set>
        Example: python3 example_usingclass.py ../../CoREAS/GP300_centered/
        
        NOTE: still ongoing work, already usable
        
        ATTENTION:
        -- adapt the paths given in the config-file so that eg the antenna reponse can be read-in correctly
        """)
        sys.exit(0)
    
    
    

    # path to folder containing the inp-file, trace files and antpos.dat 
    eventfolder = sys.argv[1]

    #create  "empty detector"
    det = detector()
    #create detector=antenna array from file defined in config file
    create_from_file(det, arrayfile)

    ## get all antennas positions
    #array = get_array(det)
    ## get all slopes
    #slopes = get_slopes(det)

    ## antenna IDs
    #antIDs = array.T[0]
    ## get only position info
    #positions = array[:,1:4] 

    event = [] # python list


    #with logger.log_to_file(eventfolder +'/ana_trigger.log', filter_level='INFO'):

    print("\nScan of events ...")
    # loop over all folder
    for path in tqdm.tqdm(glob.glob(eventfolder+"/*/")):
        if os.path.isdir(path): # only pick event folders
            logger.debug("... Reading Event from:"+ path)

            
            # loop over all antenna positions in event
            i=0
            trigger_any=[]
            trigger_xy=[]
                

            for file in glob.glob(path+"/*.hdf5"):
                ## find antenna position and its slope per ID - works
                ID = int(file.split('/')[-1].split('.hdf5')[0].split('table_')[-1])
                #pos_ant = find_antennaposition(det, ID)
                #pos_slope = find_antennaslope(det, ID)
                
                if i==0: # just get the first antenna to readin meta info            
                    #### EFIELD just as example how to read in
                    f = Table.read(file, path='efield') 
                    #print("\n simulated position ", f.meta["position"])
                    
                    # create shower object and set attributes
                    testshower = sim_shower()
                    loadInfo_toShower(testshower, f.meta)
                    param = testshower.get_all() # get all parameters, all call them separately
                    logger.info("   SUMMARY EVENT: ShowerID = "+  str(param[0])+ " primary = "+ str(param[1])+ " energy/eV = "+ str(param[2]) + " zenith/deg = "+ str(param[3])+ " azimuth/deg = "+ str(param[4])+ " injectionheight/m = "+ str(param[5]) )
                        
                    event.append(testshower)
                i+=1
                    
                    
                #### VOLTAGES - add trigger info to shower object
                try:
                    g = Table.read(file, path='voltages') 
                    
                    # info: trigger = [any_aggr, xz_aggr, thr_aggr, any_cons, xy_cons, thr_cons]
                    # Here: only ask for aggressive value for triggering
                    if g.meta["trigger"][1] ==1:
                        trigger_any.append(ID)
                    if g.meta["trigger"][2] ==1:
                        trigger_xy.append(ID)
                except IOError:
                    logger.error("Voltages not computed for antenna: "+ str(ID) +" in "+path)
                    
                # check whether voltages exists, if not compute voltage
                
                # check whether antenna ID position and slope already exits, otherwise load to detector
                
                
                
            ## EXAMPLE: Trigger Analysis
            if len(trigger_any)>5 or len(trigger_xy)>5:
                logger.info("   => shower triggers (aggr): any =" + str(len(trigger_any)) + " xy = " + str(len(trigger_xy)))
                # add trigger info to class
                event[-1].add_trigger(1)
            else:
                testshower.add_trigger(0)
                
                    
            # plot full array (gray), simulated position with voltages and mark triggers (red) for each event, save as png
            ## could be similar to example_plot_2D
            
        else: 
            continue
    
    
    ###### START ANALYSIS ################### 
    print("\nStart an analysis ...")    
    #print(event[0].showerID, event[0].trigger)
    Event_ID=map(lambda i: i.showerID, event)

    ### Calculate Ratio of detected events
    trigger=list(map(lambda i: i.get_trigger(), event))
    print("\n"+str(sum(trigger))+" out of "+str(len(trigger))+" events detected --> "+str(100.* sum(trigger)/len(trigger))+"% detection rate"+"\n")    

    ### find triggered events
    trigger=np.asarray(trigger)
    index = np.where(trigger==1)[0]

    # parameters
    energy=np.asarray(list(map(lambda i: i.get_energy(), event)))
    zenith=np.asarray(list(map(lambda i: i.get_zenith(), event)))
    azimuth=np.asarray(list(map(lambda i: i.get_azimuth(), event)))
    primary=np.asarray(list(map(lambda i: i.get_primary(), event)))

    # Plot
    plt.rcParams.update({'figure.figsize':(12,5)})

    fig, (ax1, ax2, ax3, ax4) = plt.subplots(1, 4)
    ax1.hist(energy, color="orange", bins=50, alpha=0.5, label='all = '+str(len(trigger)))
    ax1.hist(energy[index], bins=50, label='trigger (aggr.) = '+str(sum(trigger)))
    ax1.legend()
    ax1.set_title('energy in eV')
    ax2.hist(zenith, color="orange", bins=50, alpha=0.5)
    ax2.hist(zenith[index], bins=50)
    ax2.set_title('zenith in deg')
    ax3.hist(azimuth, color="orange", bins=50, alpha=0.5)
    ax3.hist(azimuth[index], bins=50)
    ax3.set_title('azimuth in deg')
    ax4.hist(primary, color="orange", bins=50, alpha=0.5)
    ax4.hist(primary[index], bins=50)
    ax4.set_title('primary')

    fig.tight_layout()
    plt.savefig(eventfolder+"/trigger_stats.png")
    print("PNG saved:" + eventfolder+"/trigger_stats.png")
    logger.info("PNG saved:" + eventfolder+"/trigger_stats.png")
    plt.show()    



    #====== end of run =======      
    logger.info("Done within "+str(time.clock()) +"s")    
    
    
    
