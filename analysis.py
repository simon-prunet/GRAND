import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from Simulations import hdf5fileinout as hdf5io
from Zernike import zernike
jmax=16 #(nmax=4)

################################################################################
# Test of Matias functions for simulation and hdf5 files handling (Valentin)
################################################################################

#Shower path
simu_path = '/Users/prunet/Documents/Code/GRAND/Simulations/Stshp_XmaxLibrary_0.03981_48.19_90_Gamma_24.hdf5'
InputFilename = simu_path

EventNumber = 0 # Only one simulation here
#Run info
RunInfo = hdf5io.GetRunInfo(InputFilename)
print("#RunInfo --> ", RunInfo)
EventName = hdf5io.GetEventName(RunInfo,EventNumber)
print("#EventName --> ", EventName)


#Antannas info
AntennaInfo = hdf5io.GetAntennaInfo(InputFilename,EventName)
print("#AntennaInfo --> ", AntennaInfo)
IDs = AntennaInfo['ID'].data
print("#IDs --> ", IDs)

#Get filter on IDs name, to get index of antennas in the star shape array
ou=[id.startswith('A') for id in IDs]

#Get X,Y coordinates of all antennas
X = hdf5io.GetXFromAntennaInfo(AntennaInfo)
Y = hdf5io.GetYFromAntennaInfo(AntennaInfo)

#Get peak Hilbert time and amplitude
peaktime, peakamplitude = hdf5io.get_peak_time_hilbert_hdf5(InputFilename, antennamax="All",antennamin=0, usetrace="efield", DISPLAY=False)

# Select data from star shape
Xref = X[ou]; Xref=Xref.data
Yref = Y[ou]; Yref=Yref.data

# Normalize coordinates
Xscale = 1.1 * (Xref.max()-Xref.min())
Xrefn = Xref / Xscale
Yscale = 1.1 * (Yref.max()-Yref.min())
Yrefn = Yref / Yscale

peaktimeref = peaktime[ou]
peakamplituderef = peakamplitude[ou]

# Compute Zernike coefficients
peaktime_coeffs = zernike.compute_zernike_coeffs(peaktimeref,jmax,Xrefn,Yrefn)
peakamplitude_coeffs = zernike.compute_zernike_coeffs(peakamplituderef,jmax,Xrefn,Yrefn)

fitted_peaktimeref = zernike.zernike_array_noll(peaktime_coeffs,Xrefn,Yrefn)
fitted_peakamplituderef = zernike.zernike_array_noll(peakamplitude_coeffs,Xrefn,Yrefn)


fig=plt.figure()
ax=fig.add_subplot(121,projection='3d')
ax.scatter3D(Xref,Yref,peaktimeref)
ax2=fig.add_subplot(122,projection='3d')
ax2.scatter3D(Xref,Yref,peaktimeref-fitted_peaktimeref)

plt.show()

fig2=plt.figure()
axx=fig2.add_subplot(121,projection='3d')
axx.scatter3D(Xref,Yref,peakamplituderef)
axx2=fig2.add_subplot(122,projection='3d')
axx2.scatter3D(Xref,Yref,fitted_peakamplituderef)

plt.show()
