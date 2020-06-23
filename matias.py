import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from Simulations import hdf5fileinout as hdf5io
from Starshape.geom import LinePlaneCollision, Simplerot
from Zernike import zernike
import numpy as np

azim=-66
elev=10
xfigsize=16
yfigsize=10
ref_sample = 368
################################################################################
# Test of Matias functions for simulation and hdf5 files handling (Valentin)
################################################################################
def TranslateTrace(trace,deltat,step=0.5):

    '''
    Use FFT to translate a trace by deltat seconds. Sampling rate is given by step in ns.
    '''

    nsamp = trace.size
    freq = np.fft.fftfreq(nsamp,d=step)
    ftrace = np.fft.fft(trace)
    phase = np.exp(-1j * 2.0*np.pi * freq * deltat)
    return np.fft.ifft(ftrace*phase)

def ComputeDisplacement(trace_ref, trace, step=0.5):
    '''
    Uses phase correlation method to compute relative displacement. Result is given in ns.
    '''
    nsamp = trace_ref.size
    ftrace_ref = np.fft.fft(trace_ref)
    nftrace_ref = np.abs(ftrace_ref)
    ftrace = np.fft.fft(trace)
    nftrace = np.abs(ftrace)
    cr = ftrace.conj()*ftrace_ref / (nftrace_ref * nftrace)
    # Remove low SNR frequencies
    mask = np.where(nftrace*nftrace_ref < 0.1*np.max(nftrace*nftrace_ref))
    cr[mask] = 0.
    icr = (np.fft.ifft(cr)).real
    # Locate peak of phase correlation
    imax = np.argmax(icr)
    if (imax > nsamp/2):
        imax -= nsamp
    # Locate max of neighboring values
    # imaxp1 = (imax+1)%nsamp # Periodic boundaries 
    # imaxm1 = (imax-1)%nsamp # Periodic boundaries
    # choice = np.argmax(np.array(icr[imaxm1],icr[imaxp1]))
    # if (choice == 0): # Secondary peak is on the left
    #     sign = -1
    #     isec = imaxm1
    # else:
    #     sign = 1
    #     isec = imaxp1
    # subdt = sign * icr[isec]/(icr[imax] + icr[isec])
    # return (imax+subdt) * step
    return imax * step

#Shower path
simu_path = '/Users/prunet/Documents/Code/GRAND/Starshape/Stshp_Proton_3.98_77.4_0.0_1.hdf5'
InputFilename = simu_path

EventNumber = 0 # Only one simulation here
#Run info
RunInfo = hdf5io.GetRunInfo(InputFilename)
print("#RunInfo --> ", RunInfo)
EventName = hdf5io.GetEventName(RunInfo,EventNumber)
print("#EventName --> ", EventName)

#Shower parameters
Zenith = hdf5io.GetEventZenith(RunInfo,EventNumber)
print("#Zenith --> ", Zenith)
#Azimuth = hdf5io.GetEventAzimuth(RunInfo,EventNumber)
#print("#Azimuth --> ", Azimuth)
XmaxDistance = hdf5io.GetEventXmaxDistance(RunInfo,EventNumber)
print("#XmaxDistance --> ", XmaxDistance)

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

# Crosscheck antennas
nou=[not k for k in ou]
Xcr = X[nou]; Xcr=Xcr.data
Ycr = Y[nou]; Ycr=Ycr.data

# Compute predictions from Zernike models
peaktimecross = peaktime[nou]
peakamplitudecross = peakamplitude[nou]

# Now compute coordinates of antennas reprojected into plane perpendicular to shower axis and passing by origin
el_angle = 180.-Zenith
mat = Simplerot(el_angle) # Uses elevation angle
na_tot = X.size
Xproj = np.zeros(na_tot)
Yproj = np.zeros(na_tot)
origin = np.zeros(3)
planeNormal = np.array([np.sin(np.deg2rad(el_angle)),0.0,np.cos(np.deg2rad(el_angle))])
xmax_pos = XmaxDistance * planeNormal # 3D position of Xmax
for i in range(na_tot):
    # 3D position of antenna
    antenna_pos = np.array([X[i],Y[i],0.0])
    # Unit vector along light ray from antenna to Xmax
    rayDirection = xmax_pos - antenna_pos; rayDirection /= np.linalg.norm(rayDirection)
    # Compute intersection of light ray and plane perpendicular to shower axis
    intersection_pos = LinePlaneCollision(planeNormal,origin,rayDirection,antenna_pos)
    # Now rotate it so that X,Y coordinates lie in the sample plane perpendicular to shower axis
    rotated_pos = np.dot(mat,intersection_pos)
    Xproj[i] = rotated_pos[0]
    Yproj[i] = rotated_pos[1]
    # print i,rotated_pos[2]

# Show scatter plot of reprojected antenna coordinates in plane perpendicular to shower axis
fig = plt.figure()
ax = fig.gca()
ax.set_title('Reprojected antenna positions')
ax.scatter(Xproj,Yproj)
plt.show()

# Extract reference and cross antennas, and normalize coordinates to unit disc
Xprojref = Xproj[ou]
Yprojref = Yproj[ou]
Xprojcross = Xproj[nou]
Yprojcross = Yproj[nou]

Xscale = 1.01*np.max((Xprojref.max(),-Xprojref.min()))
Xprojrefn = Xprojref / Xscale
Xprojcrossn = Xprojcross / Xscale
Yscale = 1.01*(Yprojref.max()-Yprojref.min())/2.
Yprojrefn = Yprojref / Yscale
Yprojcrossn = Yprojcross / Yscale

fig = plt.figure()
ax = fig.gca()
ax.set_title('Reprojected antenna positions (normalized)')
ax.scatter(Xprojrefn,Yprojrefn)
plt.show()

# Now extract the Efield traces, and store them per component as an ntime,nantennas vector
# First get traces for one antenna to determine the number of time bins
tmp = hdf5io.GetAntennaEfield(InputFilename,EventName,IDs[0])
ntime = tmp.shape[0]
tmp = hdf5io.GetAntennaFilteredVoltage(InputFilename,EventName,IDs[0])
ntimeV = tmp.shape[0]

# Efield and time arrays
Efieldx = np.zeros((ntime,na_tot))
Efieldy = np.zeros((ntime,na_tot))
Efieldz = np.zeros((ntime,na_tot))
time    = np.zeros((ntime,na_tot))

# Voltage and time arrays
Voltagex = np.zeros((ntimeV,na_tot))
Voltagey = np.zeros((ntimeV,na_tot))
Voltagez = np.zeros((ntimeV,na_tot))
timeV    = np.zeros((ntimeV,na_tot))

# Displacement vectors
Efieldx_dt = np.zeros(na_tot)
Efieldy_dt = np.zeros(na_tot)
Efieldz_dt = np.zeros(na_tot)

for i in range(na_tot):

    # E field
    trace = hdf5io.GetAntennaEfield(InputFilename,EventName,IDs[i])
    time[:,i] = trace[:,0]
    Efieldx[:,i] = trace[:,1]
    Efieldy[:,i] = trace[:,2]
    Efieldz[:,i] = trace[:,3]
    # Antenna voltages
    traceV = hdf5io.GetAntennaFilteredVoltage(InputFilename,EventName,IDs[i])
    timeV[:,i]    = traceV[:,0]
    Voltagex[:,i] = traceV[:,1]
    Voltagey[:,i] = traceV[:,2]
    Voltagez[:,i] = traceV[:,3]

    # Use Hilbert envelop peak time and FFT to translate all traces such that their peaks are aligned
    # Efieldx_dt[i] = ComputeDisplacement(Efieldx[:,0],Efieldx[:,i])
    # Efieldx[:,i] = TranslateTrace(Efieldx[:,i],Efieldx_dt[i])
    # Efieldy_dt[i] = ComputeDisplacement(Efieldy[:,0],Efieldy[:,i])
    # Efieldy[:,i] = TranslateTrace(Efieldy[:,i],Efieldy_dt[i])
    # Efieldz_dt[i] = ComputeDisplacement(Efieldz[:,0],Efieldz[:,i])
    # Efieldz[:,i] = TranslateTrace(Efieldz[:,i],Efieldz_dt[i])
    Efieldx[:,i] = TranslateTrace(Efieldx[:,i],-peaktime[i]+time[ref_sample,i])
    Efieldy[:,i] = TranslateTrace(Efieldy[:,i],-peaktime[i]+time[ref_sample,i])
    Efieldz[:,i] = TranslateTrace(Efieldz[:,i],-peaktime[i]+time[ref_sample,i])
    Voltagex[:,i] = TranslateTrace(Voltagex[:,i],-peaktime[i]+timeV[ref_sample,i])
    Voltagey[:,i] = TranslateTrace(Voltagey[:,i],-peaktime[i]+timeV[ref_sample,i])
    Voltagez[:,i] = TranslateTrace(Voltagez[:,i],-peaktime[i]+timeV[ref_sample,i])

na_ref = Xprojrefn.size
na_cross = Xprojcrossn.size

Efieldx_ref = Efieldx[:,ou]
Efieldy_ref = Efieldy[:,ou]
Efieldz_ref = Efieldz[:,ou]

Voltagex_ref = Voltagex[:,ou]
Voltagey_ref = Voltagey[:,ou]
Voltagez_ref = Voltagez[:,ou]

Efieldx_cross = Efieldx[:,nou]
Efieldy_cross = Efieldy[:,nou]
Efieldz_cross = Efieldz[:,nou]

Voltagex_cross = Voltagex[:,nou]
Voltagey_cross = Voltagey[:,nou]
Voltagez_cross = Voltagez[:,nou]

# Display scatter plots of Y field trace for some time bins

# for t in range(365,375):

#     fig=plt.figure()
#     ax = fig.gca()
#     pl=ax.scatter(Xprojrefn,Yprojrefn,c=Efieldz[t,:na_ref])
#     fig.colorbar(pl,ax=ax)
#     plt.show()

# Now compute Zernike transform of field component traces per time bin


def TimeInterpolateTraces(tmin,tmax,nmax=20,mmax=3,regul=0):
    js = zernike.compute_j_list(nmax,mmax)
    ncoeff = len(js)
    Efieldx_coeffs = np.zeros((ntime,ncoeff))
    Efieldy_coeffs = np.zeros((ntime,ncoeff))
    Efieldz_coeffs = np.zeros((ntime,ncoeff))

    fitted_Efieldx_ref = np.zeros((ntime,na_ref))
    fitted_Efieldy_ref = np.zeros((ntime,na_ref))
    fitted_Efieldz_ref = np.zeros((ntime,na_ref))
    fitted_Efieldx_cross = np.zeros((ntime,na_cross))
    fitted_Efieldy_cross = np.zeros((ntime,na_cross))
    fitted_Efieldz_cross = np.zeros((ntime,na_cross))

    print (Xprojrefn.min(),Xprojrefn.max(),Yprojrefn.min(),Yprojrefn.max())

    for t in range(tmin,tmax):
        Efieldy_coeffs[t,:] = zernike.compute_zernike_coeffs(Efieldy[t,:na_ref],Xprojrefn,Yprojrefn,js,regul=regul)
        fitted_Efieldy_ref[t,:] = zernike.zernike_array_noll(Efieldy_coeffs[t,:],Xprojrefn,Yprojrefn,js)
        fitted_Efieldy_cross[t,:] = zernike.zernike_array_noll(Efieldy_coeffs[t,:],Xprojcrossn,Yprojcrossn,js)
    return fitted_Efieldy_ref, fitted_Efieldy_cross

def FourierInterpolateTraces(bandpass=None,nmax=20,mmax=3,regul=0,step=0.5):
    '''
    This routine computes modulus and phase of FFT of traces for all antennas, 
    then do a 2D interpolation of both modulus and phase to predict FFT of trace at 
    any given position. Finally performs inverse FFT to get back to time traces.
    bandpass, if not none, must give the min and max frequencies (in Hz)
    step is the time sampling in nano seconds.
    '''
    js = zernike.compute_j_list(nmax,mmax)
    ncoeff = len(js)
    # Create FFT arrays 
    nsampV = Voltagex.shape[0]
    nsampE = Efieldx.shape[0]
    FFT_Voltagex = np.zeros((nsampV,na_ref),dtype=np.complex128)
    FFT_Voltagey = np.zeros((nsampV,na_ref),dtype=np.complex128)
    FFT_Voltagez = np.zeros((nsampV,na_ref),dtype=np.complex128)
    freqsV = np.fft.fftfreq(nsampV,d=step*1e-9) # Frequency array in Hz

    # Begin by computing all reference antenna trace FFTs
    for i in range(na_ref):
        FFT_Voltagex[:,i] = np.fft.fft(Voltagex[:,i])
        FFT_Voltagey[:,i] = np.fft.fft(Voltagey[:,i])
        FFT_Voltagez[:,i] = np.fft.fft(Voltagez[:,i])

    # If bandpass is defined, compute frequency index
    if bandpass is not None:
        fmin = bandpass[0]
        fmax = bandpass[1]
        indmin = np.where(freqsV > fmin)[0][0]
        indmax = np.where(freqsV > fmax)[0][0]
    else:
        indmin = 0
        indmax = (np.where(freqsV<0))[0][0]

    # Define amplitude coefficient arrays 
    n_valid_freqsV = indmax-indmin
    FFT_Voltagex_amp_coeffs = np.zeros((nsampV,ncoeff))
    FFT_Voltagey_amp_coeffs = np.zeros((nsampV,ncoeff))
    FFT_Voltagez_amp_coeffs = np.zeros((nsampV,ncoeff))

    # Same for phase
    FFT_Voltagex_phi_coeffs = np.zeros((nsampV,ncoeff))
    FFT_Voltagey_phi_coeffs = np.zeros((nsampV,ncoeff))
    FFT_Voltagez_phi_coeffs = np.zeros((nsampV,ncoeff))


    return (FFT_Voltagex,FFT_Voltagey,FFT_Voltagez)

    
    # for ifreq in range(indmin,indmax):
    #    FFT_Voltagey_amp_coeffs[ifreq,:] = zernike.compute_zernike_coeffs(np.abs(FFT_Voltagey[ifreq,:na_ref]),Xprojrefn,Yprojrefn,js,regul=regul)

    # 
    # fig=plt.figure(figsize=(xfigsize,yfigsize))
    # ax=fig.add_subplot(121)
    # pl=ax.scatter(Xprojrefn,Yprojrefn,c=Efieldy[t,:na_ref])
    # plt.colorbar(pl,ax=ax)
    # axx=fig.add_subplot(122)
    # pll=axx.scatter(Xprojrefn,Yprojrefn,c=fitted_Efieldy_ref[t,:])
    # plt.colorbar(pll,ax=axx)
    # plt.show()

# fig=plt.figure()
# plt.plot(time[:,0],Efieldy[:,0])
# plt.plot(time[:,0],fitted_Efieldy_ref[:,0])
# plt.show()
