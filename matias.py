import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from Simulations import hdf5fileinout as hdf5io
from Starshape.geom import LinePlaneCollision, Simplerot
from Zernike import zernike
import numpy as np


# nmax_ref = 20
# mmax_ref = 4

class Starshape:

    def __init__(self,simu_path='/Users/prunet/Documents/Code/GRAND/Starshape-Zernike/Stshp_Proton_3.98_77.4_0.0_1.hdf5',normalized_coordinates_path='/Users/prunet/Documents/Code/GRAND/Starshape-Zernike/TestInput.inpnorm'):


        # Reference sample
        self.ref_sample = 0

        #Shower path
        self.simu_path = simu_path
        self.normalized_coordinates_path = normalized_coordinates_path

        InputFilename = simu_path

        self.EventNumber = 0 # Only one simulation here
        #Run info
        self.RunInfo = hdf5io.GetRunInfo(InputFilename)
        print("#RunInfo --> ", self.RunInfo)
        self.EventName = hdf5io.GetEventName(RunInfo,EventNumber)
        print("#EventName --> ", self.EventName)

        #Shower parameters
        self.Zenith = hdf5io.GetEventZenith(RunInfo,EventNumber)
        print("#Zenith --> ", self.Zenith)
        #Azimuth = hdf5io.GetEventAzimuth(RunInfo,EventNumber)
        #print("#Azimuth --> ", Azimuth)
        self.XmaxDistance = hdf5io.GetEventXmaxDistance(RunInfo,EventNumber)
        print("#XmaxDistance --> ", self.XmaxDistance)

        #Antannas info
        self.AntennaInfo = hdf5io.GetAntennaInfo(InputFilename,EventName)
        print("#AntennaInfo --> ", self.AntennaInfo)
        self.IDs = self.AntennaInfo['ID'].data
        print("#IDs --> ", self.IDs)

        #Get filter on IDs name, to get index of antennas in the star shape array
        ou=[id.startswith('A') for id in IDs]

        #Get X,Y coordinates of all antennas
        self.X = hdf5io.GetXFromAntennaInfo(AntennaInfo)
        self.Y = hdf5io.GetYFromAntennaInfo(AntennaInfo)
        self.na_tot = self.X.size

        #Get peak Hilbert time and amplitude
        self.peaktime, self.peakamplitude = hdf5io.get_peak_time_hilbert_hdf5(InputFilename, antennamax="All",antennamin=0, usetrace="efield", DISPLAY=False)

        # Select data from star shape
        Xref = X[ou]; self.Xref=Xref.data
        Yref = Y[ou]; self.Yref=Yref.data

        # Crosscheck antennas
        nou=[not k for k in ou]
        Xcross = X[nou]; self.Xcross=Xcross.data
        Ycross = Y[nou]; self.Ycross=Ycross.data

        # Compute predictions from Zernike models
        self.peaktimecross = peaktime[nou]
        self.peakamplitudecross = peakamplitude[nou]

        # # Now compute coordinates of antennas reprojected into plane perpendicular to shower axis and passing by origin
        # el_angle = 180.-Zenith
        # mat = Simplerot(el_angle) # Uses elevation angle
        # na_tot = X.size
        # Xproj = np.zeros(na_tot)
        # Yproj = np.zeros(na_tot)
        # origin = np.zeros(3)
        # planeNormal = np.array([np.sin(np.deg2rad(el_angle)),0.0,np.cos(np.deg2rad(el_angle))])
        # xmax_pos = XmaxDistance * planeNormal # 3D position of Xmax
        # for i in range(na_tot):
        #     # 3D position of antenna
        #     antenna_pos = np.array([X[i],Y[i],0.0])
        #     # Unit vector along light ray from antenna to Xmax
        #     rayDirection = xmax_pos - antenna_pos; rayDirection /= np.linalg.norm(rayDirection)
        #     # Compute intersection of light ray and plane perpendicular to shower axis
        #     intersection_pos = LinePlaneCollision(planeNormal,origin,rayDirection,antenna_pos)
        #     # Now rotate it so that X,Y coordinates lie in the sample plane perpendicular to shower axis
        #     rotated_pos = np.dot(mat,intersection_pos)
        #     Xproj[i] = rotated_pos[0]
        #     Yproj[i] = rotated_pos[1]
        #     # print i,rotated_pos[2]

        # Get normalized positions from separate file
        arr = np.loadtxt(normalized_coordinates_path)
        self.Xprojn = arr[:,0]
        self.Yprojn = arr[:,1]
        self.Xprojrefn = Xprojn[ou]
        self.Yprojrefn = Yprojn[ou]
        self.Xprojcrossn = Xprojn[nou]
        self.Yprojcrossn = Yprojn[nou]

        # Show scatter plot of ground coordinates of antennas
        fig=plt.figure()
        ax=fig.gca()
        ax.set_title('Original (ground) antenna positions')
        ax.set_aspect('equal')
        ax.scatter(self.X,self.Y)
        plt.show()

        # Show scatter plot of reprojected antenna coordinates in plane perpendicular to shower axis
        fig = plt.figure()
        ax = fig.gca()
        ax.set_title('Reprojected antenna positions (normalized)')
        ax.set_aspect('equal')
        ax.scatter(self.Xprojn,self.Yprojn)
        plt.show()

        # Compute numbers of ref and cross antennas after filtering on radius: na_ref, na_cross
        self.na_ref = self.Xprojrefn.size
        self.na_cross = self.Xprojcrossn.size
        # Update ou, nou masks to include filtering on radius. They will be defined as list of valid indices
        # instead of boolean mask
        # ou = (np.flatnonzero(ou))[radrefn_mask]
        # nou = (np.flatnonzero(nou))[radcrossn_mask]


        # Now extract the Efield traces, and store them per component as an ntime,nantennas vector
        # First get traces for one antenna to determine the number of time bins
        tmp = hdf5io.GetAntennaEfield(InputFilename,EventName,IDs[0])
        self.ntime = tmp.shape[0]
        tmp = hdf5io.GetAntennaFilteredVoltage(InputFilename,EventName,IDs[0])
        self.ntimeV = tmp.shape[0]

        # Efield and time arrays
        self.Efieldx = np.zeros((ntime,na_tot))
        self.Efieldy = np.zeros((ntime,na_tot))
        self.Efieldz = np.zeros((ntime,na_tot))
        self.time    = np.zeros((ntime,na_tot))

        # Voltage and time arrays
        self.Voltagex = np.zeros((ntimeV,na_tot))
        self.Voltagey = np.zeros((ntimeV,na_tot))
        self.Voltagez = np.zeros((ntimeV,na_tot))
        self.timeV    = np.zeros((ntimeV,na_tot))

        # Displacement vectors
        #self.Efieldx_dt = np.zeros(na_tot)
        #self.Efieldy_dt = np.zeros(na_tot)
        #self.Efieldz_dt = np.zeros(na_tot)

        for i in range(self.na_tot):

            # E field
            trace = hdf5io.GetAntennaEfield(InputFilename,EventName,IDs[i])
            self.time[:,i] = trace[:,0]
            self.Efieldx[:,i] = trace[:,1]
            self.Efieldy[:,i] = trace[:,2]
            self.Efieldz[:,i] = trace[:,3]
            # Antenna voltages
            traceV = hdf5io.GetAntennaFilteredVoltage(InputFilename,EventName,IDs[i])
            self.timeV[:,i]    = traceV[:,0]
            self.Voltagex[:,i] = traceV[:,1]
            self.Voltagey[:,i] = traceV[:,2]
            self.Voltagez[:,i] = traceV[:,3]

            # Use Hilbert envelop peak time and FFT to translate all traces such that their peaks are aligned
            # Efieldx_dt[i] = ComputeDisplacement(Efieldx[:,0],Efieldx[:,i])
            # Efieldx[:,i] = TranslateTrace(Efieldx[:,i],Efieldx_dt[i])
            # Efieldy_dt[i] = ComputeDisplacement(Efieldy[:,0],Efieldy[:,i])
            # Efieldy[:,i] = TranslateTrace(Efieldy[:,i],Efieldy_dt[i])
            # Efieldz_dt[i] = ComputeDisplacement(Efieldz[:,0],Efieldz[:,i])
            # Efieldz[:,i] = TranslateTrace(Efieldz[:,i],Efieldz_dt[i])
            # Efieldx[:,i] = TranslateTrace(Efieldx[:,i],-peaktime[i]+time[ref_sample,i])
            # Efieldy[:,i] = TranslateTrace(Efieldy[:,i],-peaktime[i]+time[ref_sample,i])
            # Efieldz[:,i] = TranslateTrace(Efieldz[:,i],-peaktime[i]+time[ref_sample,i])
            self.Voltagex[:,i] = self.TranslateTrace(self.Voltagex[:,i],-self.peaktime[i]+self.timeV[self.ref_sample,i])
            self.Voltagey[:,i] = self.TranslateTrace(self.Voltagey[:,i],-self.peaktime[i]+self.timeV[self.ref_sample,i])
            self.Voltagez[:,i] = self.TranslateTrace(self.Voltagez[:,i],-self.peaktime[i]+self.timeV[self.ref_sample,i])
            #Voltagey[:,i] = TranslateTrace(Voltagey[:,i],-peaktime[i]+timeV[ref_sample,i])


        self.Efieldx_ref = self.Efieldx[:,ou]
        self.Efieldy_ref = self.Efieldy[:,ou]
        self.Efieldz_ref = self.Efieldz[:,ou]

        self.Voltagex_ref = self.Voltagex[:,ou]
        self.Voltagey_ref = self.Voltagey[:,ou]
        self.Voltagez_ref = self.Voltagez[:,ou]

        self.Efieldx_cross = self.Efieldx[:,nou]
        self.Efieldy_cross = self.Efieldy[:,nou]
        self.Efieldz_cross = self.Efieldz[:,nou]

        self.Voltagex_cross = self.Voltagex[:,nou]
        self.Voltagey_cross = self.Voltagey[:,nou]
        self.Voltagez_cross = self.Voltagez[:,nou]

################################################################################
# Test of Matias functions for simulation and hdf5 files handling (Valentin)
################################################################################
    def TranslateTrace(self,trace,deltat,step=0.5):

        '''
        Use FFT to translate a trace by deltat seconds. Sampling rate is given by step in ns.
        '''

        nsamp = trace.size
        freq = np.fft.fftfreq(nsamp,d=step)
        ftrace = np.fft.fft(trace)
        phase = np.exp(-1j * 2.0*np.pi * freq * deltat)
        return np.fft.ifft(ftrace*phase)

    def ComputeDisplacement(self,trace_ref, trace, step=0.5):
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

    def TimeInterpolateTraces(self,tmin,tmax,nmax=20,mmax=3,regul=0):
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

    def PrepareFourierInterpolateTraces(self,bandpass=None,nmax_amp=20,nmax_phi=3,mmax_amp=4,mmax_phi=3,regul_amp=0.,step=0.5):
        '''
        This routine computes modulus and phase of FFT of traces for all antennas, 
        then do a 2D interpolation of both modulus and phase to predict FFT of trace at 
        any given position. Finally performs inverse FFT to get back to time traces.
        bandpass, if not none, must give the min and max frequencies (in Hz)
        step is the time sampling in nano seconds.
        '''
        self.bandpass = bandpass
        self.nmax_amp = nmax_amp
        self.nmax_phi = nmax_phi
        self.mmax_amp = mmax_amp
        self.mmax_phi = mmax_phi
        self.regul_amp = regul_amp
        self.step = step

        self.js_amp = zernike.compute_j_list(self.nmax_amp,self.mmax_amp)
        self.js_phi = zernike.compute_j_list(self.nmax_phi,self.mmax_phi)
        self.ncoeff_amp = len(self.js_amp)
        self.ncoeff_phi = len(self.js_phi)

        # Create FFT arrays 
        self.nsampV = self.Voltagex.shape[0]
        self.nsampE = self.Efieldx.shape[0]
        self.FFT_Voltagex_ref = np.zeros((self.nsampV,self.na_ref),dtype=np.complex128)
        self.FFT_Voltagey_ref = np.zeros((self.nsampV,self.na_ref),dtype=np.complex128)
        self.FFT_Voltagez_ref = np.zeros((self.nsampV,self.na_ref),dtype=np.complex128)
        self.FFT_Voltagex_cross = np.zeros((self.nsampV,self.na_cross),dtype=np.complex128)
        self.FFT_Voltagey_cross = np.zeros((self.nsampV,self.na_cross),dtype=np.complex128)
        self.FFT_Voltagez_cross = np.zeros((self.nsampV,self.na_cross),dtype=np.complex128)

        self.freqsV = np.fft.fftfreq(self.nsampV,d=self.step*1e-9) # Frequency array in Hz

        # Begin by computing all reference antenna trace FFTs
        for i in np.arange(self.na_ref):
            self.FFT_Voltagex_ref[:,i] = np.fft.fft(self.Voltagex_ref[:,i])
            self.FFT_Voltagey_ref[:,i] = np.fft.fft(self.Voltagey_ref[:,i])
            self.FFT_Voltagez_ref[:,i] = np.fft.fft(self.Voltagez_ref[:,i])
        for i in np.arange(self.na_cross):
            self.FFT_Voltagex_cross[:,i] = np.fft.fft(self.Voltagex_cross[:,i])
            self.FFT_Voltagey_cross[:,i] = np.fft.fft(self.Voltagey_cross[:,i])
            self.FFT_Voltagez_cross[:,i] = np.fft.fft(self.Voltagez_cross[:,i])

        # If bandpass is defined, compute frequency index
        if self.bandpass is not None:
            self.fmin = self.bandpass[0]
            self.fmax = self.bandpass[1]
            self.indmin = np.where(self.freqsV > self.fmin)[0][0]
            self.indmax = np.where(self.freqsV > self.fmax)[0][0]
        else:
            self.indmin = 0
            self.indmax = (np.where(self.freqsV<0))[0][0]

        # Define amplitude coefficient arrays 
        self.n_valid_freqsV = self.indmax-self.indmin
        self.FFT_Voltagex_ref_amp_coeffs = np.zeros((self.nsampV,self.ncoeff_amp))
        self.FFT_Voltagey_ref_amp_coeffs = np.zeros((self.nsampV,self.ncoeff_amp))
        self.FFT_Voltagez_ref_amp_coeffs = np.zeros((self.nsampV,self.ncoeff_amp))

        ##### TO BE CONTINUED ##### 


        # Same for phase
        FFT_Voltagex_ref_phi_coeffs = np.zeros((nsampV,ncoeff_phi))
        FFT_Voltagey_ref_phi_coeffs = np.zeros((nsampV,ncoeff_phi))
        FFT_Voltagez_ref_phi_coeffs = np.zeros((nsampV,ncoeff_phi))

        # For all valid (positive) frequencies, compute Zernike coefficients for amplitude and phase
        # Amplitude
        for i in np.arange(indmin,indmax):
            print("Processing amplitude, frequency bin {0:d}\n",format(i))
            if (i==0 or i==nsampV/2): # Zero mode and Nyquist are real
                FFT_Voltagex_ref_amp_coeffs[i,:] = zernike.compute_zernike_coeffs(FFT_Voltagex_ref[i,:],Xprojrefn,Yprojrefn,js_amp,regul=regul_amp)
                FFT_Voltagey_ref_amp_coeffs[i,:] = zernike.compute_zernike_coeffs(FFT_Voltagey_ref[i,:],Xprojrefn,Yprojrefn,js_amp,regul=regul_amp)
                FFT_Voltagez_ref_amp_coeffs[i,:] = zernike.compute_zernike_coeffs(FFT_Voltagez_ref[i,:],Xprojrefn,Yprojrefn,js_amp,regul=regul_amp)
            else:
                FFT_Voltagex_ref_amp_coeffs[i,:] = zernike.compute_zernike_coeffs(np.abs(FFT_Voltagex_ref[i,:]),Xprojrefn,Yprojrefn,js_amp,regul=regul_amp)
                FFT_Voltagey_ref_amp_coeffs[i,:] = zernike.compute_zernike_coeffs(np.abs(FFT_Voltagey_ref[i,:]),Xprojrefn,Yprojrefn,js_amp,regul=regul_amp)
                FFT_Voltagez_ref_amp_coeffs[i,:] = zernike.compute_zernike_coeffs(np.abs(FFT_Voltagez_ref[i,:]),Xprojrefn,Yprojrefn,js_amp,regul=regul_amp)
            # Copy over negative frequencies
            if (i>0 and i<nsampV/2):
                FFT_Voltagex_ref_amp_coeffs[-i,:] = FFT_Voltagex_ref_amp_coeffs[i,:]
                FFT_Voltagey_ref_amp_coeffs[-i,:] = FFT_Voltagey_ref_amp_coeffs[i,:]
                FFT_Voltagez_ref_amp_coeffs[-i,:] = FFT_Voltagez_ref_amp_coeffs[i,:]

        # Phase. Use weighting scheme in fit according to amplitude (phase gets bogus when amplitude is too low at large radii)
        for i in np.arange(indmin,indmax):
            print("Processing phase, frequency bin {0:d}\n".format(i))
            if (i>0 and i<nsampV/2):
                FFT_Voltagex_ref_phi_coeffs[i,:] = zernike.compute_zernike_coeffs(np.unwrap(np.angle(FFT_Voltagex_ref[i,:])),Xprojrefn,Yprojrefn,js_phi,weights=np.abs(FFT_Voltagex_ref[i,:]))
                FFT_Voltagey_ref_phi_coeffs[i,:] = zernike.compute_zernike_coeffs(np.unwrap(np.angle(FFT_Voltagey_ref[i,:])),Xprojrefn,Yprojrefn,js_phi,weights=np.abs(FFT_Voltagey_ref[i,:]))
                FFT_Voltagez_ref_phi_coeffs[i,:] = zernike.compute_zernike_coeffs(np.unwrap(np.angle(FFT_Voltagez_ref[i,:])),Xprojrefn,Yprojrefn,js_phi,weights=np.abs(FFT_Voltagez_ref[i,:]))
                # Copy over negative frequencies. Phase gets a minus sign so that the Fourier coefficient gets complex conjugated.
                FFT_Voltagex_ref_phi_coeffs[-i,:] = -FFT_Voltagex_ref_phi_coeffs[i,:]
                FFT_Voltagey_ref_phi_coeffs[-i,:] = -FFT_Voltagey_ref_phi_coeffs[i,:]
                FFT_Voltagez_ref_phi_coeffs[-i,:] = -FFT_Voltagez_ref_phi_coeffs[i,:]

        return (FFT_Voltagex_ref,FFT_Voltagey_ref,FFT_Voltagez_ref,FFT_Voltagex_ref_amp_coeffs,FFT_Voltagey_ref_amp_coeffs,FFT_Voltagez_ref_amp_coeffs,FFT_Voltagex_ref_phi_coeffs,FFT_Voltagey_ref_phi_coeffs,FFT_Voltagez_ref_phi_coeffs)

    
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

def sample_disc(n=1000):

    u = np.random.uniform(size=n)
    r = np.sqrt(u)
    p = np.random.uniform(high=2.*np.pi,size=n)
    return (r*np.cos(p), r*np.sin(p))

# fig=plt.figure()
# plt.plot(time[:,0],Efieldy[:,0])
# plt.plot(time[:,0],fitted_Efieldy_ref[:,0])
# plt.show()
