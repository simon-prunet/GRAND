import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from Simulations import hdf5fileinout as hdf5io
from Zernike import zernike
nmax=9
mmax=3
regul=0
azim=-66
elev=10
xfigsize=16
yfigsize=10
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
Xscale = 1.01*(Xref.max()-Xref.min())/2.
Xrefn = Xref / Xscale
Yscale = 1.01*(Yref.max()-Yref.min())/2.
Yrefn = Yref / Yscale

peaktimeref = peaktime[ou]
peakamplituderef = peakamplitude[ou]

# Compute list of js from nmax, mmax
js = zernike.compute_j_list(nmax,mmax)
# Compute Zernike coefficients
peaktime_coeffs = zernike.compute_zernike_coeffs(peaktimeref,Xrefn,Yrefn,js,regul=regul)
peakamplitude_coeffs = zernike.compute_zernike_coeffs(peakamplituderef,Xrefn,Yrefn,js,regul=regul)

fitted_peaktimeref = zernike.zernike_array_noll(peaktime_coeffs,Xrefn,Yrefn,js)
fitted_peakamplituderef = zernike.zernike_array_noll(peakamplitude_coeffs,Xrefn,Yrefn,js)

# Fit quality on reference antennas

fig=plt.figure(figsize=(xfigsize,yfigsize))
fig.canvas.set_window_title('Peak time on star shape')
ax=fig.add_subplot(121,projection='3d',azim=azim,elev=elev)
ax.scatter3D(Xref,Yref,peaktimeref)
ax.set_title('Original peak time')
ax2=fig.add_subplot(122,projection='3d',azim=azim,elev=elev)
ax2.scatter3D(Xref,Yref,peaktimeref-fitted_peaktimeref)
ax2.set_title('Residual error in peak time')
fig.tight_layout()
fig.savefig('Peak_time_star_shape.png',dpi=300)
plt.show()


fig2=plt.figure(figsize=(xfigsize,yfigsize))
fig2.canvas.set_window_title('Peak amplitude on star shape')
axx=fig2.add_subplot(131,projection='3d',azim=azim,elev=elev)
axx.scatter3D(Xref,Yref,peakamplituderef)
axx.set_title('Original peak amplitude')
axx2=fig2.add_subplot(132,projection='3d',azim=azim,elev=elev)
axx2.scatter3D(Xref,Yref,fitted_peakamplituderef)
axx2.set_title('Fitted peak amplitude')
axx3=fig2.add_subplot(133,projection='3d',azim=azim,elev=elev)
axx3.scatter3D(Xref,Yref,peakamplituderef-fitted_peakamplituderef)
axx3.set_title('Residual error in peak amplitude')
fig2.tight_layout()
fig2.savefig('Peak_amplitude_star_shape.png',dpi=300)
plt.show()

# Fit quality on crosscheck antennas
nou=[not k for k in ou]
Xcr = X[nou]; Xcr=Xcr.data
Ycr = Y[nou]; Ycr=Ycr.data
# Normalize coordinates
Xcrn = Xcr / Xscale
Ycrn = Ycr / Yscale
# Compute predictions from Zernike models
peaktimecross = peaktime[nou]
peakamplitudecross = peakamplitude[nou]

pred_peaktimecross = zernike.zernike_array_noll(peaktime_coeffs,Xcrn,Ycrn,js)
pred_peakamplitudecross = zernike.zernike_array_noll(peakamplitude_coeffs,Xcrn,Ycrn,js)

fig3=plt.figure(figsize=(xfigsize,yfigsize))
fig3.canvas.set_window_title('Peak time on test antennas')
axxx=fig3.add_subplot(121,projection='3d',azim=azim,elev=elev)
axxx.scatter3D(Xcr,Ycr,peaktimecross)
axxx.set_title('Original peak time')
axxx2=fig3.add_subplot(122,projection='3d',azim=azim,elev=elev)
axxx2.scatter3D(Xcr,Ycr,peaktimecross-pred_peaktimecross)
axxx2.set_title('Residual peak time error')
fig3.tight_layout()
fig3.savefig('Peak_time_test_antennas.png',dpi=300)
plt.show()

fig4=plt.figure(figsize=(xfigsize,yfigsize))
fig4.canvas.set_window_title('Peak amplitude on test antennas')
axxxx=fig4.add_subplot(121,projection='3d',azim=azim,elev=elev)
axxxx.scatter3D(Xcr,Ycr,peakamplitudecross)
axxxx.set_title('Original peak amplitudes')
axxxx2=fig4.add_subplot(122,projection='3d',azim=azim,elev=elev)
axxxx2.scatter3D(Xcr,Ycr,peakamplitudecross-pred_peakamplitudecross)
axxxx2.set_title('Residual peak amplitude error')
fig4.tight_layout()
fig4.savefig('Peak_amplitude_test_antennas.png',dpi=300)
plt.show()

fig5=plt.figure(figsize=(xfigsize/1.5,yfigsize))
ay=fig5.add_subplot(111)
ay.set_title('Peak amplitude relative residual error on star shape')
rel_error = (peakamplituderef-fitted_peakamplituderef)/(fitted_peakamplituderef)
pl=ay.scatter(Xref,Yref,c=rel_error)
fig5.colorbar(pl,ax=ay)
fig5.tight_layout()
fig5.savefig('Peak_amplitude_relative_error_star_shape.png',dpi=300)
plt.show()

fig6=plt.figure(figsize=(xfigsize/1.5,yfigsize))
ayy=fig6.add_subplot(111)
ayy.set_title('Peak time residual error on star shape')
error = (peaktimeref-fitted_peaktimeref)
pl=ayy.scatter(Xref,Yref,c=error)
fig6.colorbar(pl,ax=ayy)
fig6.tight_layout()
fig6.savefig('Peak_time_error_star_shape.png',dpi=300)
plt.show()

fig7=plt.figure(figsize=(xfigsize/1.5,yfigsize))
ayyy=fig7.add_subplot(111)
ayyy.set_title('Peak amplitude relative residual error on test antennas')
rel_error = (peakamplitudecross-pred_peakamplitudecross)/(pred_peakamplitudecross)
pl=ayyy.scatter(Xcr,Ycr,c=rel_error)
fig7.colorbar(pl,ax=ayyy)
fig7.tight_layout()
fig7.savefig('Peak_amplitude_relative_error_test_antennas.png',dpi=300)
plt.show()

fig8=plt.figure(figsize=(xfigsize/1.5,yfigsize))
ayyyy=fig8.add_subplot(111)
ayyyy.set_title('Peak time residual error on test antennas')
error = (peaktimecross-pred_peaktimecross)
pl=ayyyy.scatter(Xcr,Ycr,c=error)
fig8.colorbar(pl,ax=ayyyy)
fig8.tight_layout()
fig8.savefig('Peak_time_error_test_antennas.png',dpi=300)
plt.show()


