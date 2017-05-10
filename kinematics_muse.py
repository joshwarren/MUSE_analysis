## ==================================================================
## Finding the kinematic misalignments
## ==================================================================
## warrenj 20150825 Routine to find PA_phot, PA_kin and Psi, the
## misalignment angle. Also to find photometry info e.g. ellipticity. 


import numpy as np # for reading files
import os
import glob # for searching for files
from astropy.io import fits as pyfits # reads fits files (is from astropy)
from find_galaxy import find_galaxy # part of mge package, fits photometry
from fit_kinematic_pa import fit_kinematic_pa # fit kinemetry
import math # for sine functions
import matplotlib.pyplot as plt # used for plotting
import matplotlib.axes as ax # for adding text onto images
from scipy.optimize import curve_fit # for fitting a gaussian
from checkcomp import checkcomp
from classify import get_R_e
cc = checkcomp()
import cPickle as pickle

#---------------------------------------------------------------------------
def spxToKpc(x, z):
	H0 = 70.4 # Mpc / kms^-1
	val =  3*10**5 * z/H0 *10**3 * x * 0.67*4.85*10**(-6)
	return val

def spxToRe(x, R_e):
	val = x * 0.67 / R_e
	return val


def kinematics(galaxy, discard=0, plots=False, D=None):
	print '  kinematics'

	analysis_dir = "%s/Data/muse/analysis" % (cc.base_dir)
	galaxiesFile = "%s/galaxies.txt" % (analysis_dir)

	output = '%s/%s/results' % (analysis_dir, galaxy)
	if D is None:
		pickleFile = open('%s/pickled/dataObj_.pkl' % (output))
		D = pickle.load(pickleFile)
		pickleFile.close()


	x_gals, y_gals, SN_kin_gals, SN_pop_gals = np.loadtxt(
		galaxiesFile, unpack=True, skiprows=1, usecols=(1,2,3,4))
	galaxy_gals = np.loadtxt(galaxiesFile, skiprows=1, usecols=(0,),dtype=str)
	i_gal = np.where(galaxy_gals==galaxy)[0][0]

	R_e = get_R_e(galaxy)
# ------------=============== Photometry =================----------
	save_to = "%s/%s/results/" % (analysis_dir,	galaxy
		) + "plots/photometry_.png"
	f = find_galaxy(D.unbinned_flux, quiet=True, plot=plots)#, 
		#galaxy=galaxy.upper(), redshift=z, sav_fig=save_to)
	#f_err = find_galaxy(galaxy_data_error, quiet=True, plot=False)
	x_gals[i_gal] = f.xpeak # f.xmed?
	y_gals[i_gal] = f.ypeak # f.ymed?

	# print "ellip: " + str(f.eps) #+ "+/-" + str(abs(f.eps-f_err.eps))
	# print "PA_photo: " + str(90-f.theta) #+ "+/-" + str(abs(f.theta-f_err.theta))


# # ------------================ Lambda_R ==================----------
# 	# distance from center
# 	R = np.sqrt(np.square(xBar)+np.square(yBar))
# 	# Angle of r from axis of rotation
# 	#ang = abs(math.asin(math.sin(math.radians(f.theta))) - np.abs(
# 	#	np.arctan((xBar)/(yBar))))
# 	#R = r * np.sin(ang)

# 	order = np.argsort(R)

# 	xBar_r = -xBar*math.cos(f.theta)-yBar*math.sin(f.theta)
# 	yBar_r = xBar*math.sin(f.theta) - yBar*math.cos(f.theta)

# 	R_m = np.sqrt(np.square(xBar_r)*(1-f.eps) + np.square(yBar_r)/(1+f.eps))
# 	h = 0.1 # discretization accuracy
# 		# discretizing the unit square
# 	x, y = np.mgrid[-f.xmed:40-2*discard-f.xmed:h, -f.ymed:40-2*discard-f.ymed:h]
# 	for i, r in enumerate(R_m):
# 			# set all points of ellipse that are inside 
# 		el = (np.square(-x*math.cos(f.theta)-y*math.sin(f.theta))*(1-f.eps) + \
# 			np.square(x*math.sin(f.theta) - y*math.cos(f.theta))/(1+f.eps))/r**2 <= 1
		
# 		if el.any():
# 			A_s = np.sum(el) * h * h
# 			A_ellipse = math.pi * r**2 * math.sqrt((1+f.eps**2)/(1-f.eps))
# 			if 0.85 * A_ellipse < A_s:
# 				R_m[i] = math.sqrt(A_s/math.pi)
# 			else:
# 				R_m[i] = np.nan

# 	order_m = np.argsort(R_m)

# 	# NB: lam is ordered in terms of increasing R.
# 	lam_num = D.flux[order_m]*R[order_m]*np.abs(np.array(
# 		D.components['stellar'].plot['vel'][order_m])) # numerator
# 	lam_den = D.flux[order_m]*R[order_m]*np.sqrt(np.square(
# 		np.array(D.components['stellar'].plot['vel'])[order_m]) + np.square(
# 		np.array(D.components['stellar'].plot['sigma'])[order_m])) # denominator

# 	# cumulative summation with mask from R_m
# 	lam_num = np.cumsum(lam_num[~np.isnan(R_m)])
# 	lam_den = np.cumsum(lam_den[~np.isnan(R_m)])

# 	lam = lam_num/lam_den
# 	plt.figure()
# 	plt.title(r"Radial $\lambda_R$ profile")
# 	plt.ylabel(r"$\lambda_R$")
# 	plt.xlabel("Radius (R_e)")
# 	x = spxToRe(R[order_m], R_e)[~np.isnan(R_m[order_m])] # Plotted as a fucntion of R not R_m
# 	order = np.argsort(x)
# 	lambda_R[i_gal2] = lam[order][-1]
# 	print 'lambda_R_MAX: ', lambda_R[i_gal2]
# 	plt.plot(x[order[5:]], lam[order[5:]])
# 	ax =plt.gca()
# 	plt.text(0.02,0.98, "Galaxy: " + galaxy.upper(), verticalalignment='top',
# 		transform=ax.transAxes)
# 	plt.savefig("%s/%s/results/%s/plots/" % (analysis_dir, galaxy, wav_range_dir
# 		) + "lambda_R_%s.png" % (wav_range), bbox_inches="tight")
# 	if plots: 
# 		plt.show()

# ------------============== Save outputs ================----------
	template = "{0:12}{1:4}{2:4}{3:8}{4:8}\n"

	with open(galaxiesFile, 'wb') as f:
		f.write(template.format("Galaxy", "x", "y", "Kin SN", "Pop SN"))

		for i in range(len(galaxy_gals)):
			f.write(template.format(galaxy_gals[i], str(int(x_gals[i])), 
				str(int(y_gals[i])), str(SN_kin_gals[i]), str(SN_pop_gals[i])))

##############################################################################

# Use of kinematics.py

if __name__ == '__main__':
	galaxy = 'eso443-g024'
	discard = 2
	wav_range = '4200-'

	kinematics(galaxy, discard=discard, wav_range=wav_range, plots=False)
