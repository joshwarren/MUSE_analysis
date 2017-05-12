## ==================================================================
## Finding the kinematic misalignments
## ==================================================================
## warrenj 20150825 Routine to find PA_phot, PA_kin and Psi, the
## misalignment angle. Also to find photometry info e.g. ellipticity. 
from checkcomp import checkcomp
cc = checkcomp()
if 'home' not in cc.device:
	import matplotlib # 20160202 JP to stop lack-of X-windows error
	matplotlib.use('Agg') 
import numpy as np # for reading files
import os
from find_galaxy import find_galaxy # part of mge package, fits photometry
from fit_kinematic_pa import fit_kinematic_pa # fit kinemetry
import matplotlib.pyplot as plt # used for plotting
import matplotlib.axes as ax # for adding text onto images
from scipy.optimize import curve_fit # for fitting a gaussian
from scipy.interpolate import interp1d
from classify import get_R_e
import cPickle as pickle

#---------------------------------------------------------------------------
def spxToKpc(x, z):
	H0 = 70.4 # Mpc / kms^-1
	val =  3*10**5 * z/H0 *10**3 * x * 0.67*4.85*10**(-6)
	return val

def spxToRe(x, R_e):
	val = x * 0.67 / R_e
	return val


def kinematics(galaxy, opt='kin', discard=0, plots=False, D=None):
	print '  kinematics'

	analysis_dir = "%s/Data/muse/analysis" % (cc.base_dir)
	galaxiesFile = "%s/galaxies.txt" % (analysis_dir)

	output = '%s/%s/%s' % (analysis_dir, galaxy, opt)
	if D is None:
		pickleFile = open('%s/pickled/dataObj.pkl' % (output))
		D = pickle.load(pickleFile)
		pickleFile.close()
	
	d = np.loadtxt(data_file, unpack=True, dtype=str)
	galaxy_gals = d[0][1:]
	x_gals, y_gals = d[1][1:].astype(int), d[2][1:].astype(int)
	SN_gals = {d[i][0]:d[i][1:].astype(float) for i in range(3,len(d))}
	i_gal = np.where(galaxy_gals==galaxy)[0][0]

	R_e = get_R_e(galaxy)
# ------------=============== Photometry =================----------
	# save_to = "%s/%s/results/" % (analysis_dir,	galaxy
	# 	) + "plots/photometry_.png"
	f = find_galaxy(D.unbinned_flux, quiet=True, plot=plots)#, 
		#galaxy=galaxy.upper(), redshift=z, sav_fig=save_to)
	#f_err = find_galaxy(galaxy_data_error, quiet=True, plot=False)
	x_gals[i_gal] = f.xpeak # f.xmed?
	y_gals[i_gal] = f.ypeak # f.ymed?

	# print "ellip: " + str(f.eps) #+ "+/-" + str(abs(f.eps-f_err.eps))
	# print "PA_photo: " + str(90-f.theta) #+ "+/-" + str(abs(f.theta-f_err.theta))


# ------------================ Lambda_R ==================----------

	beta = np.tan(D.yBar/D.xBar) # Angle between point and center of galaxy and RA 
									# axis.
	pa = 90 - f.theta # Position angle
	res = 0.2 # arcsec. Spatial resolution of MUSE

	# Semi-major axis
	a = res* np.sqrt(((D.xBar - f.xpeak)**2 + (D.yBar - f.ypeak)**2) * \
		(np.sin(beta + pa)**2 + np.cos(beta + pa)**2/(1 - f.eps)**2))

	R_m = a * np.sqrt(1 - f.eps)
	R_m_sort = R_m.argsort() # argsort of R_m_sort will be the inverse opperation 
	
	# Area of ellipse
	A_ellipse = np.pi * a**2 * (1 - f.eps)
	# Area sampled
	A_s = np.cumsum(D.n_spaxels_in_bin[R_m_sort])[R_m_sort.argsort()] * res**2

	R_m[A_ellipse > A_s] = np.sqrt(A_s[A_ellipse>A_s]/np.pi)
	# R_max occurs when 0.85*A_ellipse = A_s
	# R_m[0.85*A_ellipse > A_s] = np.nan

	# NB: numerator and denominator are in R_m order
	numerator = np.cumsum(D.flux[R_m_sort] * 
		np.sqrt(D.xBar**2 + D.yBar**2)[R_m_sort] * 
		np.abs(D.components['stellar'].plot['vel'][R_m_sort]))

	denominator = np.cumsum(D.flux[R_m_sort] * 
		np.sqrt(D.xBar**2 + D.yBar**2)[R_m_sort] * 
		np.sqrt(D.components['stellar'].plot['vel']**2 + 
			D.components['stellar'].plot['sigma']**2)[R_m_sort])

	lambda_R = numerator[R_m_sort.argsort()]/denominator[R_m_sort.argsort()]
	lambda_R[np.isnan(R_m)] = np.nan

	lambda_Re = interp1d(R_m[~np.isnan(R_m)], lambda_R[~np.isnan(R_m)],
		bounds_error=False, fill_value=(np.nan, 
			lambda_R[R_m_sort][~np.isnan(lambda_R[R_m_sort])][-1]))(R_e)

	print 'lambda_Re: ', lambda_Re

	fig, ax = plt.subplots()
	ax.set_title(r"Radial $\lambda_R$ profile")
	ax.set_ylabel(r"$\lambda_R$")
	ax.set_xlabel(r"Radius ($R_m$/$R_e$)")
	ax.plot(R_m[R_m_sort][10:]/R_e, lambda_R[R_m_sort][10:])
	# ax.text(0.02,0.98, "Galaxy: " + galaxy.upper())#, verticalalignment='top',
		# transform=ax.transAxes)
	plt.savefig("%s/plots/lambda_R.png" % (output), bbox_inches="tight")
	if plots: 
		plt.show()
# ------------============== Save outputs ================----------
	temp = "{0:12}{1:4}{2:4}"+''.join(['{%i:%i}'%(i+3,len(t)+1) for i, t in 
		enumerate(SN_gals.keys())])+'\n'

	SN_titles = list(SN_gals.keys())
	with open(data_file, 'w') as f:
		f.write(temp.format("Galaxy", "x", "y", *(s for s in SN_titles)))
		for i in range(len(galaxy_gals)):
			f.write(temp.format(galaxy_gals[i], str(int(x_gals[i])), 
				str(int(y_gals[i])), *(str(round(SN_gals[s][i],2)) for s in SN_titles)))

##############################################################################

# Use of kinematics.py

if __name__ == '__main__':
	galaxy = 'ic1459'
	discard = 0

	kinematics(galaxy, opt='kin', discard=discard, plots=False)