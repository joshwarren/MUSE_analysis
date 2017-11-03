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
from scipy.interpolate import interp1d
from classify import get_R_e
import cPickle as pickle
from plot_results_muse import set_lims

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
	galaxiesFile2 = "%s/galaxies2.txt" % (analysis_dir)

	output = '%s/%s/%s' % (analysis_dir, galaxy, opt)
	if D is None:
		pickleFile = open('%s/pickled/dataObj.pkl' % (output))
		D = pickle.load(pickleFile)
		pickleFile.close()
	
	d = np.loadtxt(galaxiesFile, unpack=True, dtype=str)
	galaxy_gals = d[0][1:]
	x_gals, y_gals = d[1][1:].astype(int), d[2][1:].astype(int)
	SN_gals = {d[i][0]:d[i][1:].astype(float) for i in range(3,len(d))}
	i_gal = np.where(galaxy_gals==galaxy)[0][0]

	galaxy_gals2 = np.loadtxt(galaxiesFile2, unpack=True, usecols=(0,), dtype=str, 
		skiprows=1)
	lambda_Re_gals, ellipticity_gals, pa_gals, star_kine_pa_gals, gas_kine_pa_gals  =  \
		np.loadtxt(galaxiesFile2, unpack=True, skiprows=1, usecols=(1,2,3,4,5))
	i_gal2 = np.where(galaxy_gals2==galaxy)[0][0]

	R_e = get_R_e(galaxy)
# ------------=============== Photometry =================----------
	# save_to = "%s/%s/results/" % (analysis_dir,	galaxy
	# 	) + "plots/photometry_.png"
	f = find_galaxy(D.unbinned_flux, quiet=True, plot=plots)#, 
		#galaxy=galaxy.upper(), redshift=z, sav_fig=save_to)
	#f_err = find_galaxy(galaxy_data_error, quiet=True, plot=False)
	x_gals[i_gal] = f.xpeak # f.xmed?
	y_gals[i_gal] = f.ypeak # f.ymed?

	pa_gals[i_gal2] = 90 - f.theta
	ellipticity_gals[i_gal2] = f.eps

	# print "ellip: " + str(f.eps) #+ "+/-" + str(abs(f.eps-f_err.eps))
	# print "PA_photo: " + str(90-f.theta) #+ "+/-" + str(abs(f.theta-f_err.theta))


# ------------================ Lambda_R ==================----------

	beta = np.tan(D.yBar/D.xBar) # Angle between point and center of 
									# galaxy and RA axis.
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

	vel = D.components['stellar'].plot['vel']
	vel_lim = set_lims(vel, symmetric=True)
	mask = (vel < vel_lim[0]) + (vel > vel_lim[1])

	sigma = D.components['stellar'].plot['sigma']
	sigma_lim = set_lims(sigma, positive=True)
	mask += (sigma < sigma_lim[0]) + (sigma > sigma_lim[1])
	
	vel[mask], sigma[mask] = np.nan, np.nan

	# NB: numerator and denominator are in R_m order
	numerator = np.nancumsum(D.flux[R_m_sort] * 
		np.sqrt(D.xBar**2 + D.yBar**2)[R_m_sort] * np.abs(vel[R_m_sort]))

	denominator = np.nancumsum(D.flux[R_m_sort] * 
		np.sqrt(D.xBar**2 + D.yBar**2)[R_m_sort] * np.sqrt(vel**2 + sigma**2)[R_m_sort])

	lambda_R = numerator[R_m_sort.argsort()]/denominator[R_m_sort.argsort()]
	lambda_R[np.isnan(R_m)] = np.nan

	lambda_Re = interp1d(R_m[~np.isnan(R_m)], lambda_R[~np.isnan(R_m)],
		bounds_error=False, fill_value=(np.nan, 
			lambda_R[R_m_sort][~np.isnan(lambda_R[R_m_sort])][-1]))(R_e)

	print 'lambda_Re: ', lambda_Re
	lambda_Re_gals[i_gal2] = lambda_Re

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
# ------------========= Stellar Kinematics ===============----------
	save_to = "%s/plots/stellar_kinematics.png" % (output)
	k = fit_kinematic_pa(D.xBar - f.xpeak, D.yBar - f.ypeak, 
		np.array(D.components['stellar'].plot['vel']), quiet=True, plot=plots, 
		sav_fig=save_to)
	star_kine_pa_gals[i_gal2] = k[0]
# ------------=========== Gas Kinematics =================----------
# NB: this is not written for gas=2 or gas=3 options. 
	if D.gas == 1:
		save_to = "%s/plots/gas_kinematics.png" % (output)
		k = fit_kinematic_pa(D.xBar - f.xpeak, D.yBar - f.ypeak, 
			np.array(D.components['Hbeta'].plot['vel']), quiet=True, plot=plots, 
			sav_fig=save_to)
		gas_kine_pa_gals[i_gal2] = k[0]
# ------------============== Save results ================----------
	temp = "{0:12}{1:4}{2:4}"+''.join(['{%i:%i}'%(i+3,len(t)+1) for i, t in 
		enumerate(SN_gals.keys())])+'\n'

	SN_titles = list(SN_gals.keys())
	with open(galaxiesFile, 'w') as f:
		f.write(temp.format("Galaxy", "x", "y", *(s for s in SN_titles)))
		for i in range(len(galaxy_gals)):
			f.write(temp.format(galaxy_gals[i], str(int(x_gals[i])), 
				str(int(y_gals[i])), *(str(round(SN_gals[s][i],2)) for s in SN_titles)))



	temp = "{0:12}{1:10}{2:6}{3:9}{4:14}{5:14}\n"
	with open(galaxiesFile2, 'w') as f:
		f.write(temp.format("Galaxy", "lambda_Re", "eps", "pa", "star_kine_pa", 
			"gas_kine_pa"))
		for i in range(len(galaxy_gals2)):
			f.write(temp.format(galaxy_gals2[i], str(round(lambda_Re_gals[i],4)), 
				str(round(ellipticity_gals[i], 3)), str(round(pa_gals[i], 3)),
				str(round(star_kine_pa_gals[i],3)), str(round(gas_kine_pa_gals[i],3))))

	return D

##############################################################################

# Use of kinematics.py

if __name__ == '__main__':
	galaxy = 'ngc1399'
	discard = 0

	kinematics(galaxy, opt='kin', discard=discard, plots=False)