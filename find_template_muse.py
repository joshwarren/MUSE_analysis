# ==================================================================
#  Find templates to use in fit
# ==================================================================
# warrenj 20150216 Process to analyse the reduced VIMOS data.
# warrenj 20160913 Ported to python

import numpy as np
import glob
from astropy.io import fits
from scipy import ndimage # for gaussian blur
from ppxf import ppxf
import ppxf_util as util
from errors2_muse import remove_anomalies, use_templates, determine_goodpixels, \
	get_stellar_templates, get_dataCubeDirectory
from checkcomp import checkcomp
cc= checkcomp()

quiet = True
c = 299792.458


def setup(galaxy, use_all_temp=False):
# ----------===============================================---------
# ----------============= Input parameters  ===============---------
# ----------===============================================---------
	dir = '%s/Data/muse'  %(cc.base_dir)
	templatesDirectory = '%s/models/miles_library' % (cc.home_dir)

	FWHM_gal = 2.3/(1+z) # MUSE documentation (R=2000 @ 4600A)
	moments = 4 # number of componants to calc with ppxf (see 
				# keyword moments in ppxf.pro for more details)
	degree = 4 # order of addative Legendre polynomial used to 
		   # correct the template continuum shape during the fit 

	# use the VIMOS fits
	data_file = "%s/Data/vimos/analysis/galaxies.txt" % (cc.base_dir)
	# different data types need to be read separetly
	z_gals, vel_gals, sig_gals = np.loadtxt(data_file, unpack=True, skiprows=1, 
		usecols=(1,2,3))
	galaxy_gals = np.loadtxt(data_file, skiprows=1, usecols=(0,),dtype=str)
	if galaxy != 'ngc1316':
		i_gal = np.where(galaxy_gals==galaxy)[0][0]
	else:
		# z should be good ngc1316 and ngc1399 are in Fornax Cluster
		i_gal = np.where(galaxy_gals=='ngc1399')[0][0]
	vel = vel_gals[i_gal]
	sig = sig_gals[i_gal]
	z = z_gals[i_gal]


## ----------=============== Miles library =================---------
	stellar_templates = get_stellar_templates(galaxy, FWHM_gal, use_all_temp=True)
	templates = stellar_templates.templates
	velscale = stellar_templates.velscale
## ----------========= Reading the spectrum  ===============---------

	dataCubeDirectory = get_dataCubeDirectory(galaxy)

	f = fits.open(dataCubeDirectory)
	galaxy_data, header = f[1].data, f[1].header
	galaxy_noise = f[2].data

	## write key parameters from header - can then be altered in future	
	CRVAL_spec = header['CRVAL3']
	CDELT_spec = header['CD3_3']
	s = galaxy_data.shape

	# Collapse to single spectrum
	gal_spec = np.nansum(galaxy_data, axis=(1,2))
	gal_noise = np.sqrt(np.nansum(galaxy_noise**2, axis=(1,2)))
## ----------========= Calibrating the spectrum  ===========---------
	lam = np.arange(s[0])*CDELT_spec + CRVAL_spec
	gal_spec, lam, cut = remove_anomalies(gal_spec, window=201, repeats=3, lam=lam, 
		set_range=set_range, return_cuts=True)
	lamRange = np.array([lam[0],lam[-1]])/(1+z)
	gal_noise = gal_noise[cut]


	## smooth spectrum to fit with templates resolution
	if FWHM_gal < stellar_templates.FWHM_tem:
		sigma = stellar_templates.FWHM_dif/2.355/CDELT_spec # Sigma difference in pixels
		gal_spec = ndimage.gaussian_filter1d(gal_spec, sigma)
		gal_noise = np.sqrt(ndimage.gaussian_filter1d(gal_noise**2, sigma))


	## rebin spectrum logarthmically
	bin_log, logLam_bin, _ = util.log_rebin(lamRange, gal_spec, velscale=velscale)
	gal_noise, logLam_bin, _ = util.log_rebin(lamRange, gal_noise**2, velscale=velscale)
	gal_noise = np.sqrt(gal_noise)

	noise = gal_noise + 0.0000000000001

	dv = (stellar_templates.logLam_template[0]-logLam_bin[0])*c # km/s
	# Find the pixels to ignore to avoid being distracted by gas emission
	#; lines or atmospheric absorbsion line.  
	goodpixels = determine_goodpixels(logLam_bin,lamRange_template,vel, z) 
	lambdaq = np.exp(logLam_bin)
	start = [vel, sig] # starting guess

	return templates, bin_log, noise, velscale, start, goodpixels, moments, \
		degree, dv, lambdaq, not quiet, quiet



## ----------===============================================---------
## ----------============== The bestfit part ===============---------
## ----------===============================================---------
def find_template(galaxy):
	print '     Finding templates to use'

	templates, bin_log, noise, velscale, start, goodpixels,	moments, degree, dv, \
		lambdaq, plot, quiet = setup(galaxy, use_all_temp=True)


	pp = ppxf(templates, bin_log, noise, velscale, start, 
			goodpixels=goodpixels, moments=moments, degree=degree, vsyst=dv, 
			lam=lambdaq, plot=not quiet, quiet=quiet)

	with open('%s/Data/muse/analysis/%s/templates.txt' % (cc.base_dir, galaxy), 'w') as f:
		for i in range(templates.shape[1]):
			if pp.weights[i] != 0.0:
				f.write(str(i) + '   ' + str(pp.weights[i]) + '\n')