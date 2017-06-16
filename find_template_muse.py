# ==================================================================
#  Find templates to use in fit
# ==================================================================
# warrenj 20150216 Process to analyse the reduced VIMOS data.
# warrenj 20160913 Ported to python

import numpy as np
from astropy.io import fits
from checkcomp import checkcomp
cc= checkcomp()
from errors2_muse import run_ppxf, set_params, remove_anomalies, get_dataCubeDirectory
c = 299792.458

def find_template(galaxy, set_range=None):
	params = set_params()
	params.gas = 0
	params.set_range = set_range
## ----------========= Reading the spectrum  ===============---------

	dataCubeDirectory = get_dataCubeDirectory(galaxy)

	f = fits.open(dataCubeDirectory)

	## write key parameters from header - can then be altered in future	
	CRVAL_spec = f[1].header['CRVAL3']
	CDELT_spec = f[1].header['CD3_3']
	s = f[1].data.shape

	# Collapse to single spectrum
	gal_spec = np.zeros(s[0])
	gal_noise = np.zeros(s[0])

	for i in xrange(s[0]):
		gal_spec[i] = np.nansum(f[1].data[i, int(s[1]/2.0-50):int(s[1]/2.0+50), 
			int(s[2]/2.0-50):int(s[2]/2.0+50)])
		gal_noise[i] = np.sqrt(np.nansum(f[2].data[i, int(s[1]/2.0-50):int(s[1]/2.0+50), 
			int(s[2]/2.0-50):int(s[2]/2.0+50)])**2)

	del f

## ----------========= Calibrating the spectrum  ===========---------
	lam = np.arange(s[0])*CDELT_spec + CRVAL_spec
	gal_spec, lam, cut = remove_anomalies(gal_spec, window=201, repeats=0, 
		lam=lam, return_cuts=True, set_range=params.set_range, n_sigma=2)
	lamRange = np.array([lam[0],lam[-1]])
	gal_noise = gal_noise[cut]


	pp = run_ppxf(galaxy, gal_spec, gal_noise, lamRange, CDELT_spec, params, 
		use_all_temp=True)
	pp.fig.savefig('%s/Data/muse/analysis/%s/find_temp.png'% (cc.base_dir, galaxy))

	with open('%s/Data/muse/analysis/%s/templates.txt' % (cc.base_dir, galaxy), 'w') as f:
		for i in range(len(pp.component)):
			if pp.weights[i] != 0.0:
				f.write(str(i) + '   ' + str(pp.weights[i]) + '\n')