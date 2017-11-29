import numpy as np 
from astropy.io import fits
from checkcomp import checkcomp
cc = checkcomp()
if cc.remote:
	import matplotlib # 20160202 JP to stop lack-of X-windows error
	matplotlib.use('Agg') # 20160202 JP to stop lack-of X-windows error
	import matplotlib.pyplot as plt # used for plotting
else:
	import matplotlib.pyplot as plt # used for plotting
from errors2_muse import get_dataCubeDirectory
from classify import get_R_e
from prefig import Prefig
Prefig()

import cPickle as pickle

opt = 'pop'
for galaxy in [
			'ic1459',
			'ic4296',
			'ngc1316',
			'ngc1399']:
	print galaxy 

	f = fits.open(get_dataCubeDirectory(galaxy))
	data_file = '%s/Data/muse/analysis/galaxies.txt' % (cc.base_dir)
	x_cent_gals, y_cent_gals = np.loadtxt(data_file, unpack=True, 
	    skiprows=1, usecols=(1,2), dtype=int)
	galaxy_gals = np.loadtxt(data_file, skiprows=1, usecols=(0,),dtype=str)
	i_gal = np.where(galaxy_gals==galaxy)[0][0]
	center = (x_cent_gals[i_gal], y_cent_gals[i_gal])

	pickleFile = open('/Data/muse/analysis/%s/%s/pickled/dataObj.pkl' % (
		galaxy,opt), 'rb')
	D = pickle.load(pickleFile)
	pickleFile.close()

	r = np.sqrt((D.xBar-center[0])**2 + (D.yBar-center[1])**2)
	s = np.argsort(r)
	lam_R = np.nancumsum((D.flux * r * 
		np.abs(D.components['stellar'].plot['vel']))[s])/np.nancumsum(
	    (D.flux * r * np.sqrt(D.components['stellar'].plot['vel']**2 + 
	    D.components['stellar'].plot['sigma']**2))[s])

	R_e = get_R_e(galaxy)
	r *= abs(f[1].header['CD1_1'])*60**2/R_e
	r = r[s]
	file = '%s/Data/muse/analysis/%s/%s/lambda_R.txt' % (cc.base_dir, 
			galaxy, opt)

	with open(file, 'w') as f2:
		for i in range(len(r)):
			f2.write('%.4f   %.4f \n' %(r[i], lam_R[i]))
