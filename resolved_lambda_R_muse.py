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
from Bin2 import Data
import disk_fit_functions_binned as dfn

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

	D = Data(galaxy, instrument='muse', opt='kin')
	opt = D.opt

	galaxy_gals, pa_gals = np.loadtxt('%s/Data/muse/analysis/galaxies2.txt' % 
		(cc.base_dir), usecols=(0,3), skiprows=1, unpack=True, dtype=str)
	i_gal = np.where(galaxy_gals == galaxy)[0][0]
	pa = float(pa_gals[i_gal])

	vel = D.components['stellar'].plot['vel']

	disk,pars=dfn.disk_fit_exp(D.xBar,D.yBar, vel, vel.uncert, verbose=False,
		grid_length=150)#,
	# 	, leeway=5., sigclip=None, grid_length=150, pa=pa)

	max_vel = np.nanmax(np.abs(disk))

	r = np.sqrt((D.xBar-center[0])**2 + (D.yBar-center[1])**2)
	m = ((vel > disk - max_vel) * (vel < disk + max_vel)) + r < max(r)/2
	vel = vel[m]
	r = r[m]

	s = np.argsort(r)
	lam_R = np.nancumsum((D.flux[m] * r * np.abs(vel))[s])/np.nancumsum(
	    (D.flux[m] * r * np.sqrt(vel**2 
	    + D.components['stellar'].plot['sigma'][m]**2))[s])

	R_e = get_R_e(galaxy)
	r *= abs(f[1].header['CD1_1'])*60**2/R_e
	r = r[s]
	file = '%s/Data/muse/analysis/%s/%s/lambda_R.txt' % (cc.base_dir, 
			galaxy, opt)

	with open(file, 'w') as f2:
		for i in range(len(r)):
			f2.write('%.4f   %.4f \n' %(r[i], lam_R[i]))
