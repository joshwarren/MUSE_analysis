## ==================================================================
## 		Stellar population within inner 1 arcsec
## ==================================================================
# The entire pipeline - making use of error2_muse and pop_muse routines.


from checkcomp import checkcomp
cc = checkcomp()
import numpy as np 
from astropy.io import fits
from errors2_muse import get_dataCubeDirectory, run_ppxf, set_params, remove_anomalies
from pop_muse import population

res = 0.2 # MUSE spatial resolution

def get_specFromAperture(galaxy, app_size=1.0, inside=True):
	f = fits.open(get_dataCubeDirectory(galaxy))
	s = f[1].data.shape

	x = np.arange(s[1]).repeat(s[2]).reshape(s[1],s[2])
	y = np.tile(np.arange(s[2]),s[1]).reshape(s[1],s[2])

	galaxy_gals = np.loadtxt('%s/Data/muse/analysis/galaxies.txt' % (cc.base_dir),
		unpack=True, skiprows=1, usecols=(0,), dtype=str)
	x_cent, y_cent = np.loadtxt('%s/Data/muse/analysis/galaxies.txt' % (cc.base_dir),
		unpack=True, skiprows=1, usecols=(1,2), dtype=int)
	i_gal = np.where(galaxy_gals==galaxy)[0][0]
	x_cent = x_cent[i_gal]
	y_cent = y_cent[i_gal]

	# Distance in arcsec from central spaxel
	d = np.sqrt((x - x_cent)**2 + (y - y_cent)**2) * res
	if inside:
		mask = d <= app_size
	else:
		mask = d > app_size

	return np.nansum(f[1].data[:, mask], axis=1), \
		np.sqrt(np.nansum(f[2].data[:, mask]**2, axis=1)), \
		np.arange(s[0])*f[1].header['CD3_3'] + f[1].header['CRVAL3']


def KDC_pop(galaxy):
	params = set_params(opt='pop')

	data_file = "%s/Data/vimos/analysis/galaxies.txt" % (cc.base_dir)
	z_gals, vel_gals, sig_gals = np.loadtxt(data_file, unpack=True, skiprows=1, 
		usecols=(1,2,3))
	galaxy_gals = np.loadtxt(data_file, skiprows=1, usecols=(0,),dtype=str)
	i_gal = np.where(galaxy_gals==galaxy)[0][0]
	vel = vel_gals[i_gal]
	sig = sig_gals[i_gal]
	z = z_gals[i_gal]


	spec, noise, lam = get_specFromAperture(galaxy, app_size=1.0)
	CD = lam[1] - lam[0]
	spec, lam, cut = remove_anomalies(spec, window=201, repeats=3, 
		lam=lam, return_cuts=True, set_range=params.set_range, n_sigma=2)
	noise = noise[cut]
	lamRange = np.array([lam[0],lam[-1]])/(1+z)

	pp = run_ppxf(galaxy, spec, noise, lamRange, CD, vel, sig, z, 1, 'kin', params, 
		produce_plot=False)

	pop = population(pp=pp, galaxy=galaxy)


	# Outside apperture
	print 'outside'
	spec, noise, lam = get_specFromAperture(galaxy, app_size=1.0, inside=False)
	CD = lam[1] - lam[0]
	spec, lam, cut = remove_anomalies(spec, window=201, repeats=3, 
		lam=lam, return_cuts=True, set_range=params.set_range, n_sigma=2)
	noise = noise[cut]
	lamRange = np.array([lam[0],lam[-1]])/(1+z)

	f = pop.fig
	ax = pop.ax

	del pop

	pp_outside = run_ppxf(galaxy, spec, noise, lamRange, CD, vel, sig, z, 1, 'kin', params,
		produce_plot=False)
	print 'pp done'
	pop_outside = population(pp=pp_outside, galaxy=galaxy)
	print 'pop_outside done'
	pop_outside.plot_probability_distribution(save=True, f=f, ax_array=ax)

	del pop_outside




if __name__ == '__main__':
	for gal in ['ic1459', 'ic4296', 'ngc1316', 'ngc1399']:
		KDC_pop(gal)


