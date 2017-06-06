


from checkcomp import checkcomp
cc = checkcomp()
import numpy as np 
from astropy.io import fits
from errors2_muse import get_dataCubeDirectory, run_ppxf, set_params
res = 0.2 # MUSE spatial resolution

def get_specWithin(galaxy, app_size=1.0):
	f = fits.open(get_dataCubeDirectory(galaxy))
	s = f[1].data.shape

	x = np.arange(s[1]).repeat(s[2]).reshape(s[1],s[2])
	y = np.tile(np.arange(s[2]),s[1]).reshape(s[1],s[2])

	galaxy_gals = np.loadtxt('%s/Data/muse/analysis/galaxies.txt' % (cc.base_dir),
		unpack=True, skiprow=1, usecols=(0,), dtype=str)
	x_cent, y_cent = np.loadtxt('%s/Data/muse/analysis/galaxies.txt' % (cc.base_dir),
		unpack=True, skiprow=1, usecols=(1,2), dtype=int)
	i_gal = np.where(galaxy_gals==galaxy)[0][0]
	x_cent = x_cent[i_gal]
	y_cent = y_cent[i_gal]

	# Distance in arcsec from central spaxel
	d = np.sqrt((x - x_cent)**2 + (y - y_cent)**2) * res

	return np.nansum(f[1].data[:, d < app_size], axis=1), \
		np.sqrt(np.nansum(f[2].data[:, d < app_size]**2, axis=1))



def KDC_pop(galaxy):
	spec, noise = get_specWithin(galaxy, app_size=1.0)
	params = set_params()

	data_file = "%s/Data/vimos/analysis/galaxies.txt" % (cc.base_dir)
	z_gals, vel_gals, sig_gals = np.loadtxt(data_file, unpack=True, skiprows=1, 
		usecols=(1,2,3))
	galaxy_gals = np.loadtxt(data_file, skiprows=1, usecols=(0,),dtype=str)
	i_gal = np.where(galaxy_gals==galaxy)[0][0]
	vel = vel_gals[i_gal]
	sig = sig_gals[i_gal]
	z = z_gals[i_gal]

	pp = run_ppxf(spec, noise, lamRange, vel, sig, z, 1, 'kin', params)


	






if __name__ == '__main__':
	KDC_pop('ic4296')


