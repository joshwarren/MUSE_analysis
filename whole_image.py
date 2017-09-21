# Routine to create a 'bin' which contains the entire image.
import numpy as np
import os
from astropy.io import fits
from errors2_muse import run_ppxf, set_params, get_dataCubeDirectory
from errors2 import apply_range
from checkcomp import checkcomp
cc = checkcomp()
import matplotlib.pyplot as plt

n_e = 100 # cm^-3
c = 299792.458 # speed of light in km/s
H0 = 70 # km/s/Mpc

def whole_image(galaxy):
	print galaxy
	if cc.device == 'glamdring':
		data_file = "%s/analysis/galaxies.txt" % (cc.base_dir)
	else:
		data_file = "%s/Data/vimos/analysis/galaxies.txt" % (cc.base_dir)
	galaxy_gals, z_gals = np.loadtxt(data_file, unpack=True, skiprows=1, 
		usecols=(0, 1), dtype=str)
	galaxy_gals = np.loadtxt(data_file, skiprows=1, usecols=(0,),dtype=str)
	i_gal = np.where(galaxy_gals==galaxy)[0][0]
	z = float(z_gals[i_gal])
	D = z*c/H0 # Mpc

	discard = 70
	f = fits.open(get_dataCubeDirectory(galaxy))
	spec = np.nansum(f[1].data[:,discard:-discard,discard:-discard], axis=(1,2))
	noise = np.sqrt(np.nansum(f[2].data[:,discard:-discard,discard:-discard], 
		axis=(1,2)))

	params = set_params(opt='pop', reps=0, temp_mismatch=True, produce_plot=True)

	lam = (np.arange(len(spec)) - (f[1].header['CRPIX3'] - 1)) * \
		f[1].header['CD3_3'] + f[1].header['CRVAL3']
	spec, lam, cut = apply_range(spec, lam=lam, return_cuts=True, 
		set_range=params.set_range)
	lamRange = np.array([lam[0],lam[-1]])
	noise = noise[cut]
	pp = run_ppxf(galaxy, spec, noise, lamRange, f[1].header['CD3_3'], params)



	pp.ax.ax2.plot(pp.lam, pp.matrix[:, pp.templatesToUse=='Hbeta'].flatten(), 'k')
	pp.fig.savefig('%s.png'%(galaxy))

	pp.noise = np.min([pp.noise, np.abs(pp.galaxy-pp.bestfit)],axis=0)

	OIII_spec = pp.matrix[:, pp.templatesToUse=='[OIII]5007d'].flatten() * \
		pp.weights[pp.templatesToUse=='[OIII]5007d']


	Hb_spec = pp.matrix[:, pp.templatesToUse=='Hbeta'].flatten() * \
		pp.weights[pp.templatesToUse=='Hbeta']
	Hb_flux = np.trapz(Hb_spec, x=pp.lam)
	Ha_flux = 2.86 * Hb_flux

	print 'From Hbeta'
	Mass = 280 * (D/10)**2 * (Ha_flux*10**-15/10**-14) * (1000/n_e) # Solar masses
	if max(OIII_spec)/np.median(pp.noise[
		(pp.lam < 5007./(1 + (pp.sol[1][0] - 300)/c)) *
		(pp.lam > 5007./(1 + (pp.sol[1][0] + 300)/c))]) > 4:

		if max(Hb_spec)/np.median(pp.noise[
			(pp.lam < 4861./(1 + (pp.sol[1][0] - 300)/c)) *
			(pp.lam > 4861./(1 + (pp.sol[1][0] + 300)/c))]) > 2.5:

			print '    %.2f log10(Solar Masses)' % (np.log10(Mass))
		else:
			print '    <%.2f log10(Solar Masses)' % (np.log10(Mass))
	else:
		print '    <%.2f log10(Solar Masses)' % (np.log10(Mass))

	print 'Direct from Halpha'
	Ha_spec = pp.matrix[:, pp.templatesToUse=='Halpha'].flatten() * \
		pp.weights[pp.templatesToUse=='Halpha']
	Ha_flux = np.trapz(Ha_spec, x=pp.lam)
	Mass = 280 * (D/10)**2 * (Ha_flux*10**-15/10**-14) * (1000/n_e) # Solar masses

	if max(OIII_spec)/np.median(pp.noise[
		(pp.lam < 5007./(1 + (pp.sol[1][0] - 300)/c)) *
		(pp.lam > 5007./(1 + (pp.sol[1][0] + 300)/c))]) > 4:

		if max(Ha_spec)/np.median(pp.noise[
			(pp.lam < 6563./(1 + (pp.sol[1][0] - 300)/c)) *
			(pp.lam > 6563./(1 + (pp.sol[1][0] + 300)/c))]) > 2.5:
			print '    %.2f log10(Solar Masses)' % (np.log10(Mass))
		else:
			print '    <%.2f log10(Solar Masses)' % (np.log10(Mass))
	else:
		print '    <%.2f log10(Solar Masses)' % (np.log10(Mass))

	print 'Balmer decrement: %.2f' % (Ha_flux/Hb_flux)


	print 'Hbeta %.2f Halpha %.2f' % (Hb_flux, Ha_flux)
	print 'max Hb %.2f Ha %.2f' %(Hb_spec.max(), Ha_spec.max())





if __name__=='__main__':
	galaxies = ['ic1459', 'ic4296', 'ngc1316', 'ngc1399']
	for g in galaxies:
		whole_image(g)