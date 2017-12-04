# Routine to create a 'bin' which contains the entire image.
from checkcomp import checkcomp
cc = checkcomp()
if cc.remote:
	import matplotlib # 20160202 JP to stop lack-of X-windows error
	matplotlib.use('Agg') # 20160202 JP to stop lack-of X-windows error
import numpy as np
import os
from astropy.io import fits
from errors2_muse import run_ppxf, set_params, get_dataCubeDirectory
from errors2 import apply_range
import matplotlib.pyplot as plt
from whole_image import get_Mass
from Bin import trapz_uncert
from global_mg_sigma import in_aperture


n_e = 100 # cm^-3
c = 299792.458 # speed of light in km/s
H0 = 70 # km/s/Mpc

def whole_image(galaxy, verbose=False):
	print galaxy
	max_reps = 100


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

	data_file = "%s/Data/muse/analysis/galaxies.txt" % (cc.base_dir)
	x_gals, y_gals = np.loadtxt(data_file, unpack=True, skiprows=1, 
		usecols=(1, 2), dtype=int)
	galaxy_gals = np.loadtxt(data_file, unpack=True, skiprows=1, 
		usecols=(0,), dtype=str)
	i_gal = np.where(galaxy_gals==galaxy)[0][0]
	centre = (x_gals[i_gal], y_gals[i_gal])

	limits_file = '%s/Data/muse/analysis/galaxies_gasMass.txt' %(cc.base_dir)
	galaxy_gals, mass, e_mass, bulmer, e_bulmer = np.loadtxt(limits_file, 
		unpack=True, dtype=str, skiprows=1)
	i_gal = np.where(galaxy_gals==galaxy)[0][0]

	max_radius = 90
	radius = float(max_radius)
	# f = fits.open(get_dataCubeDirectory(galaxy))
	f = fits.open(get_dataCubeDirectory(galaxy)[:-5]+'2.fits')
	while radius > 2:
		mask = in_aperture(centre[0], centre[1], radius, instrument='muse')

		spec = f[1].data
		noise = f[2].data

		spec[np.isnan(spec)] = 0
		noise[np.isnan(noise)] = 0

		spec = np.einsum('ijk,jk->i', spec, mask)/np.sum(mask)
		noise = np.sqrt(np.einsum('ijk,jk->i', noise**2, mask))/np.sum(mask)

		if radius == max_radius:
			reps = max_reps
			params = set_params(opt='pop', reps=reps, temp_mismatch=True, 
				produce_plot=False)
		

		lam = (np.arange(len(spec)) - (f[1].header['CRPIX3'] - 1)) * \
			f[1].header['CD3_3'] + f[1].header['CRVAL3']
		spec, lam, cut = apply_range(spec, lam=lam, return_cuts=True, 
			set_range=params.set_range)
		lamRange = np.array([lam[0],lam[-1]])
		noise = noise[cut]
		pp = run_ppxf(galaxy, spec, noise, lamRange, f[1].header['CD3_3'], 
			params)

		# pp.ax.ax2.plot(pp.lam, pp.matrix[:, 
		# 	pp.templatesToUse=='Hbeta'].flatten(), 'k')
		# pp.fig.savefig('%s.png'%(galaxy))

		# pp.noise = np.min([pp.noise, np.abs(pp.galaxy-pp.bestfit)],axis=0)

		OIII_spec = pp.matrix[:, pp.templatesToUse=='[OIII]5007d'].flatten(
			) * pp.weights[pp.templatesToUse=='[OIII]5007d']


		Hb_spec = pp.matrix[:, pp.templatesToUse=='Hbeta'].flatten() * \
			pp.weights[pp.templatesToUse=='Hbeta']
		Hb_flux = np.trapz(Hb_spec, x=pp.lam)
		# Ha_flux = 2.86 * Hb_flux

		# print 'From Hbeta'
		# Mass = get_mass(Ha_flux, D, instrument='muse') # Solar masses
		# if max(OIII_spec)/np.median(pp.noise[
		# 	(pp.lam < 5007./(1 + (pp.sol[1][0] - 300)/c)) *
		# 	(pp.lam > 5007./(1 + (pp.sol[1][0] + 300)/c))]) > 4:

		# 	if max(Hb_spec)/np.median(pp.noise[
		# 		(pp.lam < 4861./(1 + (pp.sol[1][0] - 300)/c)) *
		# 		(pp.lam > 4861./(1 + (pp.sol[1][0] + 300)/c))]) > 2.5:

		# 		print '    %.2f log10(Solar Masses)' % (np.log10(Mass))
		# 	else:
		# 		print '    <%.2f log10(Solar Masses)' % (np.log10(Mass))
		# else:
		# 	print '    <%.2f log10(Solar Masses)' % (np.log10(Mass))

		Ha_spec = pp.matrix[:, pp.templatesToUse=='Halpha'].flatten() * \
			pp.weights[pp.templatesToUse=='Halpha']
		Ha_flux = np.trapz(Ha_spec, x=pp.lam)
		if reps == max_reps:
			Hb_spec_uncert = pp.MCgas_uncert_spec[
				pp.templatesToUse[pp.component!=0]=='Hbeta', :].flatten()
			Hb_flux_uncert = trapz_uncert(Hb_spec_uncert, x=pp.lam)

			Ha_spec_uncert = pp.MCgas_uncert_spec[
				pp.templatesToUse[pp.component!=0]=='Halpha', :].flatten()
			Ha_flux_uncert = trapz_uncert(Ha_spec_uncert, x=pp.lam)

		if max(OIII_spec)/np.median(pp.noise[
			(pp.lam < 5007./(1 + (pp.sol[1][0] - 300)/c)) *
			(pp.lam > 5007./(1 + (pp.sol[1][0] + 300)/c))]) > 4:

			if max(Ha_spec)/np.median(pp.noise[
				(pp.lam < 6563./(1 + (pp.sol[1][0] - 300)/c)) *
				(pp.lam > 6563./(1 + (pp.sol[1][0] + 300)/c))]) > 2.5:

				if reps==max_reps:
					Mass = get_Mass(Ha_flux, D, instrument='muse')
					e_Mass = get_Mass(Ha_flux_uncert, D, instrument='muse')
					
					mass[i_gal] = str(round(np.log10(Mass),4))
					e_mass[i_gal] =  str(round(np.abs(e_Mass/Mass/
						np.log(10)), 4))
					if verbose:
						print '%s +/- %s log10(Solar Masses)' % (
							mass[i_gal], e_mass[i_gal])
						# fig, ax = plt.subplots(2)
						# pp.ax = ax[0]
						# from ppxf import create_plot
						# fig, ax = create_plot(pp).produce 
						# ax.set_xlim([4800, 4900])
						# ax.legend()

						# pp.ax = ax[1]
						# from ppxf import create_plot
						# fig, ax = create_plot(pp).produce 
						# ax.set_xlim([6500, 6600])

						# fig.savefig('%s.png'%(galaxy))
					radius = -1
				else: # Repeat but calculate uncert
					reps = max_reps

				if max(Hb_spec)/np.median(pp.noise[
					(pp.lam < 4861./(1 + (pp.sol[1][0] - 300)/c)) *
					(pp.lam > 4861./(1 + (pp.sol[1][0] + 300)/c))]) > 2.5:
					b = Ha_flux/Hb_flux
					e_bulmer[i_gal] = str(round(b * np.sqrt(
						(Ha_flux_uncert/Ha_flux)**2 + 
						(Hb_flux_uncert/Hb_flux)**2), 2))
					bulmer[i_gal] = str(round(b, 2))
				else:
					b = Ha_flux/Hb_flux
					e_bulmer[i_gal] = str(round(b * np.sqrt(
						(Ha_flux_uncert/Ha_flux)**2 + 
						(Hb_flux_uncert/Hb_flux)**2), 2))
					bulmer[i_gal] = '<'+str(round(b, 2))
			else:
				if radius == max_radius:
					Mass = get_Mass(Ha_flux, D, instrument='muse')
					e_Mass = get_Mass(Ha_flux_uncert, D, instrument='muse')

					mass[i_gal] = '<'+str(round(np.log10(Mass),4))
					e_mass[i_gal] =  str(round(np.abs(e_Mass/Mass/
						np.log(10)), 4))
					b = Ha_flux/Hb_flux
					e_bulmer[i_gal] = str(round(b * np.sqrt(
						(Ha_flux_uncert/Ha_flux)**2 + 
						(Hb_flux_uncert/Hb_flux)**2), 2))
					bulmer[i_gal] = '<'+str(round(b, 2))
					if verbose:
						print '<%s +/- %s log10(Solar Masses)' % (
							mass[i_gal], e_mass[i_gal])

				radius -= 5
				reps = 0


		else:
			# if radius == max_radius:
			Mass = get_Mass(Ha_flux, D, instrument='muse')
			e_Mass = get_Mass(Ha_flux_uncert, D, instrument='muse')

			mass[i_gal] = '>'+str(round(np.log10(Mass),4))
			e_mass[i_gal] =  str(round(np.abs(e_Mass/Mass/
				np.log(10)), 4))
			b = Ha_flux/Hb_flux
			e_bulmer[i_gal] = str(round(b * np.sqrt(
				(Ha_flux_uncert/Ha_flux)**2 + 
				(Hb_flux_uncert/Hb_flux)**2), 2))
			bulmer[i_gal] = '<'+str(round(b, 2))
			if verbose:
				print '>%s +/- %s log10(Solar Masses)' % (
				mass[i_gal], e_mass[i_gal])
		radius -= 5
		reps = 0

		params = set_params(opt='pop', reps=reps, temp_mismatch=True, 
			produce_plot=False)

	temp = "{0:12}{1:10}{2:10}{3:10}{4:10}\n"
	with open(limits_file, 'w') as l:
		l.write(temp.format('Galaxy', 'Mass', 'e_Mass', 'Bul_dec', 
			'e_Bul_dec'))
		for i in range(len(galaxy_gals)):
			l.write(temp.format(galaxy_gals[i], mass[i], e_mass[i], bulmer[i],
				e_bulmer[i]))

if __name__=='__main__':
	# galaxies = ['ic1459', 'ic4296', 'ngc1316', 'ngc1399']
	galaxies = ['ngc1316', 'ngc1399']
	for g in galaxies:
		whole_image(g, verbose=True)