## ==================================================================
## 		Plot stellar population maps
## ==================================================================
## warrenj 20170331 Routine to plot the stellar populations found by pop.py
## on Glamdring. 

from checkcomp import checkcomp
cc = checkcomp()
if cc.remote:
	import matplotlib # 20160202 JP to stop lack-of X-windows error
	matplotlib.use('Agg') # 20160202 JP to stop lack-of X-windows error
	import matplotlib.pyplot as plt # used for plotting
else:
	import matplotlib.pyplot as plt # used for plotting
import cPickle as pickle
from plot_velfield_nointerp import plot_velfield_nointerp 
import numpy as np 
import os
from plot_results_muse import set_lims, add_
from errors2_muse import get_dataCubeDirectory, set_params, run_ppxf
from pop import get_absorption, population
from classify import get_R_e
from errors2 import apply_range
from astropy.io import fits
from prefig import Prefig
Prefig(size=(16*2,12*3), transparent=False)

def plot_stellar_pop(galaxy, method='median', opt='pop', D=None, overplot={},
	gradient=True):
	print 'Plotting stellar population'

	if cc.device == 'glamdring': 
		vin_dir = '%s/analysis_muse/%s/%s/pop' % (cc.base_dir, galaxy, opt)
		data_file = '%s/analysis_muse/galaxies.txt' % (cc.base_dir)
	else: 
		vin_dir = '%s/Data/muse/analysis/%s/%s/pop' % (cc.base_dir, galaxy, opt)
		data_file =  "%s/Data/muse/analysis/galaxies.txt" % (cc.base_dir)

	file_headings = np.loadtxt(data_file, dtype=str)[0]
	col = np.where(file_headings=='SN_%s' % (opt))[0][0]
	x_cent_gals, y_cent_gals, SN_target_gals = np.loadtxt(data_file, unpack=True, 
		skiprows=1, usecols=(1,2, col), dtype='int,int,float')
	galaxy_gals = np.loadtxt(data_file, skiprows=1, usecols=(0,),dtype=str)
	i_gal = np.where(galaxy_gals==galaxy)[0][0]
	SN_target=SN_target_gals[i_gal]
	center = (x_cent_gals[i_gal], y_cent_gals[i_gal])

	data_file =  "%s/Data/vimos/analysis/galaxies.txt" % (cc.base_dir)
	z_gals = np.loadtxt(data_file, unpack=True, skiprows=1, 
		usecols=(1))
	galaxy_gals = np.loadtxt(data_file, skiprows=1, usecols=(0,),dtype=str)
	i_gal = np.where(galaxy_gals==galaxy)[0][0]
	z = z_gals[i_gal]

	# Load pickle file from pickler.py
	out_dir = '%s/Data/muse/analysis' % (cc.base_dir)
	output = "%s/%s/%s" % (out_dir, galaxy, opt)
	out_plots = "%s/plots/population" % (output)
	if not os.path.exists(out_plots): os.makedirs(out_plots)
	 
	if D is None and gradient != 'only':
		pickle_file = '%s/pickled' % (output)
		pickleFile = open("%s/dataObj.pkl" % (pickle_file), 'rb')
		D = pickle.load(pickleFile)
		pickleFile.close()

	f = fits.open(get_dataCubeDirectory(galaxy))
	header = f[1].header
	if not gradient: f.close()

	if gradient != 'only':
		age = np.zeros(D.number_of_bins)
		met = np.zeros(D.number_of_bins)
		alp = np.zeros(D.number_of_bins)
		unc_age = np.zeros(D.number_of_bins)
		unc_met = np.zeros(D.number_of_bins)
		unc_alp = np.zeros(D.number_of_bins)

		if method == 'median':
			for i in xrange(D.number_of_bins):
				ag, me, al = np.loadtxt('%s/%i.dat' % (vin_dir, i), 
					unpack=True)
				
				age[i] = ag[0]
				unc_age[i] = ag[1]
				met[i] = me[0]
				unc_met[i] = me[1]
				alp[i] = al[0]
				unc_alp[i] = al[1]

			title = '%s median' %(galaxy.upper())
			u_title = '%s standard deviation' % (galaxy.upper())


			

		elif method == 'mostlikely':
			# from peakdetect import peakdetect

			age1 = np.zeros(D.number_of_bins)
			met1 = np.zeros(D.number_of_bins)
			alp1 = np.zeros(D.number_of_bins)

			age2 = np.zeros(D.number_of_bins)
			met2 = np.zeros(D.number_of_bins)
			alp2 = np.zeros(D.number_of_bins)
			for i in xrange(D.number_of_bins):
				ag, me, al = np.loadtxt('%s/distribution/%i.dat' % (vin_dir, 
					i), unpack=True)

				for plot, unc_plot, pop in zip([age,met,alp],
					[unc_age,unc_met,unc_alp], [ag,me,al]):

					hist = np.histogram(pop, bins=40)
					x = (hist[1][0:-1]+hist[1][1:])/2
					hist = hist[0]
					# peaks = np.array(peakdetect(hist, x_axis=x, lookahead=4)[0])
					# plot[i] = peaks[np.argmax(peaks[:,1]), 0]
					plot[i] = x[np.argmax(hist)]

					gt_fwhm = hist >= np.max(hist)/2
					unc_plot[i] = np.max(x[gt_fwhm]) - np.min(x[gt_fwhm])

			title = '%s mostlikely' % (galaxy.upper())
			u_title = '%s FWHM' % (galaxy.upper())

	if gradient:
		figs = {}
		axs = {}
		rad = {}
		rad_err = {}
		for i in ['age', 'met', 'alp']:
			fig, ax = plt.subplots()
			figs[i] = fig
			axs[i] = ax
			rad[i] = []
			rad_err[i] = []

	if gradient != 'only':
		# Age
		ax = plot_velfield_nointerp(D.x, D.y, D.bin_num, 
			D.xBar, D.yBar, age, header, nodots=True, colorbar=True, label='Age (Gyrs)', 
			vmin=0, vmax=15, title=title + ' Age', cmap='gnuplot2', 
			flux_unbinned=D.unbinned_flux, signal_noise=D.SNRatio, 
			signal_noise_target=SN_target, center=center, redshift=z)
		if overplot:
			for o, color in overplot.iteritems():
				add_(o, color, ax, galaxy)
		plt.gcf().savefig('%s/Age.png' % (out_plots))
		plt.close()

		plot_velfield_nointerp(D.x, D.y, D.bin_num, 
			D.xBar, D.yBar, unc_age, header, nodots=True, colorbar=True, label='Age (Gyrs)', 
			vmin=0, vmax=15, title=u_title+' Age', cmap='gnuplot2', 
			flux_unbinned=D.unbinned_flux, signal_noise=D.SNRatio, close=True,
			signal_noise_target=SN_target, center=center, save='%s/Age_uncert.png'%(out_plots))

		# Metalicity
		ax = plot_velfield_nointerp(D.x, D.y, D.bin_num, 
			D.xBar, D.yBar, met, header, nodots=True, colorbar=True,
			label='Metalicity [Z/H]', vmin=-2.25, vmax=0.67, title=title+' Metalicity', 
			cmap='gnuplot2', flux_unbinned=D.unbinned_flux, 
			signal_noise=D.SNRatio, signal_noise_target=SN_target, center=center)
		if overplot:
			for o, color in overplot.iteritems():
				add_(o, color, ax, galaxy)
		plt.gcf().savefig('%s/Metalicity.png' % (out_plots))
		plt.close()

		plot_velfield_nointerp(D.x, D.y, D.bin_num, D.xBar, D.yBar, unc_met, header, 
			nodots=True, colorbar=True, label='Metalicity', vmin=0, vmax=0.67+2.25, 
			title=u_title+' Metalicity [Z/H]', cmap='gnuplot2', 
			flux_unbinned=D.unbinned_flux, signal_noise=D.SNRatio, 
			signal_noise_target=SN_target, center=center, 
			save='%s/Metalicity_uncert.png'%(out_plots), close=True)

		# Alpha
		ax = plot_velfield_nointerp(D.x, D.y, D.bin_num, D.xBar, D.yBar, alp, header, 
			nodots=True, colorbar=True, label='Element Ratio [alpha/Fe]', vmin=-0.3, vmax=0.5, 
			title=title+' Alpha Enhancement', cmap='gnuplot2', flux_unbinned=D.unbinned_flux, 
			signal_noise=D.SNRatio, signal_noise_target=SN_target, center=center)
		if overplot:
			for o, color in overplot.iteritems():
				add_(o, color, ax, galaxy)
		plt.gcf().savefig('%s/Alpha.png' % (out_plots))
		plt.close()


		plot_velfield_nointerp(D.x, D.y, D.bin_num, D.xBar, D.yBar, unc_alp, header, 
			nodots=True, colorbar=True, label='Element Ratio [alpha/Fe]', vmin=0, 
			vmax=0.5+0.3, title=u_title+' Alpha Enhancement', cmap='gnuplot2', 
			flux_unbinned=D.unbinned_flux, signal_noise=D.SNRatio, 
			signal_noise_target=SN_target, center=center,
			save='%s/Alpha_uncert.png'%(out_plots), close=True)


		# Detailed (no clip on color axis)
		out_plots = "%s/plots/population_detail" % (output)
		if not os.path.exists(out_plots): os.makedirs(out_plots)
		# Age
		vmin, vmax = set_lims(age)
		ax = plot_velfield_nointerp(D.x, D.y, D.bin_num, 
			D.xBar, D.yBar, age, header, nodots=True, colorbar=True, label='Age (Gyrs)', 
			title=title + ' Age', cmap='gnuplot2', vmin=vmin, vmax=vmax, 
			flux_unbinned=D.unbinned_flux, signal_noise=D.SNRatio, 
			signal_noise_target=SN_target, center=center, redshift=z)
		if overplot:
			for o, color in overplot.iteritems():
				add_(o, color, ax, galaxy)
		plt.gcf().savefig('%s/Age.png' % (out_plots))
		plt.close()

		vmin, vmax = set_lims(unc_age)
		plot_velfield_nointerp(D.x, D.y, D.bin_num, 
			D.xBar, D.yBar, unc_age, header, nodots=True, colorbar=True, label='Age (Gyrs)', 
			title=u_title+' Age', cmap='gnuplot2', vmin=vmin, vmax=vmax,
			flux_unbinned=D.unbinned_flux, signal_noise=D.SNRatio, close=True,
			signal_noise_target=SN_target, center=center, save='%s/Age_uncert.png'%(out_plots))

		# Metalicity
		vmin, vmax = set_lims(met)
		ax = plot_velfield_nointerp(D.x, D.y, D.bin_num, 
			D.xBar, D.yBar, met, header, nodots=True, colorbar=True,
			label='Metalicity [Z/H]', title=title+' Metalicity', vmin=vmin, vmax=vmax,
			cmap='gnuplot2', flux_unbinned=D.unbinned_flux, 
			signal_noise=D.SNRatio, signal_noise_target=SN_target, center=center)
		if overplot:
			for o, color in overplot.iteritems():
				add_(o, color, ax, galaxy)
		plt.gcf().savefig('%s/Metalicity.png' % (out_plots))
		plt.close()

		vmin, vmax = set_lims(unc_met)
		plot_velfield_nointerp(D.x, D.y, D.bin_num, D.xBar, D.yBar, unc_met, header, 
			nodots=True, colorbar=True, label='Metalicity', vmin=vmin, vmax=vmax, 
			title=u_title+' Metalicity [Z/H]', cmap='gnuplot2', 
			flux_unbinned=D.unbinned_flux, signal_noise=D.SNRatio, 
			signal_noise_target=SN_target, center=center, 
			save='%s/Metalicity_uncert.png'%(out_plots),close=True)

		# Alpha
		vmin, vmax = set_lims(alp)
		ax = plot_velfield_nointerp(D.x, D.y, D.bin_num, D.xBar, D.yBar, alp, header, 
			nodots=True, colorbar=True, label='Element Ratio [alpha/Fe]', vmin=vmin, vmax=vmax,
			title=title+' Alpha Enhancement', cmap='gnuplot2', flux_unbinned=D.unbinned_flux, 
			signal_noise=D.SNRatio, signal_noise_target=SN_target, center=center)
		if overplot:
			for o, color in overplot.iteritems():
				add_(o, color, ax, galaxy)
		plt.gcf().savefig('%s/Alpha.png' % (out_plots))
		plt.close()

		vmin, vmax = set_lims(unc_alp)
		plot_velfield_nointerp(D.x, D.y, D.bin_num, D.xBar, D.yBar, unc_alp, header, 
			nodots=True, colorbar=True, label='Element Ratio [alpha/Fe]', 
			title=u_title+' Alpha Enhancement', cmap='gnuplot2', vmin=vmin, vmax=vmax,
			flux_unbinned=D.unbinned_flux, signal_noise=D.SNRatio, 
			signal_noise_target=SN_target, center=center,
			save='%s/Alpha_uncert.png'%(out_plots), close=True)

		if gradient:
			r = np.sqrt((D.xBar - center[0])**2 + (D.yBar - center[1])**2)
			for i in ['age', 'met', 'alp']:
				if i=='age': 
					y = np.log10(eval(i))
					y_err =  np.abs(eval('unc_'+i)/np.array(eval(i))/
						np.log(10))
				else: 
					y = eval(i)
					y_err = eval('unc_'+i)
				axs[i].errorbar(r, y, yerr=y_err, fmt='.', c='k')

				params, cov = np.polyfit(r, y, 1, w=1/y_err, cov=True)
				axs[i].plot(r, np.poly1d(params)(r), '--k')
				# params, residuals, _, _, _ = numpy.polyfit(r, y, 1, w=1/y_err, 
				# 	full=True)
				# chi2 = residuals / (len(r) - 2)
				figs[i].text(0.15, 0.84, r'grad: %.3f $\pm$ %.3f'%(params[0], 
					np.sqrt(np.diag(cov))[0]))

	if gradient:
		out_plots = "%s/plots/population" % (output)

		index = np.zeros((150,150,2))
		for i in range(index.shape[0]):
			for j in range(index.shape[1]):
				index[i,j,:] = np.array([i,j]) - center

		step_size = 12
		annuli = np.arange(step_size, 100, step_size).astype(float)

		age_rad = np.zeros(len(annuli))
		met_rad = np.zeros(len(annuli))
		alp_rad = np.zeros(len(annuli))

		age_err_rad = np.zeros(len(annuli))
		met_err_rad = np.zeros(len(annuli))
		alp_err_rad = np.zeros(len(annuli))

		for i, a in enumerate(annuli):
			params = set_params(reps=0, opt='pop', gas=1, produce_plot=False)

			mask = (np.sqrt(index[:,:,0]**2 + index[:,:,1]**2) < a) * (
				np.sqrt(index[:,:,0]**2 + index[:,:,1]**2) > a - step_size)

			spec = np.nansum(f[1].data[:,mask], axis=1)
			noise = np.sqrt(np.nansum(f[2].data[:,mask]**2, axis=1))

			lam = np.arange(len(spec))*header['CD3_3'] + header['CRVAL3']
			spec, lam, cut = apply_range(spec, lam=lam, 
				set_range=params.set_range, return_cuts=True)
			lamRange = np.array([lam[0],lam[-1]])
			noise = noise[cut]

			pp = run_ppxf(galaxy, spec, noise, lamRange, header['CD3_3'], 
				params)

			pop = population(pp=pp, instrument='muse', method=method)

			for i in ['age', 'met', 'alp']:
				if i=='met': i2='metallicity'
				elif i=='alp': i2 = 'alpha'
				else: i2=i
				rad[i].append(getattr(pop, i2))
				rad_err[i].append(getattr(pop, 'unc_'+i))

		annuli *= abs(header['CD1_1'])*(60**2)


		gradient_file = '%s/galaxies_pop_gradients.txt' % (out_dir)
		ageRe, ageG, e_ageG, metRe, metG, e_metG, alpRe, alpG, e_alpG = \
			np.loadtxt(gradient_file, usecols=(1,2,3,4,5,6,7,8,9), 
			unpack=True, skiprows=1)
		galaxy_gals = np.loadtxt(gradient_file, usecols=(0,), unpack=True, 
			skiprows=1, dtype=str)
		i_gal = np.where(galaxy_gals == galaxy)[0][0]

		R_e = get_R_e(galaxy)

		for i in ['age', 'met', 'alp']:
			axs[i].set_xlabel('Radius (arcsec)')

			if i=='age': 
				y = np.log10(rad[i])
				y_err = np.abs(np.array(rad_err[i])/np.array(rad[i])/
					np.log(10))
			else: 
				y = np.array(rad[i])
				y_err = np.array(rad_err[i])
			axs[i].errorbar(annuli, y, yerr=y_err, 
				fmt='x', c='r')

			params, cov = np.polyfit(annuli, y, 1, w=1/y_err, cov=True)
			axs[i].plot(annuli, np.poly1d(params)(annuli), 
				'-r')
			# params, residuals, _, _, _ = numpy.polyfit(annuli, y, 1, 
			# 	w=1/y_err, full=True)
			# chi2 = residuals / (len(annuli) - 2)
			figs[i].text(0.15, 0.8, r'grad: %.3f $\pm$ %.3f'%(params[0], 
				np.sqrt(np.diag(cov))[0]), color='r')

			if i =='age':
				axs[i].set_ylabel('log(Age (Gyr))')

				ageG[i_gal], e_ageG[i_gal] = params[0],np.sqrt(np.diag(cov))[0]
				ageRe[i_gal] = np.poly1d(params)(R_e)
			elif i == 'met':
				axs[i].set_ylabel('Metalicity [Z/H]')

				metG[i_gal], e_metG[i_gal] = params[0],np.sqrt(np.diag(cov))[0]
				metRe[i_gal] = np.poly1d(params)(R_e)
			elif i == 'alp':
				axs[i].set_ylabel('Alpha Enhancement [alpha/Fe]')

				alpG[i_gal], e_alpG[i_gal] = params[0],np.sqrt(np.diag(cov))[0]
				alpRe[i_gal] = np.poly1d(params)(R_e)
			figs[i].savefig('%s/%s_grad.png' % (out_plots, i))
			plt.close(i)

		temp = "{0:12}{1:7}{2:7}{3:7}{4:7}{5:7}{6:7}{7:7}{8:7}{9:7}\n"
		with open(gradient_file, 'w') as f:
			f.write(temp.format('Galaxy', 'ageRe', 'ageG', 'e_ageG', 'metRe', 
				'metG', 'e_metG', 'alpRe', 'alpG', 'e_alpG'))
			for i in range(len(galaxy_gals)):
				f.write(temp.format(galaxy_gals[i], str(round(ageRe[i],1)),
					str(round(ageG[i],3)), str(round(e_ageG[i],3)), 
					str(round(metRe[i],1)), str(round(metG[i],3)), 
					str(round(e_metG[i],3)), str(round(alpRe[i],1)),
					str(round(alpG[i],3)), str(round(e_alpG[i],3))))

	return D








##############################################################################

# Use of plot_stellar_pop.py

if __name__ == '__main__':
	for g in ['ic1459','ic4296', 'ngc1316','ngc1399']:
		plot_stellar_pop(g, gradient='only', method='mostlikely')

