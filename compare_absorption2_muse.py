# New routine to compare MUSE observations with Rampazzo

from checkcomp import checkcomp
cc=checkcomp()
if cc.remote:
	import matplotlib # 20160202 JP to stop lack-of X-windows error
	matplotlib.use('Agg') # 20160202 JP to stop lack-of X-windows error
	import matplotlib.pyplot as plt # used for plotting
else:
	import matplotlib.pyplot as plt # used for plotting
import numpy as np 
from prefig import Prefig
Prefig()
import cPickle as pickle
from astropy.io import fits
from errors2_muse import get_dataCubeDirectory
from classify import get_R_e
from errors2_muse import run_ppxf, set_params,apply_range
from glob import glob
from pop import get_absorption


def compare_absortion(galaxy):
	f = fits.open(get_dataCubeDirectory(galaxy))
	header = f[1].header
	s = f[1].data.shape

	lines = ['H_beta', 'Fe5015', 'Mg_b', #'Fe5270', 'Fe5335', 
		'Fe5406', 'Fe5709', 'Fe5782', 
		'NaD', 'TiO1']
	color = ['purple',   'k',  'orange',   #'g',       'b',      
		'c',  'lightblue',  'grey',
		 'r',  'gold']

	R_e = get_R_e(galaxy)
	apertures = np.array([1.5, 2.5, 10, R_e/10, R_e/8, R_e/4, R_e/2]) # arcsec

	data_file =  "%s/Data/vimos/analysis/galaxies.txt" % (cc.base_dir)
	file_headings = np.loadtxt(data_file, dtype=str)[0]
	z_gals = np.loadtxt(data_file, unpack=True, skiprows=1, usecols=(1,), dtype=float)
	galaxy_gals = np.loadtxt(data_file, skiprows=1, usecols=(0,),dtype=str)
	i_gal = np.where(galaxy_gals==galaxy)[0][0]
	z = z_gals[i_gal]

	data_file =  "%s/Data/muse/analysis/galaxies.txt" % (cc.base_dir)
	x_cent_gals, y_cent_gals = np.loadtxt(data_file, unpack=True, skiprows=1, 
		usecols=(1,2), dtype=int)
	galaxy_gals = np.loadtxt(data_file, skiprows=1, usecols=(0,),dtype=str)
	i_gal = np.where(galaxy_gals==galaxy)[0][0]	
	center = np.array([x_cent_gals[i_gal], y_cent_gals[i_gal]])

	index = np.zeros((s[1],s[2],2))
	for i in range(s[1]):
		for j in range(s[2]):
			index[i,j,:] = np.array([i,j]) - center

	params = set_params(reps=0, opt='pop')
	fig, ax = plt.subplots()
	for a in apertures:
		mask = np.sqrt(index[:,:,0]**2 + index[:,:,1]**2) * header['CD3_3'] < a

		spec = np.nansum(f[1].data[:,mask], axis=1)
		noise = np.sqrt(np.nansum(f[2].data[:,mask]**2, axis=1))
		lam = np.arange(len(spec))*header['CD3_3'] + header['CRVAL3']
		spec, lam, cut = apply_range(spec, lam=lam, set_range=params.set_range, 
			return_cuts=True)
		lamRange = np.array([lam[0],lam[-1]])
		noise = noise[cut]

		pp = run_ppxf(galaxy, spec, noise, lamRange, header['CD3_3'], params, 
			produce_plot=False)

		absorp, uncert = get_absorption(lines, pp=pp, instrument='muse')
		
		for i, l in enumerate(lines):
			ax.errorbar(a, absorp[l], yerr=uncert[l], color=color[i], fmt='x')
	# for i, l in enumerate(lines):
	# 	ax.errorbar(np.nan, np.nan, color=color[i], fmt='x', label=l)
	ax.legend(facecolor='w')


	# Rampazzo_file = '%s/Data/lit_absorption/Rampazzo_aperture.txt' % (cc.base_dir)
	Rampazzo_file = '%s/Data/lit_absorption/J_A+A_433_497_table9.txt' % (cc.base_dir)
	file_headings = np.loadtxt(Rampazzo_file, dtype=str)[0]

	Rampazzo_corr_file = '%s/Data/lit_absorption/Rampazzo_lick_corrections.txt' % (
		cc.base_dir)
	alpha, beta = np.loadtxt(Rampazzo_corr_file, unpack=True, usecols=(1,2), skiprows=1)
	index_name = np.loadtxt(Rampazzo_corr_file, unpack=True, usecols=(0,), skiprows=1, 
		dtype=str)
	for i, l in enumerate(lines):
		col = np.where(file_headings==l)[0][0]

		try:
			col2 = np.where(file_headings==l)[0][1]
		except:
			try:
				col2 = np.where(file_headings=='_'+l)[0][0]
			except:
				col2 = np.where(file_headings=='e_'+l)[0][0]
		R_obs, R_err = np.loadtxt(Rampazzo_file, unpack=True, skiprows=5, 
			usecols=(col,col2))
		R_galaxies = np.loadtxt(Rampazzo_file, unpack=True, skiprows=5, usecols=(0,), 
			dtype=str)

		mask = R_galaxies==galaxy.upper()

		order = np.argsort(apertures)
		a = alpha[index_name==l]
		b = beta[index_name==l]

		ax.errorbar(apertures[order], (R_obs[mask][order]-b)/a, yerr=R_err[mask][order]/a, 
			color=color[i])

	fig.savefig('%s/Data/lit_absorption/Rampazzo_aperture_muse_%s.png' % (cc.base_dir, 
		galaxy))

if __name__=='__main__':
	for g in ['ic1459','ic4296']:
		print g
		compare_absortion(g)