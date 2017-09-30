## ==================================================================
## 		Plot the absorption indices
## ==================================================================
import cPickle as pickle
import matplotlib.pyplot as plt 
from plot_velfield_nointerp import plot_velfield_nointerp 
import numpy as np 
import os
from plot_results_muse import set_lims, add_
from checkcomp import checkcomp
cc = checkcomp()
from errors2_muse import get_dataCubeDirectory
from astropy.io import fits
from prefig import Prefig
Prefig()


def plot_absorption(galaxy, opt='pop', D=None, uncert=True, overplot={}):
	# Find lines:
	lines = [#'G4300', 'Fe4383', 'Ca4455', 'Fe4531', 
		'H_beta', 'Fe5015', 
		# 'Mg_1', 'Mg_2', 
		'Mg_b', 'Fe5270', 'Fe5335', 'Fe5406', 'Fe5709', 'Fe5782', 'NaD', 'TiO1', 
		'TiO2']
	# limits = {#'G4300', 'Fe4383', 'Ca4455', 'Fe4531', 
	# 	'H_beta':[1.0,2.9], 'Fe5015':[3.5,5.9], 'Mg_b':[3.1,4.7]}
	titles={'H_beta':r'H$_\beta$', 'Fe5015':'Fe5015', 'Mg_b':r'Mg$_b$', 
		'Fe5270':'Fe5270', 'Fe5335':'Fe5335', 'Fe5406':'Fe5406', 
		'Fe5709':'Fe5709', 'Fe5782':'Fe5782', 'NaD':'NaD', 'TiO1':'TiO1', 
		'TiO2':'TiO2'}

	print 'Absorption lines'

	# Load pickle file from pickler.py
	out_dir = '%s/Data/muse/analysis' % (cc.base_dir)
	output = "%s/%s/%s" % (out_dir, galaxy, opt)
	out_plots = "%s/plots/absorption" % (output)
	if not os.path.exists(out_plots): os.makedirs(out_plots)
	 
	if D is None:
		out_pickle = '%s/pickled' % (output)
		pickleFile = open("%s/dataObj.pkl" % (out_pickle), 'rb')
		D = pickle.load(pickleFile)
		pickleFile.close()

	f = fits.open(get_dataCubeDirectory(galaxy))
	header = f[1].header
	f.close()

	data_file =  "%s/galaxies.txt" % (out_dir)
	# different data types need to be read separetly
	file_headings = np.loadtxt(data_file, dtype=str)[0]
	col = np.where(file_headings=='SN_%s' % (opt))[0][0]

	x_cent_gals, y_cent_gals, SN_target_gals = np.loadtxt(data_file, unpack=True, 
		skiprows=1, usecols=(1,2,col), dtype='int,int,float')
	galaxy_gals = np.loadtxt(data_file, skiprows=1, usecols=(0,),dtype=str)
	i_gal = np.where(galaxy_gals==galaxy)[0][0]
	SN_target=SN_target_gals[i_gal]-10
	center = (x_cent_gals[i_gal], y_cent_gals[i_gal])

	data_file =  "%s/Data/vimos/analysis/galaxies.txt" % (cc.base_dir)
	z_gals = np.loadtxt(data_file, unpack=True, skiprows=1, 
		usecols=(1))
	galaxy_gals = np.loadtxt(data_file, skiprows=1, usecols=(0,),dtype=str)
	i_gal = np.where(galaxy_gals==galaxy)[0][0]
	z = z_gals[i_gal]


	for i, line in enumerate(lines):
		print "    " + line

		if uncert:
			ab_line, ab_uncert = D.absorption_line(line, uncert=True)
		else:
			ab_line = D.absorption_line(line)

		abmin, abmax = set_lims(ab_line, positive=True)

		# if line in limits.keys():
		# 	abmin = limits[line][0]
		# 	abmax = limits[line][1]

		saveTo = '%s/%s.png' % (out_plots, line)
		ax = plot_velfield_nointerp(D.x, D.y, 
			D.bin_num, D.xBar, D.yBar, ab_line, header, 
			vmin=abmin, vmax=abmax, nodots=True, colorbar=True, 
			label='Index strength ('+r'$\AA$'+')', title=titles[line], 
			cmap='gnuplot2', redshift=z, center=center, 
			flux_unbinned=D.unbinned_flux, signal_noise=D.SNRatio, 
			signal_noise_target=SN_target, save=saveTo)
		ax.saveTo = saveTo
		if overplot:
			for o, color in overplot.iteritems():
				add_(o, color, ax, galaxy, close=True)

		if uncert:
			abmin, abmax = set_lims(ab_uncert)

			plot_velfield_nointerp(D.x, D.y, D.bin_num, D.xBar, D.yBar, 
				ab_uncert, header, vmin=abmin, vmax=abmax, nodots=True, 
				colorbar=True, label='Index strength ('+r'$\AA$'+')', 
				title=titles[line], cmap='gnuplot2', 
				flux_unbinned=D.unbinned_flux, signal_noise=D.SNRatio, 
				signal_noise_target=SN_target, close=True,
				save='%s/%s_uncert.png' % (out_plots, line))
	return D








##############################################################################

# Use of plot_absorption.py

if __name__ == '__main__':
	plot_absorption('ic1459')

