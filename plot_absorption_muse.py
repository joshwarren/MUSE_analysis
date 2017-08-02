## ==================================================================
## 		Plot the absorption indices
## ==================================================================
import cPickle as pickle
import matplotlib.pyplot as plt 
from plot_velfield_nointerp import plot_velfield_nointerp 
import numpy as np 
import os
from plot_results_muse import set_lims
from checkcomp import checkcomp
cc = checkcomp()
from errors2_muse import get_dataCubeDirectory
from astropy.io import fits
from prefig import Prefig


def plot_absorption(galaxy, opt='pop', D=None, uncert=True):
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
	out_plots = "%s/plots" % (output)
	if not os.path.exists(out_plots): os.makedirs(out_plots)
	 
	if D is None:
		out_pickle = '%s/pickled' % (output)
		pickleFile = open("%s/dataObj.pkl" % (out_pickle), 'rb')
		D = pickle.load(pickleFile)
		pickleFile.close()

	f = fits.open(get_dataCubeDirectory(galaxy))
	header = f[1].header
	f.close()

	# Set up figure and subplots
	Prefig(size=(16*2,12*np.ceil(len(lines)/2.0)), transparent=False)
	f, ax_array = plt.subplots(int(np.ceil(len(lines)/2.0)), 2, sharex='col', 
		sharey='row', frameon=False)
	if uncert:
		f_uncert, ax_array_uncert = plt.subplots(int(np.ceil(len(lines)/2.0)), 2, 
			sharex='col', sharey='row', frameon=False)

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

		ax_array[int(np.floor(i/2)),i%2] = plot_velfield_nointerp(D.x, D.y, 
			D.bin_num, D.xBar, D.yBar, ab_line, header, 
			vmin=abmin, vmax=abmax, nodots=True, colorbar=True, 
			label='Index strength ('+r'$\AA$'+')', title=titles[line], 
			ax=ax_array[int(np.floor(i/2)),i%2], cmap='gnuplot2', 
			flux_unbinned=D.unbinned_flux, signal_noise=D.SNRatio, 
			signal_noise_target=30)

		if uncert:
			abmin, abmax = set_lims(ab_uncert)

			ax_array_uncert[int(np.floor(i/2)),i%2] = plot_velfield_nointerp(D.x, 
				D.y, D.bin_num, D.xBar, D.yBar, ab_uncert, header, 
				vmin=abmin, vmax=abmax, nodots=True, colorbar=True, 
				label='Index strength ('+r'$\AA$'+')', title=titles[line], 
				ax=ax_array_uncert[int(np.floor(i/2)),i%2], cmap='gnuplot2', 
				flux_unbinned=D.unbinned_flux, signal_noise=D.SNRatio, 
				signal_noise_target=30)


	# f.set_size_inches(8.5,int(np.ceil(len(lines)/2.0))*1.8)

	print 'Saving plot'

	saveTo = "%s/absorption.pdf" % (out_plots)
	# f.tight_layout()
	ax_array[0,1].set_xlabel('')
	ax_array[0,0].set_xlabel('')
	ax_array[0,1].set_ylabel('')
	ax_array[1,1].set_ylabel('')
	f.suptitle(galaxy.upper())
	f.savefig(saveTo)#, bbox_inches="tight")


	if uncert:
		# f_uncert.set_size_inches(8.5,int(np.ceil(len(lines)/2.0))*1.8)

		saveTo = "%s/absorption_uncert.pdf" % (out_plots)
		# f_uncert.tight_layout()
		ax_array_uncert[0,1].set_xlabel('')
		ax_array_uncert[0,0].set_xlabel('')
		ax_array_uncert[0,1].set_ylabel('')
		ax_array_uncert[1,1].set_ylabel('')
		f_uncert.suptitle(galaxy.upper() + ' Uncertainties')
		f_uncert.savefig(saveTo)#, bbox_inches="tight")

	saveTo = '%s/SNR.png' % (out_plots)
	plot_velfield_nointerp(D.x, D.y, D.bin_num, D.xBar, D.yBar, D.SNRatio, header, 
		nodots=True, colorbar=True, label='S/N', show_bin_numbers=True,
		title='Signal to Noise Ratio', flux_unbinned=D.unbinned_flux, save=saveTo, 
		close=True)

	return D








##############################################################################

# Use of plot_absorption.py

if __name__ == '__main__':
	plot_absorption('ic1459')

