from checkcomp import checkcomp
cc = checkcomp()
if 'home' not in cc.device:
	import matplotlib # 20160202 JP to stop lack-of X-windows error
	matplotlib.use('Agg') # 20160202 JP to stop lack-of X-windows error
import cPickle as pickle
import matplotlib.pyplot as plt
import numpy as np
from plot_velfield_nointerp import plot_velfield_nointerp
from plot_results import add_
from fits.io import fits
from errors2_muse import get_dataCubeDirectory
import os



def BPT(galaxy, D=None, opt='kin'):
	print '   BPT'

	analysis_dir = "%s/Data/muse/analysis" % (cc.base_dir)
	galaxiesFile = "%s/galaxies.txt" % (analysis_dir)
	x_cent_gals, y_cent_gals = np.loadtxt(galaxiesFile, unpack=True, skiprows=1, 
		usecols=(1,2), dtype=int)
	galaxy_gals = np.loadtxt(galaxiesFile, skiprows=1, usecols=(0,),dtype=str)
	i_gal = np.where(galaxy_gals==galaxy)[0][0]
	center = (x_cent_gals[i_gal], y_cent_gals[i_gal])

	output = '%s/%s/%s' % (analysis_dir, galaxy, opt)


	if D is None:
		pickleFile = open('%s/pickled/dataObj.pkl' % (output))
		D = pickle.load(pickleFile)
		pickleFile.close()

	for bin in D.bin:
		for e in bin.e_line.itervalues():
			e.__threshold__ = 0
# ------------=============== BPT diagram =================----------
	fig, ax = plt.subplots()
	y = np.log(D.e_line['[OIII]5007d'].flux/D.e_line['Hbeta'].flux)
	x = np.log(D.e_line['[SII]6716'].flux/D.e_line['Halpha'].flux)

	y_err = np.sqrt((D.e_line['[OIII]5007d'].flux.uncert/D.e_line['[OIII]5007d'].flux)**2 +
		(D.e_line['Hbeta'].flux.uncert/D.e_line['Hbeta'].flux)**2)
	x_err = np.sqrt((D.e_line['[SII]6716'].flux.uncert/D.e_line['Halpha'].flux)**2 +
		(D.e_line['Hbeta'].flux.uncert/D.e_line['Hbeta'].flux)**2)

	Seyfert2 = (0.72/(x - 0.32) + 1.30 < y) * (1.89 * x + 0.76 < y)
	LINER = ((0.72/(x - 0.32) + 1.30 < y) * (1.89 * x + 0.76 > y))# + (x > 0.32)
	SF = (y < 0.72/(x - 0.32) + 1.30)# - (x > 0.32)
	# o = (1.89 * x + 0.76 > y) * (x > 0.32)

	# ax.scatter(x[Seyfert2], y[Seyfert2], c='r')
	# ax.scatter(x[LINER], y[LINER], c='g')
	# ax.scatter(x[SF], y[SF], c='b')
	ax.errorbar(x[Seyfert2], y[Seyfert2], yerr=y_err[Seyfert2], xerr=x_err[Seyfert2], c='r',
		fmt='.')
	ax.errorbar(x[LINER], y[LINER], yerr=y_err[LINER], xerr=x_err[LINER], c='g', fmt='.')
	ax.errorbar(x[SF], y[SF], yerr=y_err[SF], xerr=x_err[SF], c='b', fmt='.')
	# ax.scatter(x[o],y[o],c='k')

	x_line1 = np.arange(-2.0, 0.1, 0.001)
	y_line1 = 0.72/(x_line1 - 0.32) + 1.30
	ax.plot(x_line1, y_line1,'k')


	y_line2 = 1.89 * x_line1 + 0.76

	m = y_line2 > y_line1
	a = np.min(x_line1[m])
	x_line2 = np.arange(a, 0.60, 0.001)
	y_line2 = 1.89 * x_line2 + 0.76

	ax.plot(x_line2, y_line2, 'k')

	ax.set_xlim([-1.2, 0.7])
	ax.set_ylim([-1.2, 1.5])

	ax.set_ylabel(r'log([OIII]/$H_\beta$)')
	ax.set_xlabel(r'log([SII]/$H_\alpha$)')

	saveTo = '%s/plots/BPT.png' % (output)
	if not os.path.exists(os.path.dirname(saveTo)):
		os.makedirs(os.path.dirname(saveTo))  
	fig.savefig(saveTo)
# ------------================= BPT map ===================----------
	m = np.ones(D.number_of_bins)*np.nan
	m[Seyfert2] = +1
	m[LINER] = 0
	m[SF] = -1

	saveTo = '%s/plots/AGN_map.png' % (output)
	ax = plot_velfield_nointerp(D.x, D.y, D.bin_num, D.xBar, D.yBar, m, vmin=-1.3, vmax=1.3,
		nodots=True, title='Map of AGN', save=saveTo, res=0.2, 
		flux_unbinned=D.unbinned_flux, center=center)

	f = fits.open(get_dataCubeDirectory(galaxy))
	header = f[1].header
	f.close()
	add_('radio', 'r', ax, galaxy, header)

	return D


if __name__=='__main__':
	BPT('ic1459')