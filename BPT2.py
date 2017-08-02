# BPT diagram and AGN map diagnotics based on ideas from 
# Davies, R.L. et.al. 2016 MNRAS 462 1616D. 


from checkcomp import checkcomp
cc = checkcomp()
if 'home' not in cc.device:
	import matplotlib # 20160202 JP to stop lack-of X-windows error
	matplotlib.use('Agg') # 20160202 JP to stop lack-of X-windows error
import cPickle as pickle
import matplotlib.pyplot as plt
import numpy as np
from plot_velfield_nointerp import plot_velfield_nointerp
from plot_results_muse import add_, set_lims
from astropy.io import fits
from errors2_muse import get_dataCubeDirectory
import os
from prefig import Prefig


def BPT(galaxy, D=None, opt='kin'):
	print '   BPT'
	Prefig(size=(16*3,12), transparent=False)

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
	# fig, ax = plt.subplots(1,3, sharey=True)
	for i, l in enumerate(['[NII]6583d', '[SII]6716', '[OI]6300d']):

		y = np.log(D.e_line['[OIII]5007d'].flux/D.e_line['Hbeta'].flux)
		x = np.log(D.e_line[l].flux/D.e_line['Halpha'].flux)

		y_err = np.sqrt((D.e_line['[OIII]5007d'].flux.uncert/
			D.e_line['[OIII]5007d'].flux)**2 +
			(D.e_line['Hbeta'].flux.uncert/D.e_line['Hbeta'].flux)**2)
		x_err = np.sqrt((D.e_line[l].flux.uncert/D.e_line[l].flux)**2 +
			(D.e_line['Halpha'].flux.uncert/D.e_line['Halpha'].flux)**2)

		large_err = (x_err > 0.5) + (y_err > 0.5)


		Seyfert2_combined = np.ones(len(x)).astype(bool)
		LINER_combined = np.ones(len(x)).astype(bool)
		SF_combined = np.ones(len(x)).astype(bool)
		x_line1 = np.arange(-2.2, 1, 0.001)
		if l == '[SII]6716':
			pass
		# 	Seyfert2 = ((0.72/(x - 0.32) + 1.30 < y) + (x > 0.32)) * (
		# 		1.89 * x + 0.76 < y) * ~large_err
		# 	LINER = ((0.72/(x - 0.32) + 1.30 < y) + (x > 0.32)) * (
		# 		1.89 * x + 0.76 > y) * ~large_err
		# 	SF = (y < 0.72/(x - 0.32) + 1.30) * (x < 0.32) * ~large_err

		# 	y_line1 = 0.72/(x_line1 - 0.32) + 1.30
		# 	m = y_line1 < 1
		# 	ax[i].plot(x_line1[m], y_line1[m],'k')

		# 	y_line2 = 1.89 * x_line1 + 0.76
		# 	m = y_line2 > y_line1
		# 	a = np.min(x_line1[m])
		# 	x_line2 = np.arange(a, 0.60, 0.001)
		# 	y_line2 = 1.89 * x_line2 + 0.76
		# 	ax[i].plot(x_line2, y_line2, 'k')

		# 	lab = '[SII]'

		# 	ax[i].set_xlim([-1.2, 0.7])


		elif l == '[NII]6583d':
			distance = np.zeros(len(x))
			Seyfert2 = ((0.61/(x - 0.47) + 1.19 < y) + (x > 0.47)) * ~large_err
			LINER = ((0.61/(x - 0.47) + 1.19 < y) + (x > 0.47)) * ~large_err
			SF = (0.61/(x - 0.47) + 1.19 > y) * (x < 0.47) * ~large_err

			y_line1 = 0.61/(x_line1 - 0.47) + 1.19
			m = y_line1 < 1
			# ax[i].plot(x_line1[m], y_line1[m],'k')

			for i in range(len(distance)):
				distance[i] = np.sqrt(np.min((x_line1[m] - x[i])**2 + 
					(y_line1[m] - y[i])**2))
			limit = 1 # if distance is greater than limit then consider it completely 
					#	in it's region.
			distance[distance > limit] = limit
			distance[SF] *= -1

			lab = '[NII]'

			# ax[i].set_xlim([-2, 1])

		elif l == '[OI]6300d':
			pass
		# 	Seyfert2 = ((y > 0.73/(x + 0.59) + 1.33) + (x > -0.59)) * (
		# 		y > 1.18 * x + 1.30) * ~large_err
		# 	LINER = ((y > 0.73/(x + 0.59) + 1.33) + (x > -0.59)) * (
		# 		y < 1.18 * x + 1.30) * ~large_err
		# 	SF = (y < 0.73/(x + 0.59) + 1.33) * (x < -0.59) * ~large_err

		# 	y_line1 = 0.73/(x_line1 + 0.59) + 1.33
		# 	m = y_line1 < 1
		# 	ax[i].plot(x_line1[m], y_line1[m],'k')

		# 	y_line2 = 1.18 * x_line1 + 1.30
		# 	m = y_line2 > y_line1
		# 	a = np.min(x_line1[m])
		# 	x_line2 = np.arange(a, 0.60, 0.001)
		# 	y_line2 = 1.18 * x_line2 + 1.30
		# 	ax[i].plot(x_line2, y_line2, 'k')

		# 	lab = '[OI]'

		# 	ax[i].set_xlim([-2.2, 0])

		# ax[i].set_ylim([-1.2, 1.5])

		Seyfert2_combined *= Seyfert2
		LINER_combined *= LINER
		SF_combined *= SF

		# ax[i].errorbar(x[Seyfert2], y[Seyfert2], yerr=y_err[Seyfert2], 
		# 	xerr=x_err[Seyfert2], c='r', fmt='.')
		# ax[i].errorbar(x[LINER], y[LINER], yerr=y_err[LINER], xerr=x_err[LINER], c='g', 
		# 	fmt='.')
		# ax[i].errorbar(x[SF], y[SF], yerr=y_err[SF], xerr=x_err[SF], c='b', fmt='.')

		# ax[0].set_ylabel(r'log([OIII]/$H_\beta$)')
		# ax[i].set_xlabel(r'log(%s/$H_\alpha$)' %( lab))

	# saveTo = '%s/plots/BPT.png' % (output)
	# if not os.path.exists(os.path.dirname(saveTo)):
	# 	os.makedirs(os.path.dirname(saveTo))  
	# fig.savefig(saveTo)
	# plt.close()
# ------------================= BPT map ===================----------

	f = fits.open(get_dataCubeDirectory(galaxy))
	header = f[1].header
	f.close()
	saveTo = '%s/plots/AGN_map2.png' % (output)
	ax = plot_velfield_nointerp(D.x, D.y, D.bin_num, D.xBar, D.yBar, distance, header, 
		vmin=-1.3, vmax=1.3, nodots=True, title='Map of AGN', save=saveTo, res=0.2, 
		flux_unbinned=D.unbinned_flux, center=center)
	ax.saveTo = saveTo
	add_('radio', 'r', ax, galaxy)
	plt.close()
# ------------============== WHaN2 diagram ================----------
	Prefig(size=(16,12*3), transparent=False)
	fig, ax = plt.subplots(3)
	x = np.log10(D.e_line['[NII]6583d'].flux/D.e_line['Halpha'].flux)
	y = np.log10(D.e_line['Halpha'].equiv_width)

	Ha_Hapeak = D.e_line['Halpha'].flux/np.max(D.e_line['Halpha'].flux)
	p1 = Ha_Hapeak <= 0.2
	p2 = Ha_Hapeak > 0.2
	c = np.array(np.sqrt((D.xBar-center[0])**2 +(D.yBar-center[1])**2))

	passive = (D.e_line['Halpha'].equiv_width<0.5) * (
		D.e_line['[NII]6583d'].equiv_width<0.5)

	ax[0].scatter(x, y,c=c)
	ax[0].scatter(x[passive],y[passive],marker='x',c='r')


	xlim = ax[0].get_xlim()
	ylim = ax[0].get_ylim()
	xlim = [-2,0.75]
	ylim = [-1,1.2]

	ax[1].scatter(x[p1], y[p1],c=c[p1], vmin=np.min(c), vmax=np.max(c))
	ax[1].scatter(x[p1*passive],y[p1*passive],marker='x',c='r')
	ax[1].text(-1.9,-0.8,r'0 $<$ H$_\alpha$/H$_\alpha^{\mathrm{peak}}$ $\leq$ 0.2')

	ax[2].scatter(x[p2], y[p2],c=c[p2], vmin=np.min(c), vmax=np.max(c))
	ax[2].scatter(x[p2*passive],y[p2*passive],marker='x',c='r')
	ax[2].text(-1.9,-0.8, r'0.2 $<$ H$_\alpha$/H$_\alpha^{\mathrm{peak}}$ $\leq$ 1')
	

	# Dianotics lines from R Cid Fernandes et.al. 2011 MNRAS 413 1687
	for a in ax:
		a.axhline(np.log10(3), ls=':', c='k') # 
		a.plot([-0.4,-0.4],[np.log10(3),ylim[1]], ls=':', c='k')
		a.plot([-0.4,xlim[1]],[np.log10(6),np.log10(6)], ls=':', c='k')
		a.set_xlim(xlim)
		a.set_ylim(ylim)
		a.set_ylabel(r'log(EW(H$_\alpha$)/$\AA$)')
		a.set_xlabel(r'log([NII]/H$_\alpha$)')
		a.text(-1.9,1,'Star Forming')
		a.text(-0.35,1,'strong AGN')
		a.text(-0.35,0.55,'weak AGN')
		a.text(-1.9,0.25,'Retired Galaxies')
	fig.suptitle('WHaN2 plot')

	fig.savefig('%s/plots/WHaN2.png' % (output))

	Prefig(size=(16,12), transparent=False)

	return D


if __name__=='__main__':
	for gal in ['ic1459','ic4296','ngc1316','ngc1399']:
		BPT(gal)