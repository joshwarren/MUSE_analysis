from checkcomp import checkcomp
cc = checkcomp()
if 'home' not in cc.device:
	import matplotlib # 20160202 JP to stop lack-of X-windows error
	matplotlib.use('Agg') # 20160202 JP to stop lack-of X-windows error
import cPickle as pickle
import matplotlib.pyplot as plt
import numpy as np
from plot_velfield_nointerp import plot_velfield_nointerp
from plot_results_muse import add_
from astropy.io import fits
from errors2_muse import get_dataCubeDirectory
import os
from prefig import Prefig
from sauron_colormap import sauron



def BPT(galaxy, D=None, opt='kin', norm="lwv"):
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

	if D.norm_method != norm:
		D.norm_method = norm
		D.find_restFrame()

	# D.__threshold__ = 3
# ------------=============== BPT diagram =================----------
	if all([l in D.e_components for l in ['[NII]6583d', '[SII]6716', 
		'[OI]6300d', 'Heta', 'Halpha', '[OIII]5007d']]):
	
		fig, ax = plt.subplots(1,3, sharey=True)
		Prefig(size=(16,12))
		fig2, ax2 = plt.subplots()
		for i, l in enumerate(['[NII]6583d', '[SII]6716', '[OI]6300d']):

			y = np.log10(D.e_line['[OIII]5007d'].flux/D.e_line['Hbeta'].flux)
			x = np.log10(D.e_line[l].flux/D.e_line['Halpha'].flux)

			y_err = np.sqrt((D.e_line['[OIII]5007d'].flux.uncert/
				D.e_line['[OIII]5007d'].flux)**2 +
				(D.e_line['Hbeta'].flux.uncert/D.e_line['Hbeta'].flux)**2)/np.log(10)
			x_err = np.sqrt((D.e_line[l].flux.uncert/D.e_line[l].flux)**2 +
				(D.e_line['Halpha'].flux.uncert/D.e_line['Halpha'].flux)**2)/np.log(10)

			large_err = (x_err > 0.5) + (y_err > 0.5)


			Seyfert2_combined = np.ones(len(x)).astype(bool)
			LINER_combined = np.ones(len(x)).astype(bool)
			SF_combined = np.ones(len(x)).astype(bool)
			x_line1 = np.arange(-2.2, 1, 0.001)
			if l == '[SII]6716':
				Seyfert2 = ((0.72/(x - 0.32) + 1.30 < y) + (x > 0.32)) * (1.89 * x + 0.76 < y) * ~large_err
				LINER = ((0.72/(x - 0.32) + 1.30 < y) + (x > 0.32)) * (1.89 * x + 0.76 > y) * ~large_err
				SF = (y < 0.72/(x - 0.32) + 1.30) * (x < 0.32) * ~large_err

				y_line1 = 0.72/(x_line1 - 0.32) + 1.30
				m = y_line1 < 1
				ax[i].plot(x_line1[m], y_line1[m],'k')

				y_line2 = 1.89 * x_line1 + 0.76
				m = y_line2 > y_line1
				a = np.min(x_line1[m])
				x_line2 = np.arange(a, 0.60, 0.001)
				y_line2 = 1.89 * x_line2 + 0.76
				ax[i].plot(x_line2, y_line2, 'k')

				lab = '[SII]'

				ax[i].set_xlim([-1.2, 0.7])


			elif l == '[NII]6583d':
				Seyfert2 = ((0.61/(x - 0.47) + 1.19 < y) + (x > 0.47)) * ~large_err
				LINER = ((0.61/(x - 0.47) + 1.19 < y) + (x > 0.47)) * ~large_err
				SF = (0.61/(x - 0.47) + 1.19 > y) * (x < 0.47) * ~large_err

				y_line1 = 0.61/(x_line1 - 0.47) + 1.19
				m = y_line1 < 1
				ax[i].plot(x_line1[m], y_line1[m],'k')
				ax2.plot(x_line1[m], y_line1[m],'k')

				lab = '[NII]'

				ax[i].set_xlim([-2, 1])
				ax2.set_xlim([-2, 1])
				ax2.set_ylim([-1.2, 1.5])

				distance = np.zeros(len(x))
				for j in range(len(distance)):
					distance[j] = np.sqrt(np.min((x_line1[m] - x[j])**2 + 
						(y_line1[m] - y[j])**2))
				limit = 1 # if distance is greater than limit then consider it completely 
						#	in it's region.
				distance[distance > limit] = limit
				distance[SF] *= -limit

				try:
					# color = sauron((distance - min(distance))/(max(distance) - min(distance)))
					# ax2.errorbar(x, y, yerr=y_err, xerr=x_err, c=color, fmt='.')
					ax2.scatter(x, y, c=distance, cmap=sauron, vmin=-limit, vmax=limit)
				except:
					print 'This did not work cotton. Galaxy: %s'%(galaxy)

				ax2.set_ylabel(r'log([OIII]/$H_\beta$)')
				ax2.set_xlabel(r'log(%s/$H_\alpha$)' % (lab))

			elif l == '[OI]6300d':
				Seyfert2 = ((y > 0.73/(x + 0.59) + 1.33) + (x > -0.59)) * (y > 1.18 * x + 1.30) * ~large_err
				LINER = ((y > 0.73/(x + 0.59) + 1.33) + (x > -0.59)) * (y < 1.18 * x + 1.30) * ~large_err
				SF = (y < 0.73/(x + 0.59) + 1.33) * (x < -0.59) * ~large_err

				y_line1 = 0.73/(x_line1 + 0.59) + 1.33
				m = y_line1 < 1
				ax[i].plot(x_line1[m], y_line1[m],'k')

				y_line2 = 1.18 * x_line1 + 1.30
				m = y_line2 > y_line1
				a = np.min(x_line1[m])
				x_line2 = np.arange(a, 0.60, 0.001)
				y_line2 = 1.18 * x_line2 + 1.30
				ax[i].plot(x_line2, y_line2, 'k')

				lab = '[OI]'

				ax[i].set_xlim([-2.2, 0])

			ax[i].set_ylim([-1.2, 1.5])

			Seyfert2_combined *= Seyfert2
			LINER_combined *= LINER
			SF_combined *= SF

			
			ax[i].errorbar(x[LINER], y[LINER], yerr=y_err[LINER], xerr=x_err[LINER], c='g', 
				fmt='.')
			ax[i].errorbar(x[Seyfert2], y[Seyfert2], yerr=y_err[Seyfert2], 
				xerr=x_err[Seyfert2], c='r', fmt='.')
			ax[i].errorbar(x[SF], y[SF], yerr=y_err[SF], xerr=x_err[SF], c='b', fmt='.')

			ax[0].set_ylabel(r'log([OIII]/$H_\beta$)')
			ax[i].set_xlabel(r'log(%s/$H_\alpha$)' %( lab))

		saveTo = '%s/plots/BPT.png' % (output)
		if not os.path.exists(os.path.dirname(saveTo)):
			os.makedirs(os.path.dirname(saveTo))  
		fig.savefig(saveTo)
		fig2.savefig('%s/plots/BPT2.png' % (output))
		plt.close()
		Prefig(size=(16,12), transparent=False)
# ------------================= BPT map ===================----------
		m = np.ones(D.number_of_bins)*np.nan
		m[Seyfert2_combined] = +1
		m[LINER_combined] = 0
		m[SF_combined] = -1

		f = fits.open(get_dataCubeDirectory(galaxy))
		header = f[1].header
		f.close()
		saveTo = '%s/plots/AGN_map.png' % (output)
		ax = plot_velfield_nointerp(D.x, D.y, D.bin_num, D.xBar, D.yBar, m, header, 
			vmin=-1.3, vmax=1.3, nodots=True, title='Map of AGN', save=saveTo, res=0.2, 
			flux_unbinned=D.unbinned_flux, center=center)
		add_('radio', 'r', ax, galaxy)
		plt.close()


# ------------================= BPT map2 ===================----------
		Prefig(size=(16,12), transparent=False)
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
	if all([l in D.e_components for l in ['[NII]6583d', 'Halpha']]):
		Prefig(size=(16,12*3), transparent=False)
		fig, ax = plt.subplots(3)
		x = np.log10(D.e_line['[NII]6583d'].flux/D.e_line['Halpha'].flux)
		y = np.log10(D.e_line['Halpha'].equiv_width)

		Ha_Hapeak = D.e_line['Halpha'].flux/np.nanmax(D.e_line['Halpha'].flux)
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
		plt.close()
# ------------=============== MEx diagram =================----------
	if all([l in D.e_components for l in ['[OIII]5007d', 'Hbeta']]):
		Prefig()
		# from Atlas3D XXXI (Section 6.2.1)
		fig, ax = plt.subplots()
		y = np.log10(D.e_line['[OIII]5007d'].flux/D.e_line['Hbeta'].flux)
		y_err = np.sqrt((D.e_line['[OIII]5007d'].flux.uncert/
			D.e_line['[OIII]5007d'].flux)**2 + (D.e_line['Hbeta'].flux.uncert/
			D.e_line['Hbeta'].flux)**2)/np.log(10)

		large_err = y_err**2 > 1
		m = ~large_err * (D.e_line['[OIII]5007d'].equiv_width < 0.8)
		ax.errorbar(D.components['stellar'].plot['sigma'][m], y[m], c='b',
			xerr = D.components['stellar'].plot['sigma'].uncert[m], yerr=y_err[m], fmt='.',
			label='EW([OIII]) < 0.8')

		m = ~large_err * (D.e_line['[OIII]5007d'].equiv_width >= 0.8)
		ax.errorbar(D.components['stellar'].plot['sigma'][m], y[m], c='r',
			xerr = D.components['stellar'].plot['sigma'].uncert[m], yerr=y_err[m], fmt='.',
			label=r'EW([OIII]) $\ge 0.8$')

		ax.legend(facecolor='w')

		x_line = np.arange(70,1000,1)
		y_line = 1.6*10**-3 * x_line + 0.33
		ax.plot(x_line, y_line, c='k')

		ax.set_xlim([0, min(max(D.components['stellar'].plot['sigma'][m]), 500)])
		ax.set_ylim([-1.2, 1.5])

		ax.axvline(70, c='k')
		ax.axhline(np.log10(0.5), 
			xmin=70./min(max(D.components['stellar'].plot['sigma'][m]), 500), c='k')
		ax.axhline(np.log10(1), 
			xmin=70./min(max(D.components['stellar'].plot['sigma'][m]), 500), c='k')


		ylim = ax.get_ylim()
		yrange = ylim[1] - ylim[0]
		ax.text(60, ylim[0] +0.96 * yrange, 'SF')
		ax.text(75, 0.55, 'Seyfert 2')
		ax.text(75, 0.15, 'LINER')
		ax.text(75, -0.23, 'Transition')

		ax.set_xlabel(r'$\sigma_\ast$')
		ax.set_ylabel(r'log [OIII]d/H$_\beta$')
		ax.set_title('Mass-excitation (MEx) diagnotics for %s' % (galaxy.upper()))

		fig.savefig('%s/plots/MEx.png' % (output))
# ------------============== SAURON diagram ===============----------
	if all([l in D.e_components for l in ['[NI]d', 'Hbeta', 
		'[OIII]5007d']]):

		# from SAURON XVI.
		fig, ax = plt.subplots()
		# y and y_err as MEx above

		x = np.log10(D.e_line['[NI]d'].flux/D.e_line['Hbeta'].flux)
		x_err = np.sqrt((D.e_line['[NI]d'].flux.uncert/
			D.e_line['[NI]d'].flux)**2 + (D.e_line['Hbeta'].flux.uncert/
			D.e_line['Hbeta'].flux)**2)/np.log(10)


		ax.errorbar(x, y, xerr = x_err, yerr = y_err, fmt='.')

		ax.set_xlim([-2.5, 0.5])
		ax.set_ylim([-1.5, 1.5])
		ax.set_xlabel(r'log [NI]d/H$_\beta$')
		ax.set_ylabel(r'log [OIII]d/H$_\beta$')
		ax.set_title('SAURON diagnotics for %s' % (galaxy.upper()))

		fig.savefig('%s/plots/SAURON_diagnoistic.png' % (output))
	return D


if __name__=='__main__':
	for gal in ['ic1459','ic4296','ngc1316','ngc1399']:
		D = BPT(gal)
