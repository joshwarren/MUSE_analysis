## Routine to compare to Atlas3d results.

from checkcomp import checkcomp
cc = checkcomp()
if 'home' not in cc.device:
	import matplotlib # 20160202 JP to stop lack-of X-windows error
	matplotlib.use('Agg') 
import numpy as np 
import matplotlib.pyplot as plt
from markers_atlas3d import marker_atlas3d
from prefig import Prefig
Prefig(transparent=False)

def angle_to_pc(galaxy, angle):
	c = 299792 #km/s
	#H = 67.8 #(km/s)/Mpc # From Planck
	H = 70.0 # value used by Bolonga group.
	file = "%s/Data/vimos/analysis/galaxies.txt" % (cc.base_dir)
	z = np.loadtxt(file, usecols=(1,), skiprows=1, unpack=True)
	gals = np.loadtxt(file, usecols=(0,), skiprows=1, unpack=True, dtype=str)
	i_gal = np.where(gals == galaxy)[0][0]
	z=z[i_gal]
	return np.radians(angle/(60.0*60.0)) * z*c/H*10**6 # pc

def compare_atlas3d():
	print 'Compare to Atlas3d/SAURON'

## ----------============== Ellipticity vs lambda_Re ==============----------
	print 'FR/SR'
	museGalaxiesFile = "%s/Data/muse/analysis/galaxies2.txt" % (cc.base_dir)
	vimosGalaxiesFile = "%s/Data/vimos/analysis/galaxies2.txt" % (cc.base_dir)
	muse_classify_file = "%s/Data/muse/analysis/galaxies_classify.txt" % (cc.base_dir)
	vimos_classify_file = "%s/Data/vimos/analysis/galaxies_classify.txt" % (cc.base_dir)


	lambda_Re_muse, ellipticity_muse =  np.loadtxt(museGalaxiesFile, unpack=True, 
		skiprows=1, usecols=(1,2))
	galaxies_muse =  np.loadtxt(museGalaxiesFile, unpack=True, skiprows=1, usecols=(0,), 
		dtype=str)

	lambda_Re_vimos, ellipticity_vimos =  np.loadtxt(vimosGalaxiesFile, unpack=True, 
		skiprows=1, usecols=(1,2))
	galaxies_vimos =  np.loadtxt(vimosGalaxiesFile, unpack=True, skiprows=1, usecols=(0,), 
		dtype=str)

	fig, ax = plt.subplots()

	# Plot Atlas3d results NB: all atlas3d tables are in alphabetical order
	atlas3d_file = '%s/Data/atlas3d/III_tableB1.dat' % (cc.base_dir)
	ellipticity_atlas, lambda_Re_atlas = np.loadtxt(atlas3d_file, unpack=True, 
		usecols=(2,7), dtype=float)
	atlas3d_file = '%s/Data/atlas3d/II_tableD1.dat' % (cc.base_dir)
	structure_atlas = np.loadtxt(atlas3d_file, unpack=True, usecols=(13,), dtype=str)
	atlas3d_file = '%s/Data/atlas3d/I_table3.dat' % (cc.base_dir)
	T_type = np.loadtxt(atlas3d_file, unpack=True, usecols=(10,))
	E = T_type < -3.5
	
	no_rot_muse = structure_atlas=='NRR/LV'
	complex_rot_muse = structure_atlas=='NRR/NF'
	KDC_atlas = (structure_atlas=='NRR/KDC') + (structure_atlas=='NRR/CRC') + \
		(structure_atlas=='RR/CRC')
	counter_rot_muse = (structure_atlas=='NRR/2s') + (structure_atlas=='RR/2s')
	regular_rot_muse = (structure_atlas=='RR/NF') + (structure_atlas=='RR/2m') + \
		(structure_atlas=='RR/KT')

	# S0s
	ax.scatter(ellipticity_atlas[no_rot_muse*~E], lambda_Re_atlas[no_rot_muse*~E], 
		marker=marker_atlas3d(0), c='lightgrey', alpha=0.5, lw=0,
		label=r'Atlas3D S0: $T>-3.5$')
	ax.scatter(ellipticity_atlas[complex_rot_muse*~E], 
		lambda_Re_atlas[complex_rot_muse*~E], marker=marker_atlas3d(1),
		c='lightgrey', alpha=0.5, lw=0)
	ax.scatter(ellipticity_atlas[KDC_atlas*~E], lambda_Re_atlas[KDC_atlas*~E], 
		marker=marker_atlas3d(2), c='lightgrey', alpha=0.5, lw=0)
	ax.scatter(ellipticity_atlas[counter_rot_muse*~E], 
		lambda_Re_atlas[counter_rot_muse*~E], marker=marker_atlas3d(3),
		c='lightgrey', alpha=0.5, lw=0)
	ax.plot(ellipticity_atlas[regular_rot_muse*~E], 
		lambda_Re_atlas[regular_rot_muse*~E], marker=marker_atlas3d(4),
		c='lightgrey', alpha=0.5, lw=0, markerfacecolor='none')

	# Ellipticals
	ax.scatter(ellipticity_atlas[no_rot_muse*E], lambda_Re_atlas[no_rot_muse*E], 
		marker=marker_atlas3d(0), c='k', alpha=0.5, lw=0, 
		label=r'Atlas3D E: $T \leq -3.5$')
	ax.scatter(ellipticity_atlas[complex_rot_muse*E], 
		lambda_Re_atlas[complex_rot_muse*E], marker=marker_atlas3d(1), c='k', alpha=0.5, 
		lw=0)
	ax.scatter(ellipticity_atlas[KDC_atlas*E], lambda_Re_atlas[KDC_atlas*E], 
		marker=marker_atlas3d(2), c='k', alpha=0.5, lw=0)
	ax.scatter(ellipticity_atlas[counter_rot_muse*E], 
		lambda_Re_atlas[counter_rot_muse*E], marker=marker_atlas3d(3), c='k', alpha=0.5, 
		lw=0)
	ax.plot(ellipticity_atlas[regular_rot_muse*E], 
		lambda_Re_atlas[regular_rot_muse*E], marker=marker_atlas3d(4), c='k', alpha=0.5, 
		lw=0, markerfacecolor='none')
	
	# Join MUSE and VIMOS
	for i_muse, g in enumerate(galaxies_muse):
		if g in galaxies_vimos:
			i_vimos = np.where(galaxies_vimos==g)[0][0]
			if i_muse == 0: # add just one label to legend
				ax.plot([ellipticity_muse[i_muse], ellipticity_vimos[i_vimos]], 
					[lambda_Re_muse[i_muse],lambda_Re_vimos[i_vimos]], 'k--', zorder=1,
					label='same galaxy in MUSE and VIMOS')
			else:
				ax.plot([ellipticity_muse[i_muse], ellipticity_vimos[i_vimos]], 
					[lambda_Re_muse[i_muse],lambda_Re_vimos[i_vimos]], 'k--', zorder=1)

	ell = np.arange(0.01,0.99,0.01)

	# Isotropic line
	# vSigma = 0.831 * np.sqrt(ell/(1 - 0.896 * ell))
	# k = 1.1
	# lambda_R = k*vSigma/np.sqrt(1 + k**2 * vSigma**2)
	# ax.plot(ell, lambda_R, 'g', label='Isotropic Line')

	# Add envolope - from footnote, page 430, SAURON X
	k = 1.1
	alpha = 0.15
	vSigma = np.sqrt((0.09 + 0.1 * ell) * ell/(1 - ell))
	lambda_R = k*vSigma/np.sqrt(1 + k**2 * vSigma**2)
	ax.plot(ell, lambda_R, 'm', label=r'$\delta = 0.7 \epsilon_\mathrm{intr}$')


	# Lines of constant intrinsic ellipticity, from Cappellari ARAA 2016 Review
	i = np.arange(0, np.pi/2, 0.01)
	for j in np.arange(0.2,1.2, 0.2):
		e_ob = ell[int(round(j*(len(ell)-1),0))]
		vSigma_ob = vSigma[int(round(j*(len(ell)-1),0))]

		ell_intr = 1 - np.sqrt(1 + e_ob * (e_ob - 2)*np.sin(i)**2)
		e = np.sqrt(1 - (1 - ell_intr)**2)
		Omega = 0.5*(np.arcsin(e) - e * np.sqrt(1 - e**2))/(
			e * np.sqrt(1 - e**2) - (1 - e**2) * np.arcsin(e))
		delta = 1 - (1 + vSigma_ob**2)/((1 - alpha * vSigma_ob**2) * Omega)
		vSigma_ob = vSigma_ob * np.sin(i)/np.sqrt(1 - delta * np.cos(i)**2)
		lambda_R = k*vSigma_ob/np.sqrt(1 + k**2 * vSigma_ob**2)
		ax.plot(ell_intr, lambda_R, 'k:', linewidth=1)


	# MUSE
	gals_muse2, regular_rot_muse, no_rot_muse, complex_rot_muse, counter_rot_muse, \
		KDC_muse = np.loadtxt(muse_classify_file, unpack=True, usecols=(0,2,3,4,5,6), 
		dtype=str, skiprows=1)
	gal_order = [np.where(galaxies_muse==g)[0][0] for g in gals_muse2]
	regular_rot_muse = (regular_rot_muse!='-')[gal_order]
	no_rot_muse = (no_rot_muse!='-')[gal_order]
	complex_rot_muse = (complex_rot_muse!='-')[gal_order]
	counter_rot_muse = (counter_rot_muse!='-')[gal_order]
	KDC_muse = (KDC_muse!='-')[gal_order]

	ax.scatter(ellipticity_muse[no_rot_muse], lambda_Re_muse[no_rot_muse], 
		marker=marker_atlas3d(0), c='b', alpha=0.5, lw=0, label='MUSE')
	ax.scatter(ellipticity_muse[complex_rot_muse], lambda_Re_muse[complex_rot_muse], 
		marker=marker_atlas3d(1), c='b', alpha=0.5, lw=0)
	ax.scatter(ellipticity_muse[KDC_muse], lambda_Re_muse[KDC_muse], 
		marker=marker_atlas3d(2), c='b', alpha=0.5, lw=0)
	ax.scatter(ellipticity_muse[counter_rot_muse], lambda_Re_muse[counter_rot_muse], 
		marker=marker_atlas3d(3), c='b', alpha=0.5, lw=0)
	ax.plot(ellipticity_muse[regular_rot_muse], lambda_Re_muse[regular_rot_muse], 
		marker=marker_atlas3d(4), c='b', alpha=0.5, lw=0, markerfacecolor='none')


	# VIMOS
	gals_vimos2, regular_rot_vimos, no_rot_vimos, complex_rot_vimos, counter_rot_vimos, \
		KDC_vimos = np.loadtxt(vimos_classify_file, unpack=True, usecols=(0,2,3,4,5,6), 
		dtype=str, skiprows=1)
	gal_order = [np.where(galaxies_vimos==g)[0][0] for g in gals_vimos2]
	regular_rot_vimos = (regular_rot_vimos!='-')[gal_order]
	no_rot_vimos = (no_rot_vimos!='-')[gal_order]
	complex_rot_vimos = (complex_rot_vimos!='-')[gal_order]
	counter_rot_vimos = (counter_rot_vimos!='-')[gal_order]
	KDC_vimos = (KDC_vimos!='-')[gal_order]

	ax.scatter(ellipticity_vimos[no_rot_vimos], lambda_Re_vimos[no_rot_vimos], 
		marker=marker_atlas3d(0), c='r', alpha=0.5, lw=0, label='MUSE')
	ax.scatter(ellipticity_vimos[complex_rot_vimos], lambda_Re_vimos[complex_rot_vimos], 
		marker=marker_atlas3d(1), c='r', alpha=0.5, lw=0)
	ax.scatter(ellipticity_vimos[KDC_vimos], lambda_Re_vimos[KDC_vimos], 
		marker=marker_atlas3d(2), c='r', alpha=0.5, lw=0)
	ax.scatter(ellipticity_vimos[counter_rot_vimos], lambda_Re_vimos[counter_rot_vimos], 
		marker=marker_atlas3d(3), c='r', alpha=0.5, lw=0)
	ax.plot(ellipticity_vimos[regular_rot_vimos], lambda_Re_vimos[regular_rot_vimos], 
		marker=marker_atlas3d(4), c='r', alpha=0.5, lw=0, markerfacecolor='none')

	ax.set_title('Atlas3D Fast/Slow Rotator Classification scheme')
	ax.set_xlabel(r'$\epsilon$')
	ax.set_ylabel(r'$\lambda_R (R_e)$')
	ax.set_xlim([0, 0.9])
	ax.set_ylim([0, 0.8])

	# Plot Slow Rotator bounds
	ax.plot([0,0.4,0.4], [0.08, 0.18, 0], 'k', label='FR/SR boundary')
	ax.legend(facecolor='w')

	# Save plot
	fig.savefig('%s/Data/muse/analysis/lambda_R_ellipticity.png' % (cc.base_dir))
	plt.close()
## ----------================ Core age vs KDC size ================----------
	print 'KDC size/age'
	muse_core_file = "%s/Data/muse/analysis/galaxies_core.txt" % (cc.base_dir)
	muse_classify_file = "%s/Data/muse/analysis/galaxies_classify.txt" % (cc.base_dir)
	vimos_core_file = "%s/Data/vimos/analysis/galaxies_core.txt" % (cc.base_dir)
	vimos_classify_file = "%s/Data/vimos/analysis/galaxies_classify.txt" % (cc.base_dir)
	# vimosGalaxiesFile = "%s/Data/vimos/analysis/galaxies_core.txt" % (cc.base_dir)
	sauron_file = '%s/Data/sauron/VIII_table8.dat' % (cc.base_dir)
	
	fig, ax = plt.subplots()
	# Sauron
	size_sauron, size_unc_sauron, age_sauron, age_unc_sauron = np.loadtxt(sauron_file,
		unpack=True, usecols=(4,5,6,7), skiprows=2) 
	fast_sauron = np.loadtxt(sauron_file, unpack=True, usecols=(8,), skiprows=2, dtype=str)
	fast_sauron = fast_sauron == 'F'

	ax.errorbar(size_sauron[fast_sauron], age_sauron[fast_sauron], fmt='.',
		xerr=size_unc_sauron[fast_sauron], yerr=age_unc_sauron[fast_sauron], color='b',
		label='Fast rotating SAURON')
	ax.errorbar(size_sauron[~fast_sauron], age_sauron[~fast_sauron], fmt='.',
		xerr=size_unc_sauron[~fast_sauron], yerr=age_unc_sauron[~fast_sauron], color='r',
		label='Slow rotating SAURON')

	# MUSE
	age_muse, age_unc_muse, OIII_eqw_muse = np.loadtxt(muse_core_file, unpack=True, 
		usecols=(1,2,7), skiprows=2)
	gals_muse1 = np.loadtxt(muse_core_file, unpack=True, usecols=(0,), 
		skiprows=2, dtype=str)
	gals_muse2, size_muse = np.loadtxt(muse_classify_file, unpack=True, usecols=(0,6), 
		dtype=str, skiprows=1)
	has_KDC = size_muse!='-'
	gals_muse2 = gals_muse2[has_KDC]
	size_muse = size_muse[has_KDC].astype(float)
	size_unc_muse = size_muse*0
	

	for i2, g in enumerate(gals_muse2):
		i1 = np.where(gals_muse1==g)[0][0]
		size_unc_muse[i2] = angle_to_pc(g, 0.5) # estimate size at 0.5" accuracy.
		size_muse[i2] = angle_to_pc(g, size_muse[i2])
		if i2 == 0:
			ax.errorbar(size_muse[i2], age_muse[i1], fmt='.',xerr=size_unc_muse[i2], 
				yerr=age_unc_muse[i1], color='g', label='MUSE')
		else:
			ax.errorbar(size_muse[i2], age_muse[i1], fmt='.',xerr=size_unc_muse[i2], 
				yerr=age_unc_muse[i1], color='g')

	# VIMOS
	age_vimos, age_unc_vimos, OIII_eqw_vimos = np.loadtxt(vimos_core_file, unpack=True, 
		usecols=(1,2,7), skiprows=2)
	gals_vimos1 = np.loadtxt(vimos_core_file, unpack=True, usecols=(0,), 
		skiprows=2, dtype=str)
	gals_vimos2, size_vimos = np.loadtxt(vimos_classify_file, unpack=True, usecols=(0,6), 
		dtype=str, skiprows=1)
	has_KDC = size_vimos!='-'
	gals_vimos2 = gals_vimos2[has_KDC]
	size_vimos = size_vimos[has_KDC].astype(float)
	size_unc_vimos = size_vimos*0
	
	count = 0
	for i2, g in enumerate(gals_vimos2):
		i1 = np.where(gals_vimos1==g)[0][0]
		size_unc_vimos[i2] = angle_to_pc(g, 0.5) # estimate size at 0.5" accuracy.
		size_vimos[i2] = angle_to_pc(g, size_vimos[i2])
		if i2 == 0:
			ax.errorbar(size_vimos[i2], age_vimos[i1], fmt='.',xerr=size_unc_vimos[i2], 
				yerr=age_unc_vimos[i1], color='c', label='VIMOS')
		else:
			ax.errorbar(size_vimos[i2], age_vimos[i1], fmt='.',xerr=size_unc_vimos[i2], 
				yerr=age_unc_vimos[i1], color='c')
		if g in gals_muse2 and count==0:
			count += 1
			i_muse1 = np.where(gals_muse1==g)[0][0]
			i_muse2 = np.where(gals_muse2==g)[0][0]
			ax.plot([size_vimos[i2], size_muse[i_muse2]], 
				[age_vimos[i1], age_muse[i_muse1]], 'k--', 
				label='same galaxy in MUSE and VIMOS')
		elif g in gals_muse2:
			count += 1
			i_muse1 = np.where(gals_muse1==g)[0][0]
			i_muse2 = np.where(gals_muse2==g)[0][0]
			ax.plot([size_vimos[i2], size_muse[i_muse2]], 
				[age_vimos[i1], age_muse[i_muse1]], 'k--', 
				label='same galaxy in MUSE and VIMOS')

	ax.legend(facecolor='w')
	ax.set_yscale('log')#, nonposy='clip', subsy=[1,2,3,4,5,6,7,8,9])
	ax.set_xlabel('KDC size (pc)')
	ax.set_ylabel('Age of central 1 arcsec (Gyrs)')
	ax.set_title('Age and size of KDCs')


	fig.savefig('%s/Data/muse/analysis/KDC_size_age.png' % (cc.base_dir))	
	plt.close()
## ----------=========== Core OIII vs radio power ============----------
	print '[OIII] vs radio power'
	GalaxiesFile = '%s/Data/galaxies_properties.txt' % (cc.base_dir)
	radio_power =  np.loadtxt(GalaxiesFile, unpack=True, skiprows=2, usecols=(1))
	galaxies =  np.loadtxt(GalaxiesFile, unpack=True, skiprows=2, usecols=(0,), dtype=str)

	fig, ax = plt.subplots()

	# MUSE
	m_gals = np.array([np.where(galaxies==g)[0][0] for g in galaxies_muse])
	ax.scatter(np.log10(OIII_eqw_muse), radio_power[m_gals], c='r', marker='x', 
		label='MUSE')

	# VIMOS
	v_gals = np.array([np.where(galaxies==g)[0][0] for g in galaxies_vimos])
	ax.scatter(np.log10(OIII_eqw_vimos), radio_power[v_gals], c='b', marker='x', 
		label='VIMOS')

	first = True
	for i_muse, g in enumerate(galaxies_muse):
		if g in galaxies_vimos:
			i_vimos = np.where(galaxies_vimos==g)[0][0]
			if first: # add just one label to legend
				ax.plot(np.log10([OIII_eqw_muse[i_muse], OIII_eqw_vimos[i_vimos]]), 
					[radio_power[m_gals[i_muse]],radio_power[v_gals[i_vimos]]], 'k--', 
					zorder=1, label='same galaxy in MUSE and VIMOS')
				first = False
			else:
				ax.plot(np.log10([OIII_eqw_muse[i_muse], OIII_eqw_vimos[i_vimos]]), 
					[radio_power[m_gals[i_muse]],radio_power[v_gals[i_vimos]]], 'k--', 
					zorder=1)

	# Atlas3D
	atlas3d_file = '%s/Data/atlas3d/XXXI_tableA1.dat' % (cc.base_dir)
	galaxies_atlas, radio_atlas = np.loadtxt(atlas3d_file, unpack=True, usecols=(0,10), 
		skiprows=2, dtype=str)
	m = np.array(['<' not in r for r in radio_atlas])
	atlas3d_file = '%s/Data/atlas3d/XXXI_tableA6.dat' % (cc.base_dir)
	galaxies_atlas2, OIII_eqw_atlas = np.loadtxt(atlas3d_file, unpack=True, 
		usecols=(0,8), skiprows=2, dtype=str)#(str,float), missing_values='-')
	m2 = np.array([r!='-' for r in OIII_eqw_atlas])

	a_gals = np.array([np.where(galaxies_atlas==g)[0][0] for g in galaxies_atlas2[m2] 
		if g in galaxies_atlas[m]])
	a_gals2 = np.array([np.where(galaxies_atlas2==g)[0][0] for g in galaxies_atlas2[m2]
		if g in galaxies_atlas[m]])

	ax.scatter(OIII_eqw_atlas[a_gals2].astype(float), radio_atlas[a_gals].astype(float), 
		marker='x', c='grey', label='Atlas3D')


	ax.axvline(np.log(0.8), color='k', linestyle=':', label='AGN limit')

	ax.legend(facecolor='w')
	ax.set_xlabel(r'log(EW [OIII]/$\mathrm{\AA}$)')
	ax.set_ylabel(r'$\log(P_{1.4 \mathrm{G Hz}} / \mathrm{W \, Hz^{-1}})$')

	fig.savefig('%s/Data/muse/analysis/OIIIew_radio.png' % (cc.base_dir))	
	plt.close()

	

if __name__=='__main__':
	compare_atlas3d()