## Routine to compare to Atlas3d results.

from checkcomp import checkcomp
cc = checkcomp()
if 'home' not in cc.device:
	import matplotlib # 20160202 JP to stop lack-of X-windows error
	matplotlib.use('Agg') 
import numpy as np 
import matplotlib.pyplot as plt
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
	ellipticity_atlas, lambda_Re_altas = np.loadtxt(atlas3d_file, unpack=True, 
		usecols=(2,7), dtype=float)
	atlas3d_file = '%s/Data/atlas3d/II_tableD1.dat' % (cc.base_dir)
	structure_altas = np.loadtxt(atlas3d_file, unpack=True, usecols=(13,), dtype=str)
	atlas3d_file = '%s/Data/atlas3d/I_table3.dat' % (cc.base_dir)
	T_type = np.loadtxt(atlas3d_file, unpack=True, usecols=(10,))
	E = T_type < -3.5
	
	no_rot_altas = structure_altas=='NRR/LV'
	complex_rot_atlas = structure_altas=='NRR/NF'
	KDC_atlas = (structure_altas=='NRR/KDC') + (structure_altas=='NRR/CRC') + \
		(structure_altas=='RR/CRC')
	counter_rot_atlas = (structure_altas=='NRR/2s') + (structure_altas=='RR/2s')
	regular_rot_atlas = (structure_altas=='RR/NF') + (structure_altas=='RR/2m') + \
		(structure_altas=='RR/KT')

	# S0s
	ax.scatter(ellipticity_atlas[no_rot_altas*~E], lambda_Re_altas[no_rot_altas*~E], 
		marker='o', c='lightgrey', alpha=0.5, lw=0)
	ax.scatter(ellipticity_atlas[complex_rot_atlas*~E], 
		lambda_Re_altas[complex_rot_atlas*~E], c='lightgrey', alpha=0.5, lw=0)
	ax.scatter(ellipticity_atlas[KDC_atlas*~E], lambda_Re_altas[KDC_atlas*~E], 
		marker='^', c='lightgrey', alpha=0.5, lw=0)
	ax.scatter(ellipticity_atlas[counter_rot_atlas*~E], 
		lambda_Re_altas[counter_rot_atlas*~E], c='lightgrey', alpha=0.5, lw=0)
	ax.scatter(ellipticity_atlas[regular_rot_atlas*~E], 
		lambda_Re_altas[regular_rot_atlas*~E], c='lightgrey', alpha=0.5, lw=0, 
		label=r'Atlas3D S0: $T>-3.5$')

	# Ellipticals
	ax.scatter(ellipticity_atlas[no_rot_altas*E], lambda_Re_altas[no_rot_altas*E], 
		marker='o', c='k', alpha=0.5, lw=0)
	ax.scatter(ellipticity_atlas[complex_rot_atlas*E], 
		lambda_Re_altas[complex_rot_atlas*E], c='k', alpha=0.5, lw=0)
	ax.scatter(ellipticity_atlas[KDC_atlas*E], lambda_Re_altas[KDC_atlas*E], marker='^', 
		c='k', alpha=0.5, lw=0)
	ax.scatter(ellipticity_atlas[counter_rot_atlas*E], 
		lambda_Re_altas[counter_rot_atlas*E], c='k', alpha=0.5, lw=0)
	ax.scatter(ellipticity_atlas[regular_rot_atlas*E], 
		lambda_Re_altas[regular_rot_atlas*E], c='k', alpha=0.5, lw=0, 
		label=r'Atlas3D E: $T \leq -3.5$')
	
	# Plot scatter
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

	ax.scatter(ellipticity_muse, lambda_Re_muse, c='b', lw=0, zorder=2, label='MUSE')
	ax.scatter(ellipticity_vimos, lambda_Re_vimos, c='r', lw=0, zorder=2, label='VIMOS')
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
	ax.scatter(np.log10(OIII_eqw_vimos), radio_power[v_gals], c='r', marker='x', 
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
				
	ax.legend(facecolor='w')
	ax.set_xlabel(r'log(EW [OIII]/$\mathrm{\AA}$)')
	ax.set_ylabel(r'$\log(P_{1.4 \mathrm{G Hz}} / \mathrm{W \, Hz^{-1}})$')

	fig.savefig('%s/Data/muse/analysis/OIIIew_radio.png' % (cc.base_dir))	
	plt.close()

	

if __name__=='__main__':
	compare_atlas3d()