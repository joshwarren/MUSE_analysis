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

def compare_atlas3d():
	print 'Compare to Atlas3d/SAURON'

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

if __name__=='__main__':
	compare_atlas3d()