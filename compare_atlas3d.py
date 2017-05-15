## Routine to compare to Atlas3d results.

from checkcomp import checkcomp
cc = checkcomp()
if 'home' not in cc.device:
	import matplotlib # 20160202 JP to stop lack-of X-windows error
	matplotlib.use('Agg') 
import numpy as np 
import matplotlib.pyplot as plt 


def compare_atlas3d():
	print 'Compare to Atlas3d/SAURON'

	museGalaxiesFile = "%s/Data/muse/analysis/galaxies2.txt" % (cc.base_dir)
	vimosGalaxiesFile = "%s/Data/vimos/analysis/galaxies2.txt" % (cc.base_dir)

	lambda_Re_muse, ellipticity_muse =  np.loadtxt(museGalaxiesFile, unpack=True, 
		skiprows=1, usecols=(1,2))

	lambda_Re_vimos, ellipticity_vimos =  np.loadtxt(vimosGalaxiesFile, unpack=True, 
		skiprows=1, usecols=(1,2))


	fig, ax = plt.subplots()

	# Plot Atlas3d results NB: all atlas3d tables are in alphabetical order
	atlas3d_file = '%s/Data/atlas3d/III_tableB1.dat' % (cc.base_dir)
	ellipticity_atlas, lambda_Re_altas = np.loadtxt(atlas3d_file, unpack=True, 
		usecols=(2,7), dtype=float)
	atlas3d_file = '%s/Data/atlas3d/II_tableD1.dat' % (cc.base_dir)
	structure_altas = np.loadtxt(atlas3d_file, unpack=True, usecols=(13,), dtype=str)
	
	no_rot_altas = structure_altas=='NRR/LV'
	complex_rot_atlas = structure_altas=='NRR/NF'
	KDC_atlas = (structure_altas=='NRR/KDC') + (structure_altas=='NRR/CRC') + \
		(structure_altas=='RR/CRC')
	counter_rot_atlas = (structure_altas=='NRR/2s') + (structure_altas=='RR/2s')
	regular_rot_atlas = (structure_altas=='RR/NF') + (structure_altas=='RR/2m') + \
		(structure_altas=='RR/KT')

	ax.scatter(ellipticity_atlas[no_rot_altas], lambda_Re_altas[no_rot_altas], marker='o',
		c='grey', alpha=0.5)
	ax.scatter(ellipticity_atlas[complex_rot_atlas], lambda_Re_altas[complex_rot_atlas], 
		c='grey', alpha=0.5)
	ax.scatter(ellipticity_atlas[KDC_atlas], lambda_Re_altas[KDC_atlas], marker='^', 
		c='grey', alpha=0.5)
	ax.scatter(ellipticity_atlas[counter_rot_atlas], lambda_Re_altas[counter_rot_atlas], 
		c='grey', alpha=0.5)
	ax.scatter(ellipticity_atlas[regular_rot_atlas], lambda_Re_altas[regular_rot_atlas], 
		c='grey', alpha=0.5)
	
	# Plot scatter
	ax.scatter(ellipticity_muse, lambda_Re_muse, c='r')
	ax.scatter(ellipticity_vimos, lambda_Re_vimos, c='g')
	ax.set_xlabel(r'$\epsilon$')
	ax.set_ylabel(r'$\lambda_R (R_e)$')
	ax.set_xlim([0, 0.9])
	ax.set_ylim([0, 0.8])

	# Plot Slow Rotator bounds
	ax.plot([0,0.4,0.4], [0.08, 0.18, 0], 'k')

	# Save plot
	fig.savefig('%s/Data/muse/analysis/lambda_R_ellipticity.png' % (cc.base_dir))

if __name__=='__main__':
	compare_atlas3d()