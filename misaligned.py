## Routine to investigate misalignments between kinematic and photometry axes.

from checkcomp import checkcomp
cc = checkcomp()
if 'home' not in cc.device:
	import matplotlib # 20160202 JP to stop lack-of X-windows error
	matplotlib.use('Agg') 
import numpy as np 
import matplotlib.pyplot as plt
from prefig import Prefig
Prefig(transparent=False)

# Arguments and output are in degrees
def min_angle_diff(a1, a2, parity='odd'):
	if parity=='odd':
		return np.degrees(np.arcsin(np.abs(np.sin(np.radians(a1 - a2)))))
	elif parity=='even':
		return np.degrees(np.arccos(np.cos(np.radians(a1 - a2))))




def misaligned():
	print 'Misaligment plots'

## ----------====== Stellar misalignment vs Gas misalignment =======----------
	museGalaxiesFile = "%s/Data/muse/analysis/galaxies2.txt" % (cc.base_dir)
	vimosGalaxiesFile = "%s/Data/vimos/analysis/galaxies2.txt" % (cc.base_dir)

	PA_phot_muse, PA_kine_star_muse, PA_kine_gas_muse =  np.loadtxt(museGalaxiesFile, 
		unpack=True, skiprows=1, usecols=(3,4,5))
	star_mis_muse = min_angle_diff(PA_kine_star_muse, PA_phot_muse)
	gas_mis_muse = min_angle_diff(PA_kine_gas_muse, PA_phot_muse)
	galaxies_muse =  np.loadtxt(museGalaxiesFile, unpack=True, skiprows=1, usecols=(0,), 
		dtype=str)

	PA_phot_vimos, PA_kine_star_vimos, PA_kine_gas_vimos =  np.loadtxt(vimosGalaxiesFile, 
		unpack=True, skiprows=1, usecols=(3,4,5))
	star_mis_vimos = min_angle_diff(PA_kine_star_vimos, PA_phot_vimos)
	gas_mis_vimos = min_angle_diff(PA_kine_gas_vimos, PA_phot_vimos)
	galaxies_vimos =  np.loadtxt(vimosGalaxiesFile, unpack=True, skiprows=1, usecols=(0,), 
		dtype=str)

	fig, ax = plt.subplots()

	ax.scatter(star_mis_muse, gas_mis_muse, c='r', marker='x', label='MUSE')
	ax.scatter(star_mis_vimos, gas_mis_vimos, c='b', marker='x', label='VIMOS')

	first = True
	for i_gal_m, g in enumerate(galaxies_muse):
		if g in galaxies_vimos:
			i_gal_v = np.where(galaxies_vimos==g)[0][0]
			if first:
				ax.plot([star_mis_vimos[i_gal_v], star_mis_muse[i_gal_m]], 
					[gas_mis_vimos[i_gal_v], gas_mis_muse[i_gal_m]], 'k--',
					label='Same galaxy')
				first = False
			else:
				ax.plot([star_mis_vimos[i_gal_v], star_mis_muse[i_gal_m]], 
					[gas_mis_vimos[i_gal_v], gas_mis_muse[i_gal_m]], 'k--')

	xlim = ax.get_xlim()
	ylim = ax.get_ylim()
	ax.plot([-10,90], [-10,90], 'k:', label='Co/Counter rotating')
	ax.axvspan(-10,15, color='k', alpha=0.4, label='Misaligned limit')
	ax.axhspan(-10, 15, color='k', alpha=0.4)
	# ax.set_xlim(xlim)
	# ax.set_ylim(ylim)
	ax.set_xlim([0,90])
	ax.set_ylim([0,90])
	ax.legend(facecolor='w', fontsize=8)
	ax.set_xlabel('Stellar misalignment')
	ax.set_ylabel('Gas misalignment')
	ax.set_aspect('equal')

	# plt.show()
	plt.close()
## ----------====== Stellar misalignment vs Radio power =======----------
	print 'Radio power vs misalignment'
	GalaxiesFile = '%s/Data/galaxies_properties.txt' % (cc.base_dir)
	radio_power =  np.loadtxt(GalaxiesFile, unpack=True, skiprows=2, usecols=(1))
	galaxies =  np.loadtxt(GalaxiesFile, unpack=True, skiprows=2, usecols=(0,), dtype=str)

	atlas3d_file = '%s/Data/atlas3d/II_tableD1.dat' % (cc.base_dir)
	star_mis_altas = np.loadtxt(atlas3d_file, unpack=True, usecols=(7,))
	galaxies_altas = np.loadtxt(atlas3d_file, unpack=True, usecols=(0,), dtype=str)

	atlas3d_file = '%s/Data/atlas3d/XXXI_tableA1.dat' % (cc.base_dir)
	galaxies_altas2, radio_atlas = np.loadtxt(atlas3d_file, unpack=True, usecols=(0,10), 
		skiprows=2, dtype=str)

	fig, ax = plt.subplots()
	# VIMOS
	v_gals = np.array([np.where(galaxies==g)[0][0] for g in galaxies_vimos])
	ax.scatter(star_mis_vimos, radio_power[v_gals], c='b', marker='x', label='VIMOS')

	# MUSE
	m_gals = np.array([np.where(galaxies==g)[0][0] for g in galaxies_muse])
	ax.scatter(star_mis_muse, radio_power[m_gals], c='r', marker='x', label='MUSE')

	# Atlas3D
	m = np.array(['<' not in r for r in radio_atlas])
	radio_atlas_ob = radio_atlas[m].astype(float)
	a_gals = np.array([np.where(galaxies_altas==g)[0][0] for g in galaxies_altas2[m]])
	ax.scatter(star_mis_altas[a_gals], radio_atlas_ob, c='grey', marker='x', 
		label='Atlas3D')

	radio_atlas_lim = np.array([r.replace('<','') for r in radio_atlas[~m]])
	a_gals = np.array([np.where(galaxies_altas==g)[0][0] for g in galaxies_altas2[~m]])
	ax.quiver(star_mis_altas[a_gals], radio_atlas_lim, np.zeros(len(a_gals)), 
		np.ones(len(a_gals))*-0.03, scale=1, color='grey', width=0.002)


	first = True
	for i_gal_m, g in enumerate(galaxies_muse):
		if g in galaxies_vimos:
			i_gal_v = np.where(galaxies_vimos==g)[0][0]
			if first:
				ax.plot([star_mis_vimos[i_gal_v], star_mis_muse[i_gal_m]], 
					[radio_power[v_gals[i_gal_v]], radio_power[m_gals[i_gal_m]]], 'k--',
					label='Same galaxy')
				first = False
			else:
				ax.plot([star_mis_vimos[i_gal_v], star_mis_muse[i_gal_m]], 
					[radio_power[v_gals[i_gal_v]], radio_power[m_gals[i_gal_m]]], 'k--')
	ax.axvspan(0,15, color='k', alpha=0.4, label='Misaligned limit')
	ax.legend(facecolor='w', fontsize=8)
	ax.set_xlim([0,90])
	ax.set_xlabel(r'$\psi_\mathrm{kin - phot}$ (deg)')
	ax.set_ylabel(r'$\log(P_{1.4 \mathrm{G Hz}} / \mathrm{W \, Hz^{-1}})$')
	# plt.show()
	plt.close()

## ----------======== Gas misalignment vs Radio power =========----------
	fig, ax = plt.subplots()
	gas_mis_vimos = min_angle_diff(PA_kine_gas_vimos, PA_kine_star_vimos, parity='even')
	gas_mis_muse = min_angle_diff(PA_kine_gas_muse, PA_kine_star_muse, parity='even')


	fig, ax = plt.subplots()
	# VIMOS
	v_gals = np.array([np.where(galaxies==g)[0][0] for g in galaxies_vimos])
	ax.scatter(gas_mis_vimos, radio_power[v_gals], c='b', marker='x', label='VIMOS')

	# MUSE
	m_gals = np.array([np.where(galaxies==g)[0][0] for g in galaxies_muse])
	ax.scatter(gas_mis_muse, radio_power[m_gals], c='r', marker='x', label='MUSE')

	first = True
	for i_gal_m, g in enumerate(galaxies_muse):
		if g in galaxies_vimos:
			i_gal_v = np.where(galaxies_vimos==g)[0][0]
			if first:
				ax.plot([gas_mis_vimos[i_gal_v], gas_mis_muse[i_gal_m]], 
					[radio_power[v_gals[i_gal_v]], radio_power[m_gals[i_gal_m]]], 'k--',
					label='Same galaxy')
				first = False
			else:
				ax.plot([gas_mis_vimos[i_gal_v], gas_mis_muse[i_gal_m]], 
					[radio_power[v_gals[i_gal_v]], radio_power[m_gals[i_gal_m]]], 'k--')
	ax.axvspan(0,30, color='k', alpha=0.4, label='Misaligned limit')
	ax.legend(facecolor='w', fontsize=8)
	ax.set_xlim([0,180])
	ax.set_xlabel(r'$\psi_\mathrm{gas - stars}$ (deg)')
	ax.set_ylabel(r'$\log(P_{1.4 \mathrm{G Hz}} / \mathrm{W \, Hz^{-1}})$')
	# plt.show()
	plt.close()


	


if __name__=='__main__':
	misaligned()