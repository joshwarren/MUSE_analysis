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
	vimos_classify_file = "%s/Data/vimos/analysis/galaxies_classify_by_eye.txt" % (
		cc.base_dir)


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
	M_k_atlas, T_type = np.loadtxt(atlas3d_file, unpack=True, usecols=(8,10))
	galaxies_atlas = np.loadtxt(atlas3d_file, unpack=True, usecols=(0), dtype=str)
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

	# Isotropic line: add envolope - from footnote, page 430, SAURON X
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
	gals_muse2, group_muse = np.loadtxt(muse_classify_file, unpack=True, usecols=(0,8), 
		dtype=str, skiprows=1)
	gal_order = [np.where(galaxies_muse==g)[0][0] for g in gals_muse2]
	a_muse = (group_muse=='a')[gal_order] # No rotation
	b_muse = (group_muse=='b')[gal_order] # Complex rotation
	c_muse = (group_muse=='c')[gal_order] # KDC
	d_muse = (group_muse=='d')[gal_order] # Counter-rotating disks
	e_muse = (group_muse=='e')[gal_order] # Regular rotator
	f_muse = (group_muse=='f')[gal_order] # Unclassified

	ax.scatter(ellipticity_muse[a_muse], lambda_Re_muse[a_muse], 
		marker=marker_atlas3d(0), c='b', lw=0, label='MUSE')
	ax.scatter(ellipticity_muse[b_muse], lambda_Re_muse[b_muse], 
		marker=marker_atlas3d(1), c='b', lw=0)
	ax.scatter(ellipticity_muse[c_muse], lambda_Re_muse[c_muse], 
		marker=marker_atlas3d(2), c='b', lw=0)
	ax.scatter(ellipticity_muse[d_muse], lambda_Re_muse[d_muse], 
		marker=marker_atlas3d(3), c='b', lw=0)
	ax.plot(ellipticity_muse[e_muse], lambda_Re_muse[e_muse], 
		marker=marker_atlas3d(4), c='b', lw=0, markerfacecolor='none')


	# VIMOS
	gals_vimos2, group_vimos = np.loadtxt(vimos_classify_file, unpack=True, usecols=(0,8), 
		dtype=str, skiprows=1)
	print group_vimos=='f'
	gal_order = [np.where(galaxies_vimos==g)[0][0] for g in gals_vimos2]
	a_vimos = (group_vimos=='a')[gal_order] # No rotation
	b_vimos = (group_vimos=='b')[gal_order] # Complex rotation
	c_vimos = (group_vimos=='c')[gal_order] # KDC
	d_vimos = (group_vimos=='d')[gal_order] # Counter-rotating disks
	e_vimos = (group_vimos=='e')[gal_order] # Regular rotator
	f_vimos = (group_vimos=='f')[gal_order] # Unclassified - have included with complex vel

	ax.scatter(ellipticity_vimos[a_vimos], lambda_Re_vimos[a_vimos], 
		marker=marker_atlas3d(0), c='r', lw=0, label='VIMOS')
	ax.scatter(ellipticity_vimos[b_vimos], lambda_Re_vimos[b_vimos], 
		marker=marker_atlas3d(1), c='r', lw=0)
	ax.scatter(ellipticity_vimos[f_vimos], lambda_Re_vimos[f_vimos], 
		marker=marker_atlas3d(1), c='r', lw=0)
	ax.scatter(ellipticity_vimos[c_vimos], lambda_Re_vimos[c_vimos], 
		marker=marker_atlas3d(2), c='r', lw=0)
	ax.scatter(ellipticity_vimos[d_vimos], lambda_Re_vimos[d_vimos], 
		marker=marker_atlas3d(3), c='r', lw=0)
	ax.plot(ellipticity_vimos[e_vimos], lambda_Re_vimos[e_vimos], 
		marker=marker_atlas3d(4), c='r', lw=0, markerfacecolor='none')

	ax.set_title('Atlas3D Fast/Slow Rotator Classification scheme')
	ax.set_xlabel(r'$\epsilon$')
	ax.set_ylabel(r'$\lambda_R (R_e)$')
	
	# Plot Slow Rotator bounds
	ax.plot([0,0.4,0.4], [0.08, 0.18, 0], 'k', label='FR/SR boundary')
	leg = plt.legend(facecolor='w')
	ax.add_artist(leg)

	# Proxy for kinematics legend
	h1, = ax.plot(np.nan, marker=marker_atlas3d(0), c='b', lw=0, 
		label='No Rotation')
	h2, = ax.plot(np.nan, marker=marker_atlas3d(1), c='b', lw=0,
		label='Complex Rotation')
	h3, = ax.plot(np.nan, marker=marker_atlas3d(2), c='b', lw=0, 
		label='KDC/CDC')
	h4, = ax.plot(np.nan, marker=marker_atlas3d(3), c='b', lw=0,
		label='Counter Rotating Disks')
	h5, = ax.plot(np.nan, marker=marker_atlas3d(4), c='b', lw=0,
		label='Regular Rotator',  markerfacecolor='none')

	plt.legend(handles=[h1,h2,h3,h4,h5], facecolor='w', loc=5)

	# Show fraction of Slow Rotators in background of plot per elliptiity bin with 
	# width 0.1
	SRfraction_atlas = []
	SRfraction_vimos = []
	SRfraction_muse = []
	expectedSRs = 0

	FR_atlas = (lambda_Re_atlas > 0.08 + ellipticity_atlas/4) + (ellipticity_atlas > 0.4)
	FR_vimos = (lambda_Re_vimos > 0.08 + ellipticity_vimos/4) + (ellipticity_vimos > 0.4)
	FR_muse = (lambda_Re_muse > 0.08 + ellipticity_muse/4) + (ellipticity_muse > 0.4)
	for i, ell in enumerate(np.arange(0,0.9,0.1)):
		ell_bin_atlas = (ellipticity_atlas >= ell ) * (ellipticity_atlas < ell+0.1)
		SRfraction_atlas.append(np.sum(~FR_atlas*ell_bin_atlas)/
			float(np.sum(ell_bin_atlas)))

		ell_bin_vimos = (ellipticity_vimos >= ell ) * (ellipticity_vimos < ell+0.1)
		SRfraction_vimos.append(np.sum(~FR_vimos*ell_bin_vimos)/
			float(np.sum(ell_bin_vimos)))

		ell_bin_muse = (ellipticity_muse >= ell ) * (ellipticity_muse < ell+0.1)
		SRfraction_muse.append(np.sum(~FR_muse*ell_bin_muse)/
			float(np.sum(ell_bin_muse)))

		expectedSRs = np.nansum((np.sum(ell_bin_vimos) * SRfraction_atlas[i], expectedSRs))

	SRfraction_atlas.insert(0,SRfraction_atlas[0])
	SRfraction_vimos.insert(0,SRfraction_vimos[0])
	SRfraction_muse.insert(0,SRfraction_muse[0])

	SRfraction_atlas = np.array(SRfraction_atlas)
	SRfraction_vimos = np.array(SRfraction_vimos)
	SRfraction_muse = np.array(SRfraction_muse)

	SRfraction_atlas[np.isnan(SRfraction_atlas)] = 0
	SRfraction_vimos[np.isnan(SRfraction_vimos)] = 0
	SRfraction_muse[np.isnan(SRfraction_muse)] = 0

	ax2 = ax.twinx()
	ax2.set_ylabel('Fraction of Slow Rotators in ellipticity bin', rotation=270, 
		labelpad=25)
	ax2.plot(np.arange(0,1,0.1), SRfraction_atlas, color='k', alpha=0.3, ls='steps')
	ax2.plot(np.arange(0,1,0.1), SRfraction_vimos, color='r', alpha=0.3, ls='steps')
	ax2.plot(np.arange(0,1,0.1), SRfraction_muse, color='b', alpha=0.3, ls='steps')
	ax2.text(0.02,0.9, 
		'Expected # of SRs in our sample \n   based on Atlas3D: %.2f/10' % (expectedSRs))

	ax2.set_ylim([0,1.05])
	ax.set_xlim([0, 0.9])
	ax.set_ylim([0, 0.8])

	# Save plot
	fig.savefig('%s/Data/muse/analysis/lambda_R_ellipticity.png' % (cc.base_dir))
	plt.close()	
## ----------============ K-band magnitude vs lambda_R ===========----------
	print 'K-band magnitude vs lambda_Re'
	Prefig(size=(16,12*1.7), transparent=False)
	fig, ax = plt.subplots(3,1, sharex=True,  gridspec_kw = {'height_ratios':[1, 1, 3]})
	ax[2].scatter(M_k_atlas[FR_atlas], lambda_Re_atlas[FR_atlas], marker='x', c='k', 
		label='Atlas3D Fast Rotators')
	ax[2].scatter(M_k_atlas[~FR_atlas], lambda_Re_atlas[~FR_atlas], marker='^', c='k', 
		label='Atlas3D Slow Rotators')

	GalaxiesFile = '%s/Data/galaxies_properties.txt' % (cc.base_dir)
	radio, M_k =  np.loadtxt(GalaxiesFile, unpack=True, skiprows=2, usecols=(1,2))
	galaxies =  np.loadtxt(GalaxiesFile, unpack=True, skiprows=2, usecols=(0,), dtype=str)

	v_gals = []
	for g in galaxies_vimos:
		i_gal = np.where(galaxies==g)[0][0]
		v_gals.append(i_gal)

	m_gals = []
	for g in galaxies_muse:
		i_gal = np.where(galaxies==g)[0][0]
		m_gals.append(i_gal)

	ax[2].scatter(M_k[v_gals][FR_vimos], lambda_Re_vimos[FR_vimos], c='r', marker='x',
		label='VIMOS Fast Rotators')
	ax[2].scatter(M_k[v_gals][~FR_vimos], lambda_Re_vimos[~FR_vimos], c='r', marker='^',
		label='VIMOS Slow Rotators')
	ax[2].scatter(M_k[m_gals][FR_muse], lambda_Re_muse[FR_muse], c='b', marker='x',
		label='MUSE Fast Rotators')
	ax[2].scatter(M_k[m_gals][~FR_muse], lambda_Re_muse[~FR_muse], c='b', marker='^',
		label='MUSE Slow Rotators')

	# Join MUSE and VIMOS
	first = True
	for i_muse, g in enumerate(galaxies_muse):
		if g in galaxies_vimos:
			i_vimos = np.where(galaxies_vimos==g)[0][0]
			if first: # add just one label to legend
				ax[2].plot([M_k[m_gals][i_muse],M_k[v_gals][i_vimos]], 
					[lambda_Re_muse[i_muse], lambda_Re_vimos[i_vimos]], 'k--', 
					zorder=1, label='same galaxy in MUSE and VIMOS')
				first = False
			else:
				ax[2].plot([M_k[m_gals][i_muse],M_k[v_gals][i_vimos]], 
					[lambda_Re_muse[i_muse], lambda_Re_vimos[i_vimos]], 'k--', 
					zorder=1)



	plt.legend(facecolor='w')

	_,_,h1= ax[1].twinx().hist(M_k_atlas[FR_atlas], histtype='step', color='k', 
		normed=True, label='Atlas3d Fast Rotators')
	plt.yticks([])
	_,_,h2 = ax[1].twinx().hist(M_k_atlas[~FR_atlas], histtype='step', color='k', 
		normed=True, linestyle='--', label='Atlas3d Slow Rotators')
	plt.yticks([])
	# _,_,h3 = ax[1].twinx().hist(M_k[v_gals][FR_vimos], histtype='step', color='r',
	# 	label='VIMOS Fast Rotators')
	# plt.yticks([])
	# _,_,h4 = ax[1].twinx().hist(M_k[v_gals][~FR_vimos], histtype='step', color='r', 
	# 	linestyle='--', label='VIMOS Slow Rotators')
	# plt.yticks([])
	# _,_,h5 = ax[1].twinx().hist(M_k[m_gals][FR_muse], histtype='step', color='b',
	# 	label='MUSE Fast Rotators')
	# plt.yticks([])
	# _,_,h6 = ax[1].twinx().hist(M_k[m_gals][~FR_muse], histtype='step', color='b',
	# 	linestyle='--', label='MUSE Slow Rotators')
	# plt.yticks([])
	for l in M_k[v_gals][FR_vimos]:
		h3 = ax[1].axvline(l, color='r', label='VIMOS Fast Rotators')
	for l in M_k[v_gals][~FR_vimos]:
		h4 = ax[1].axvline(l, color='r', ls='--', label='VIMOS Slow Rotators')
	for l in M_k[m_gals][FR_muse]:
		h5 = ax[1].axvline(l, color='b', label='MUSE Fast Rotators')
	for l in M_k[m_gals][~FR_muse]:
		h6 = ax[1].axvline(l, color='b', label='MUSE Slow Rotators')

	plt.legend(handles=[h1[0],h2[0],h3,h4,h5,h6], facecolor='w', loc='center left')

	ax[2].set_xlabel(r'$M_k \mathrm{(mag)}$')
	ax[2].set_ylabel(r'$\lambda_R (R_e)$')
	ax[1].set_ylabel('Normalised histogram')
	ax[0].set_ylabel('Fraction of SR')
	ax[2].invert_xaxis()


	# Show fraction of Slow Rotators in background of plot per M_k bin with width 0.5
	SRfraction_atlas = []
	SRfraction_vimos = []
	SRfraction_muse = []
	expectedSRs = 0

	# step = 1#0.5
	# steps = np.arange(-27,-21,step)
	# step = [2,1,1,1,1]#,1]
	# steps = [-27,-25,-24,-23,-22]#,-21]

	order = np.argsort(M_k_atlas)
	steps = [-27.]+list(M_k_atlas[order][~FR_atlas[order]][::-5][::-1]+0.1)
	step = np.diff(steps+[-21.])

	for i, M in enumerate(steps):
		# ax[2].axvline(M) # Show bins in main plot
		M_k_bin_atlas = (M_k_atlas >= M ) * (M_k_atlas < M + step[i])
		SRfraction_atlas.append(np.sum(~FR_atlas * M_k_bin_atlas)/
			float(np.sum(M_k_bin_atlas)))

		M_k_bin_vimos = (M_k[v_gals] >= M ) * (M_k[v_gals] < M + step[i])
		SRfraction_vimos.append(np.sum(~FR_vimos * M_k_bin_vimos)/
			float(np.sum(M_k_bin_vimos)))

		M_k_bin_muse = (M_k[m_gals] >= M ) * (M_k[m_gals] < M + step[i])
		SRfraction_muse.append(np.sum(~FR_muse * M_k_bin_muse)/
			float(np.sum(M_k_bin_muse)))
		expectedSRs = np.nansum((np.sum(M_k_bin_vimos) * SRfraction_atlas[i], expectedSRs))

	SRfraction_atlas.insert(0,SRfraction_atlas[0])
	SRfraction_vimos.insert(0,SRfraction_vimos[0])
	SRfraction_muse.insert(0,SRfraction_muse[0])

	SRfraction_atlas = np.array(SRfraction_atlas)
	SRfraction_vimos = np.array(SRfraction_vimos)
	SRfraction_muse = np.array(SRfraction_muse)

	SRfraction_atlas[np.isnan(SRfraction_atlas)] = 0
	SRfraction_vimos[np.isnan(SRfraction_vimos)] = 0
	SRfraction_muse[np.isnan(SRfraction_muse)] = 0

	steps.append(-21)
	# ax[0].plot(np.arange(-27,-21+step,step), SRfraction_atlas, color='k', ls='steps--')
	# ax[0].plot(np.arange(-27,-21+step,step), SRfraction_vimos, color='r', ls='steps--')
	# ax[0].plot(np.arange(-27,-21+step,step), SRfraction_muse, color='b', ls='steps--')
	ax[0].plot(steps, SRfraction_atlas, color='k', ls='steps--')
	ax[0].plot(steps, SRfraction_vimos, color='r', ls='steps--')
	ax[0].plot(steps, SRfraction_muse, color='b', ls='steps--')
	ax[0].set_ylim([-0.05,1.05])
	ax[0].text(-21,0.75, 
		'Expected # of SRs in our sample \n   based on Atlas3D: %.2f/10' % (expectedSRs))

	fig.suptitle('K-band magnitude distribution for F/S rotators')
	fig.savefig('%s/Data/muse/analysis/lambda_R_M_k.png' % (cc.base_dir))
	plt.close()

	Prefig(transparent=False)
## ----------============== ellipticity vs M_k ===============----------
	print 'Ellipticity vs M_k'

	fig, ax = plt.subplots()
	ax.scatter(M_k_atlas[FR_atlas], ellipticity_atlas[FR_atlas], color='k', marker='x',
		label='Atlas3d Fast Rotators')
	ax.scatter(M_k_atlas[~FR_atlas], ellipticity_atlas[~FR_atlas], color='k', marker='^',
		label='Atlas3d Slow Rotators')

	ax.scatter(M_k[v_gals][FR_vimos], ellipticity_vimos[FR_vimos], color='r', marker='x',
		label='VIMOS Fast Rotators')
	ax.scatter(M_k[v_gals][~FR_vimos], ellipticity_vimos[~FR_vimos], color='r', marker='^',
		label='VIMOS Slow Rotators')

	M_k_steps = [-27,-25,-24,-23,-22,-21]
	ell_steps = np.arange(0,0.9,0.1)
	expectedSRs = 0
	for i in range(len(M_k_steps)-1):
		ax.axvline(M_k_steps[i], c='k', alpha=0.5)
		for j in range(len(ell_steps)-1):
			M = (M_k[v_gals] >= M_k_steps[i]) * (M_k[v_gals] < M_k_steps[i+1])
			E = (ellipticity_vimos >= ell_steps[i]) * (ellipticity_vimos < ell_steps[i+1])

			n_vimos_in_bin = np.sum(M*E)

			M = (M_k_atlas >= M_k_steps[i]) * (M_k_atlas < M_k_steps[i+1])
			E = (ellipticity_atlas >= ell_steps[i]) * (ellipticity_atlas < ell_steps[i+1])

			expectedSRs = np.nansum([expectedSRs, n_vimos_in_bin*np.sum(~FR_atlas[M*E])/
				float(np.sum(M*E))])
	for e in ell_steps:
		ax.axhline(e, c='k', alpha=0.5)

	ax.text(-21.5, 0.8, 
		'Expected # of SRs in our sample \n   based on Atlas3D: %.2f/10' % (expectedSRs))

	ax.set_title('K-band magnitude to elliticity relationship')
	ax.set_xlabel(r'$M_k \mathrm{(mag)}$')
	ax.set_ylabel(r'$\epsilon$')
	ax.invert_xaxis()

	plt.legend(facecolor='w')
	fig.savefig('%s/Data/muse/analysis/ellipticity_M_k.png' % (cc.base_dir))

## ----------============== Radio power vs M_k ===============----------
	print 'M_k vs Radio power'

	atlas3d_file = '%s/Data/atlas3d/XXXI_tableA1.dat' % (cc.base_dir)
	galaxies_atlas2, radio_atlas = np.loadtxt(atlas3d_file, unpack=True, usecols=(0,10), 
		skiprows=2, dtype=str)
	m = np.array(['<' not in r for r in radio_atlas])
	m2 = np.array([g in galaxies_atlas2 for g in galaxies_atlas])

	fig, ax = plt.subplots()

	ax.scatter(M_k[v_gals][FR_vimos], radio[v_gals][FR_vimos], color='b', 
		label='VIMOS Fast Rotators')
	ax.scatter(M_k[v_gals][~FR_vimos], radio[v_gals][~FR_vimos], color='c',
		label = 'VIMOS Slow Rotators')

	ax.scatter(M_k_atlas[m2][m*FR_atlas[m2]], radio_atlas[m*FR_atlas[m2]].astype(float), 
		color='r', label='Atlas3d Fast Rotators')
	ax.scatter(M_k_atlas[m2][m*~FR_atlas[m2]], radio_atlas[m*~FR_atlas[m2]].astype(float), 
		color='orange', label='Atlas3d Slow Rotators')
	ax.invert_xaxis()
	plt.legend(facecolor='w')

	ax.set_title('K-band magnitude vs Radio Power')
	ax.set_xlabel(r'$M_k \mathrm{(mag)}$')
	ax.set_ylabel(r'$\log P_{1.4\mathrm{GHz}}$')

	fig.savefig('%s/Data/muse/analysis/radio_power_M_k.png' % (cc.base_dir))
	plt.close()
## ----------================ Core age vs KDC size ================----------
	print 'KDC size/age'
	muse_core_file = "%s/Data/muse/analysis/galaxies_core.txt" % (cc.base_dir)
	vimos_core_file = "%s/Data/vimos/analysis/galaxies_core.txt" % (cc.base_dir)
	# vimosGalaxiesFile = "%s/Data/vimos/analysis/galaxies_core.txt" % (cc.base_dir)
	sauron_file = '%s/Data/sauron/VIII_table8.dat' % (cc.base_dir)
	
	fig, ax = plt.subplots()
	# Sauron
	size_sauron, size_unc_sauron, age_sauron, age_unc_sauron = np.loadtxt(sauron_file,
		unpack=True, usecols=(4,5,6,7), skiprows=2) 
	fast_sauron = np.loadtxt(sauron_file, unpack=True, usecols=(8,), skiprows=2, dtype=str)
	fast_sauron = fast_sauron == 'F'

	ax.errorbar(size_sauron[fast_sauron], age_sauron[fast_sauron], fmt='.',
		xerr=size_unc_sauron[fast_sauron], yerr=age_unc_sauron[fast_sauron], color='k',
		label='Fast rotating SAURON')
	ax.errorbar(size_sauron[~fast_sauron], age_sauron[~fast_sauron], fmt='.',
		xerr=size_unc_sauron[~fast_sauron], yerr=age_unc_sauron[~fast_sauron], 
		color='lightgrey', label='Slow rotating SAURON')

	# MUSE
	age_muse, age_unc_muse, OIII_eqw_muse = np.loadtxt(muse_core_file, unpack=True, 
		usecols=(1,2,7), skiprows=2)
	gals_muse1 = np.loadtxt(muse_core_file, unpack=True, usecols=(0,), 
		skiprows=2, dtype=str)
	gals_muse2, size_muse = np.loadtxt(muse_classify_file, unpack=True, usecols=(0,5), 
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
				yerr=age_unc_muse[i1], color='b', label='MUSE')
		else:
			ax.errorbar(size_muse[i2], age_muse[i1], fmt='.',xerr=size_unc_muse[i2], 
				yerr=age_unc_muse[i1], color='b')

	# VIMOS
	age_vimos, age_unc_vimos, OIII_eqw_vimos = np.loadtxt(vimos_core_file, unpack=True, 
		usecols=(1,2,7), skiprows=2)
	gals_vimos1 = np.loadtxt(vimos_core_file, unpack=True, usecols=(0,), 
		skiprows=2, dtype=str)
	gals_vimos2, size_vimos = np.loadtxt(vimos_classify_file, unpack=True, usecols=(0,5), 
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
				yerr=age_unc_vimos[i1], color='r', label='VIMOS')
		else:
			ax.errorbar(size_vimos[i2], age_vimos[i1], fmt='.',xerr=size_unc_vimos[i2], 
				yerr=age_unc_vimos[i1], color='r')
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