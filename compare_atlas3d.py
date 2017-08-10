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
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.coordinates import match_coordinates_sky

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

class galaxy(object):
	def __init__(self, **kwargs):
		self.name = ''
		self.lambda_Re = np.nan
		self.ellipticity = np.nan
		self.T_type = np.nan
		self.structure = ''
		self.kin_group = ''
		
		for k, v in kwargs.iteritems():
			setattr(self, k, v)

	@property
	def FR(self):
		return bool((self.lambda_Re > 0.08 + self.ellipticity/4) + (
			self.ellipticity > 0.4))
	@property
	def E(self):
		return bool(self.T_type < -3.5)
	@property
	def a(self): # No Rotation
		return bool((self.kin_group == 'a') + (self.structure=='NRR/LV'))
	@property
	def b(self): # Complex Rotation
		return bool((self.kin_group == 'b') + (self.structure=='NRR/NF'))
	@property
	def c(self): # KDC
		return bool((self.kin_group == 'c') + ((self.structure=='NRR/KDC') + 
			(self.structure=='NRR/CRC') + (self.structure=='RR/CRC')))
	@property
	def d(self): # Counter-rotating disks
		return bool((self.kin_group == 'd') + ((self.structure=='NRR/2s') + 
			(self.structure=='RR/2s')))
	@property
	def e(self): # Regular Rotator
		return bool((self.kin_group == 'e') + ((self.structure=='RR/NF') + 
			(self.structure=='RR/2m') + (self.structure=='RR/KT')))
	@property
	def f(self): # Unclassified - have included with complex vel
		return bool(self.kin_group == 'f')


class galaxy_list(list):
	def __init__(self):
		pass

	def add_galaxy(self, g):
		if g.name not in self.names:
			self.append(g)
	def add_radio(self, galaxies, radio):
		for g in self:
			if g.name in galaxies:
				if '<' in radio[galaxies==g.name][0]:
					g.radio_limit = float(radio[galaxies==g.name][0].strip('<'))
					g.radio = np.nan
				else: 
					g.radio_limit = np.nan
					g.radio = float(radio[galaxies==g.name][0])
			else:
				g.radio = np.nan
				g.radio_limit = np.nan

	def add_first(self, first_coords, first_radio):
		first_coords = SkyCoord(first_coords, unit=(u.hourangle, u.deg))
		for g in self:
			idx, ang,_ = match_coordinates_sky(g.coords, first_coords)
			if ang.arcsec < 30:
				g.first_index = idx
				g.distance_to_first = ang.arcsec
				g.first_radio = np.log10((first_radio[idx]/1000*u.Jy * \
					4*np.pi*(g.distance*10**6*u.pc)**2).decompose().value) # in W/Hz
				g.first_coords = first_coords[idx]
			else:
				g.first_index = np.nan
				g.distance_to_first = np.nan
				g.first_radio = np.nan
				g.first_coords = np.nan
		for g in self:
			if np.sum([g2.first_index==g.first_index for g2 in self]) > 1:
				if g.distance_to_first != np.min([g2.distance_to_first for g2 
					in self if g2.first_index == g.first_index]):
					g.first_index = np.nan
					g.distance_to_first = np.nan
					g.first_radio = np.nan
					g.first_coords = np.nan
	def create_galaxy(self, name, **kwargs):
		if name not in self.names:
			kwargs['name'] = name
			self.append(galaxy(**kwargs))
		else:
			raise 'Galaxy named %s already in the galaxy_list object' % (name)

	def __getattr__(self, attr):
		return np.array([getattr(g, attr) for g in self])

	@property
	def names(self):
		return [g.name for g in self]



def compare_atlas3d():
	print 'Compare to Atlas3d/SAURON'

## ----------============== Ellipticity vs lambda_Re ==============----------
	print 'FR/SR'
	
	fig, ax = plt.subplots()

	# Plot Atlas3d results NB: all atlas3d tables are in alphabetical order
	atlas3d_file = '%s/Data/atlas3d/I_table3.dat' % (cc.base_dir)
	RA_atlas, dec_atlas, distance_atlas, M_k_atlas, T_type = np.loadtxt(
		atlas3d_file, unpack=True, usecols=(1,2,7,8,10))
	galaxies_atlas = np.loadtxt(atlas3d_file, unpack=True, usecols=(0), dtype=str)
	atlas3d_file = '%s/Data/atlas3d/III_tableB1.dat' % (cc.base_dir)
	ellipticity_atlas, lambda_Re_atlas = np.loadtxt(atlas3d_file, unpack=True, 
		usecols=(2,7), dtype=float)
	atlas3d_file = '%s/Data/atlas3d/II_tableD1.dat' % (cc.base_dir)
	structure_atlas = np.loadtxt(atlas3d_file, unpack=True, usecols=(13,), dtype=str)
	
	atlas_gals = galaxy_list()
	atlas_gals.color = 'k'
	atlas_gals.label = 'Atlas3D'
	for i, n in enumerate(galaxies_atlas):
		atlas_gals.create_galaxy(n, 
			M_k = M_k_atlas[i], 
			distance = distance_atlas[i],
			coords = SkyCoord(str(RA_atlas[i])+' '+str(dec_atlas[i]), unit='deg'),
			ellipticity = ellipticity_atlas[i],
			lambda_Re = lambda_Re_atlas[i],
			structure = structure_atlas[i])

	# Circle selected Atlas3D galaxies
	atlas_selected_file = '%s/Data/atlas3d/selected_galaxies.txt' % (cc.base_dir)
	select_atlas_galaxies = np.loadtxt(atlas_selected_file, unpack=True, usecols=(0,), 
		dtype=str, skiprows=1)
	VNESS_flux, VLSS_flux = np.loadtxt(atlas_selected_file, unpack=True, usecols=(1,3), 
		skiprows=1)
	for g in atlas_gals:
		if g.name in select_atlas_galaxies:
			i = np.where(select_atlas_galaxies == g.name)[0][0]
			g.VNESS_27 = VNESS_flux[i]
			g.VLSS_27 = VLSS_flux[i]
			g.selected_27 = True
		else:
			g.VLSS_27 = np.nan
			g.VNESS_27 = np.nan
			g.selected_27 = False
	ax.scatter(atlas_gals.ellipticity[atlas_gals.selected_27], 
		atlas_gals.lambda_Re[atlas_gals.selected_27], marker='o', edgecolor='r',
		facecolors='none', s=150, lw=1, label='Atlas3D subsample from 2.7GHz')

	# S0s
	ax.scatter(atlas_gals.ellipticity[atlas_gals.a*~atlas_gals.E], 
		atlas_gals.lambda_Re[atlas_gals.a*~atlas_gals.E], 
		marker=marker_atlas3d(0), c='lightgrey', alpha=0.5, lw=0,
		label=r'Atlas3D S0: $T>-3.5$')
	ax.scatter(atlas_gals.ellipticity[atlas_gals.b*~atlas_gals.E], 
		atlas_gals.lambda_Re[atlas_gals.b*~atlas_gals.E], 
		marker=marker_atlas3d(1), c='lightgrey', alpha=0.5, lw=0)
	ax.scatter(atlas_gals.ellipticity[atlas_gals.c*~atlas_gals.E], 
		atlas_gals.lambda_Re[atlas_gals.c*~atlas_gals.E], 
		marker=marker_atlas3d(2), c='lightgrey', alpha=0.5, lw=0)
	ax.scatter(atlas_gals.ellipticity[atlas_gals.d*~atlas_gals.E], 
		atlas_gals.lambda_Re[atlas_gals.d*~atlas_gals.E], 
		marker=marker_atlas3d(3), c='lightgrey', alpha=0.5, lw=0)
	ax.plot(atlas_gals.ellipticity[atlas_gals.e*~atlas_gals.E], 
		atlas_gals.lambda_Re[atlas_gals.e*~atlas_gals.E], 
		marker=marker_atlas3d(4), c='lightgrey', alpha=0.5, lw=0, 
		markerfacecolor='none')

	# Ellipticals
	ax.scatter(atlas_gals.ellipticity[atlas_gals.a*atlas_gals.E], 
		atlas_gals.lambda_Re[atlas_gals.a*atlas_gals.E], 
		marker=marker_atlas3d(0), c='k', alpha=0.5, lw=0, 
		label=r'Atlas3D E: $T \leq -3.5$')
	ax.scatter(atlas_gals.ellipticity[atlas_gals.b*atlas_gals.E], 
		atlas_gals.lambda_Re[atlas_gals.b*atlas_gals.E], 
		marker=marker_atlas3d(1), c='k', alpha=0.5, lw=0)
	ax.scatter(atlas_gals.ellipticity[atlas_gals.c*atlas_gals.E], 
		atlas_gals.lambda_Re[atlas_gals.c*atlas_gals.E], 
		marker=marker_atlas3d(2), c='k', alpha=0.5, lw=0)
	ax.scatter(atlas_gals.ellipticity[atlas_gals.d*atlas_gals.E], 
		atlas_gals.lambda_Re[atlas_gals.d*atlas_gals.E], 
		marker=marker_atlas3d(3), c='k', alpha=0.5, lw=0)
	ax.plot(atlas_gals.ellipticity[atlas_gals.e*atlas_gals.E], 
		atlas_gals.lambda_Re[atlas_gals.e*atlas_gals.E], 
		marker=marker_atlas3d(4), c='k', alpha=0.5, lw=0, markerfacecolor='none')
	

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
	museGalaxiesFile = "%s/Data/muse/analysis/galaxies2.txt" % (cc.base_dir)
	muse_classify_file = "%s/Data/muse/analysis/galaxies_classify.txt" % (cc.base_dir)

	lambda_Re_muse, ellipticity_muse =  np.loadtxt(museGalaxiesFile, unpack=True, 
		skiprows=1, usecols=(1,2))
	galaxies_muse =  np.loadtxt(museGalaxiesFile, unpack=True, skiprows=1, 
		usecols=(0,), dtype=str)
	

	gals_muse2, group_muse = np.loadtxt(muse_classify_file, unpack=True, 
		usecols=(0,8), dtype=str, skiprows=1)

	muse_gals = galaxy_list()
	muse_gals.color = 'b'
	muse_gals.label = 'MUSE'
	for i in range(len(galaxies_muse)):
		muse_gals.create_galaxy(galaxies_muse[i], 
			lambda_Re=lambda_Re_muse[i],
			ellipticity=ellipticity_muse[i])
		i2 = np.where(gals_muse2 == galaxies_muse[i])[0][0]
		muse_gals[i].kin_group = group_muse[i2]

	# ax.scatter(muse_gals.ellipticity[muse_gals.a], muse_gals.lambda_Re[muse_gals.a], 
	# 	marker=marker_atlas3d(0), c='b', lw=0, label='MUSE')
	# ax.scatter(muse_gals.ellipticity[muse_gals.b], muse_gals.lambda_Re[muse_gals.b], 
	# 	marker=marker_atlas3d(1), c='b', lw=0)
	# ax.scatter(muse_gals.ellipticity[muse_gals.c], muse_gals.lambda_Re[muse_gals.c], 
	# 	marker=marker_atlas3d(2), c='b', lw=0)
	# ax.scatter(muse_gals.ellipticity[muse_gals.d], muse_gals.lambda_Re[muse_gals.d], 
	# 	marker=marker_atlas3d(3), c='b', lw=0)
	# ax.plot(muse_gals.ellipticity[muse_gals.e], muse_gals.lambda_Re[muse_gals.e], 
	# 	marker=marker_atlas3d(4), c='b', lw=0, markerfacecolor='none')


	# VIMOS
	vimosGalaxiesFile = "%s/Data/vimos/analysis/galaxies2.txt" % (cc.base_dir)
	vimos_classify_file = "%s/Data/vimos/analysis/galaxies_classify_by_eye.txt" % (
		cc.base_dir)

	lambda_Re_vimos, ellipticity_vimos =  np.loadtxt(vimosGalaxiesFile, unpack=True, 
		skiprows=1, usecols=(1,2))
	galaxies_vimos =  np.loadtxt(vimosGalaxiesFile, unpack=True, skiprows=1, 
		usecols=(0,), dtype=str)
	
	gals_vimos2, group_vimos = np.loadtxt(vimos_classify_file, unpack=True, 
		usecols=(0,8), dtype=str, skiprows=1)

	vimos_gals = galaxy_list()
	vimos_gals.color = 'r'
	vimos_gals.label = 'VIMOS'
	for i in range(len(galaxies_vimos)):
		vimos_gals.create_galaxy(galaxies_vimos[i], lambda_Re=lambda_Re_vimos[i],
			ellipticity=ellipticity_vimos[i])
		i2 = np.where(gals_vimos2 == galaxies_vimos[i])[0][0]
		vimos_gals[i].kin_group = group_vimos[i2]

	for gals in [vimos_gals, muse_gals]:
		ax.scatter(gals.ellipticity[gals.a], gals.lambda_Re[gals.a], 
			marker=marker_atlas3d(0), c=gals.color, lw=0, label=gals.label)
		ax.scatter(gals.ellipticity[gals.b], gals.lambda_Re[gals.b], 
			marker=marker_atlas3d(1), c=gals.color, lw=0)
		ax.scatter(gals.ellipticity[gals.f], gals.lambda_Re[gals.f], 
			marker=marker_atlas3d(1), c=gals.color, lw=0)
		ax.scatter(gals.ellipticity[gals.c], gals.lambda_Re[gals.c], 
			marker=marker_atlas3d(2), c=gals.color, lw=0)
		ax.scatter(gals.ellipticity[gals.d], gals.lambda_Re[gals.d], 
			marker=marker_atlas3d(3), c=gals.color, lw=0)
		ax.plot(gals.ellipticity[gals.e], gals.lambda_Re[gals.e], 
			marker=marker_atlas3d(4), c=gals.color, lw=0, markerfacecolor='none')

	# Join MUSE and VIMOS
	label=True
	for i_muse, g in enumerate(muse_gals):
		if g.name in vimos_gals.names:
			i_vimos = np.where(vimos_gals.name==g.name)[0][0]
			if label: # add just one label to legend
				ax.plot([muse_gals.ellipticity[i_muse], 
					vimos_gals.ellipticity[i_vimos]], 
					[muse_gals.lambda_Re[i_muse],vimos_gals.lambda_Re[i_vimos]], 
					'k--', zorder=1, label='same galaxy in MUSE and VIMOS')
				label = False
			else:
				ax.plot([muse_gals.ellipticity[i_muse], 
					vimos_gals.ellipticity[i_vimos]], 
					[muse_gals.lambda_Re[i_muse],vimos_gals.lambda_Re[i_vimos]], 
					'k--', zorder=1)

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
	expectedSRs_err = 0

	for i, ell in enumerate(np.arange(0,0.9,0.1)):
		ell_bin_atlas = (atlas_gals.ellipticity >= ell ) * (
			atlas_gals.ellipticity < ell+0.1)
		SRfraction_atlas.append(np.sum(~atlas_gals.FR*ell_bin_atlas)/
			float(np.sum(ell_bin_atlas)))
		ell_uncert_bin_atlas = np.sqrt(SRfraction_atlas[i] * (
			1 - SRfraction_atlas[i])/np.sum(ell_bin_atlas))

		ell_bin_vimos = (vimos_gals.ellipticity >= ell ) * (
			vimos_gals.ellipticity < ell+0.1)
		SRfraction_vimos.append(np.sum(~vimos_gals.FR*ell_bin_vimos)/
			float(np.sum(ell_bin_vimos)))

		ell_bin_muse = (muse_gals.ellipticity >= ell ) * (
			muse_gals.ellipticity < ell+0.1)
		SRfraction_muse.append(np.sum(~muse_gals.FR*ell_bin_muse)/
			float(np.sum(ell_bin_muse)))

		expectedSRs = np.nansum((np.sum(ell_bin_vimos) * SRfraction_atlas[i], 
			expectedSRs))
		expectedSRs_err = np.sqrt(expectedSRs_err**2 + (
			ell_uncert_bin_atlas * np.sum(ell_bin_vimos))**2)


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
		"Expected # of SRs in our sample \n   based on Atlas3D: (%.2f+/-%.2f)/10" % (
			expectedSRs, expectedSRs_err))
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

	GalaxiesFile = '%s/Data/galaxies_properties.txt' % (cc.base_dir)
	radio, M_k =  np.loadtxt(GalaxiesFile, unpack=True, skiprows=2, usecols=(1,2))
	galaxies =  np.loadtxt(GalaxiesFile, unpack=True, skiprows=2, usecols=(0,), dtype=str)

	for g in vimos_gals:
		if g.name in galaxies:
			i_gal = np.where(galaxies==g.name)[0][0]
			g.M_k = M_k[i_gal]
			g.radio = radio[i_gal]
		else:
			g.M_k = np.nan
			g.radio = np.nan

	for g in muse_gals:
		if g.name in galaxies:
			i_gal = np.where(galaxies==g.name)[0][0]
			g.M_k = M_k[i_gal]
			g.radio = radio[i_gal]
		else:
			g.M_k = np.nan
			g.radio = np.nan

	for gals in [vimos_gals, muse_gals, atlas_gals]:
		ax[2].scatter(gals.M_k[gals.FR], gals.lambda_Re[gals.FR], c=gals.color, 
			marker='x', label='%s Fast Rotators' % (gals.label))
		ax[2].scatter(gals.M_k[~gals.FR], gals.lambda_Re[~gals.FR], c=gals.color, 
			marker='^', label='%s Slow Rotators' % (gals.label))
	ax[2].scatter(atlas_gals.M_k[atlas_gals.selected_27], 
		atlas_gals.lambda_Re[atlas_gals.selected_27], marker='o', edgecolor='r',
		facecolors='none', s=150, lw=1, label='Atlas3D subsample from 2.7GHz')


	# Join MUSE and VIMOS
	label = True
	for i_muse, g in enumerate(muse_gals):
		if g.name in vimos_gals.names:
			i_vimos = np.where(vimos_gals.name==g.name)[0][0]
			if label: # add just one label to legend
				ax[2].plot([muse_gals.M_k[i_muse], vimos_gals.M_k[i_vimos]], 
					[muse_gals.lambda_Re[i_muse],vimos_gals.lambda_Re[i_vimos]], 
					'k--', zorder=1, label='same galaxy in MUSE and VIMOS')
				label = False
			else:
				ax[2].plot([muse_gals.M_k[i_muse], vimos_gals.M_k[i_vimos]], 
					[muse_gals.lambda_Re[i_muse],vimos_gals.lambda_Re[i_vimos]], 
					'k--', zorder=1)

	plt.legend(facecolor='w')

	_,_,h1= ax[1].twinx().hist(atlas_gals.M_k[atlas_gals.FR], histtype='step', color='k', 
		normed=True, label='Atlas3d Fast Rotators')
	plt.yticks([])
	_,_,h2 = ax[1].twinx().hist(atlas_gals.M_k[~atlas_gals.FR], histtype='step', color='k', 
		normed=True, linestyle='--', label='Atlas3d Slow Rotators')
	plt.yticks([])
	h = [h1[0],h2[0]]
	for gals in [vimos_gals, muse_gals]:
		label = True
		for l in gals.M_k[gals.FR]:
			if label:
				h.append(ax[1].axvline(l, color=gals.color, 
					label='%s Fast Rotators' % (gals.label)))
				label = False
			else:
				ax[1].axvline(l, color=gals.color)
		label = True
		for l in gals.M_k[~gals.FR]:
			if label:
				h.append(ax[1].axvline(l, color=gals.color, ls='--',
					label='%s Slow Rotators' % (gals.label)))
				label = False
			else:
				ax[1].axvline(l, color=gals.color, ls='--')
	plt.legend(handles=h, facecolor='w', loc='center left')


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
	expectedSRs_err = 0

	# step = 1#0.5
	# steps = np.arange(-27,-21,step)
	# step = [2,1,1,1,1]#,1]
	# steps = [-27,-25,-24,-23,-22]#,-21]

	order = np.argsort(atlas_gals.M_k)
	steps = [-27.]+list(atlas_gals.M_k[order][~atlas_gals.FR[order]][::-5][::-1]+0.1)
	step = np.diff(steps+[-21.])

	for i, M in enumerate(steps):
		# ax[2].axvline(M) # Show bins in main plot
		M_k_bin_atlas = (atlas_gals.M_k >= M ) * (atlas_gals.M_k < M + step[i])
		SRfraction_atlas.append(np.sum(~atlas_gals.FR * M_k_bin_atlas)/
			float(np.sum(M_k_bin_atlas)))
		M_k_uncert_bin_atlas = np.sqrt(SRfraction_atlas[i] * (1 - SRfraction_atlas[i])/
			np.sum(M_k_bin_atlas))


		M_k_bin_vimos = (vimos_gals.M_k >= M ) * (vimos_gals.M_k < M + step[i])
		SRfraction_vimos.append(np.sum(~vimos_gals.FR * M_k_bin_vimos)/
			float(np.sum(M_k_bin_vimos)))

		M_k_bin_muse = (muse_gals.M_k >= M ) * (muse_gals.M_k < M + step[i])
		SRfraction_muse.append(np.sum(~muse_gals.FR * M_k_bin_muse)/
			float(np.sum(M_k_bin_muse)))
		expectedSRs = np.nansum((np.sum(M_k_bin_vimos) * SRfraction_atlas[i], 
			expectedSRs))
		expectedSRs_err = np.sqrt(expectedSRs_err**2 + (
			M_k_uncert_bin_atlas * np.sum(M_k_bin_vimos))**2)

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
	ax[0].plot(steps, SRfraction_atlas, color='k', ls='steps--')
	ax[0].plot(steps, SRfraction_vimos, color='r', ls='steps--')
	ax[0].plot(steps, SRfraction_muse, color='b', ls='steps--')
	ax[0].set_ylim([-0.05,1.05])
	ax[0].text(-21,0.75, 
		"Expected # of SRs in our sample \n   based on Atlas3D: (%.2f+/-%.2f)/10" % (
			expectedSRs, expectedSRs_err))

	fig.suptitle('K-band magnitude distribution for F/S rotators')
	fig.savefig('%s/Data/muse/analysis/lambda_R_M_k.png' % (cc.base_dir))
	plt.close()
	Prefig(transparent=False)
## ----------============== ellipticity vs M_k ===============----------
	print 'Ellipticity vs M_k'

	fig, ax = plt.subplots()
	for gals in [vimos_gals, atlas_gals]:
		ax.scatter(gals.M_k[gals.FR], gals.ellipticity[gals.FR], color=gals.color, 
			marker='x', label='%s Fast Rotators' % (gals.label))
		ax.scatter(gals.M_k[~gals.FR], gals.ellipticity[~gals.FR], color=gals.color, 
			marker='^', label='%s Slow Rotators' % (gals.label))
	ax.scatter(atlas_gals.M_k[atlas_gals.selected_27], 
		atlas_gals.ellipticity[atlas_gals.selected_27], marker='o', edgecolor='r',
		facecolors='none', s=150, lw=1, label='Atlas3D subsample from 2.7GHz')


	M_k_steps = [-27,-25,-24,-23,-22,-21]
	ell_steps = np.arange(0,0.9,0.1)
	expectedSRs = 0
	expectedSRs_err = 0
	for i in range(len(M_k_steps)-1):
		ax.axvline(M_k_steps[i], c='k', alpha=0.5)
		for j in range(len(ell_steps)-1):
			for gals in [vimos_gals, atlas_gals]:
				M = (gals.M_k >= M_k_steps[i]) * (gals.M_k < M_k_steps[i+1])
				E = (gals.ellipticity >= ell_steps[j]) * (
					gals.ellipticity < ell_steps[j+1])
				if gals.label == 'VIMOS':
					n_vimos_in_bin = np.sum(M*E)
				elif gals.label == 'Atlas3D':
					SRfraction_atlas_bin = np.sum(~gals.FR[M*E])/float(np.sum(M*E))

					expectedSRs = np.nansum([expectedSRs, 
						n_vimos_in_bin*SRfraction_atlas_bin])
					expectedSRs_err = np.sqrt(np.nansum(
						[expectedSRs_err**2, ((SRfraction_atlas_bin * (
						1 - SRfraction_atlas_bin)/np.sum(M*E)) * n_vimos_in_bin)**2]))

	for e in ell_steps:
		ax.axhline(e, c='k', alpha=0.5)

	ax.text(-21.5, 0.8, 
		"Expected # of SRs in our sample \n   based on Atlas3D: (%.2f+/-%.2f)/10" % (
			expectedSRs, expectedSRs_err))

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
	atlas_gals.add_radio(galaxies_atlas2, radio_atlas)

	fig, ax = plt.subplots()

	for gals in [vimos_gals, atlas_gals]:
		ax.scatter(gals.M_k[gals.FR], gals.radio[gals.FR], color=gals.color, marker='x',
			label='%s Fast Rotators' % (gals.label))
		ax.scatter(gals.M_k[~gals.FR], gals.radio[~gals.FR], color=gals.color, marker='^',
			label = '%s Slow Rotators' % (gals.label))
	ax.scatter(atlas_gals.M_k[atlas_gals.selected_27], 
		atlas_gals.radio[atlas_gals.selected_27], marker='o', edgecolor='r',
		facecolors='none', s=150, lw=1, label='Atlas3D subsample from 2.7GHz')

	ax.invert_xaxis()
	plt.legend(facecolor='w')

	ax.set_title('K-band magnitude vs Radio Power')
	ax.set_xlabel(r'$M_k \mathrm{(mag)}$')
	ax.set_ylabel(r'$\log P_{1.4\mathrm{GHz}}$')

	fig.savefig('%s/Data/muse/analysis/radio_power_M_k.png' % (cc.base_dir))
	plt.close('all')

## ----------=========== Radio power (FIRST) vs M_k ===============----------
	print 'Radio Power (with FIRST) vs M_K'
	fig, ax = plt.subplots()

	for gals in [vimos_gals, atlas_gals]:
		ax.scatter(gals.M_k[gals.FR], gals.radio[gals.FR], color=gals.color, marker='x',
			label='%s Fast Rotators' % (gals.label))
		ax.scatter(gals.M_k[~gals.FR], gals.radio[~gals.FR], color=gals.color, marker='^',
			label = '%s Slow Rotators' % (gals.label))
	ax.scatter(atlas_gals.M_k[atlas_gals.selected_27], 
		atlas_gals.radio[atlas_gals.selected_27], marker='o', edgecolor='r',
		facecolors='none', s=150, lw=1, label='Atlas3D subsample from 2.7GHz')
	

	# Add the FIRST radio power to Atlas galaxies
	first_file = '%s/Data/atlas3d/atlas_in_first.txt' % (cc.base_dir)
	RA_h, RA_m, RA_s, dec_d, dec_m, dec_s, radio_first = np.loadtxt(first_file,
		unpack=True, usecols=(5,6,7,8,9,10,13), skiprows=4, 
		dtype=str)#'int,int,float,int,int,float,float')
	radio_first = radio_first.astype(float)

	coords_first = []
	for i in range(len(RA_h)):
		coords_first.append(RA_h[i]+' '+RA_m[i]+' '+RA_s[i]+' '+dec_d[i]+' '+
			dec_m[i]+' '+dec_s[i])
	atlas_gals.add_first(coords_first, radio_first)

	ax.scatter(atlas_gals.M_k[atlas_gals.FR], atlas_gals.first_radio[atlas_gals.FR],
		color='orange', marker='x', label='FIRST Fast Rotators')
	ax.scatter(atlas_gals.M_k[~atlas_gals.FR], atlas_gals.first_radio[~atlas_gals.FR], 
		color='orange', marker='^', label='FIRST Slow Roatators')
	ax.scatter(atlas_gals.M_k[atlas_gals.selected_27], 
		atlas_gals.first_radio[atlas_gals.selected_27], marker='o', edgecolor='r',
		facecolors='none', s=150, lw=1)

	for i, g in enumerate(atlas_gals):
		ax.plot([g.M_k, g.M_k], [g.first_radio, g.radio], c='k', ls='--')
		# if ~np.isnan(g.first_radio):
		# 	ax.text(g.M_k, g.first_radio, g.name, )

	ax.invert_xaxis()
	plt.legend(facecolor='w')

	ax.set_title('K-band magnitude vs FIRST Radio Power')
	ax.set_xlabel(r'$M_k \mathrm{(mag)}$')
	ax.set_ylabel(r'$\log P_{1.4\mathrm{GHz}}$')

	fig.savefig('%s/Data/muse/analysis/first_radio_power_M_k.png' % (cc.base_dir))
	plt.close()
## ----------================ Core age vs KDC size ================----------
	print 'KDC size/age'
	muse_core_file = "%s/Data/muse/analysis/galaxies_core.txt" % (cc.base_dir)
	vimos_core_file = "%s/Data/vimos/analysis/galaxies_core.txt" % (cc.base_dir)
	sauron_file = '%s/Data/sauron/VIII_table8.dat' % (cc.base_dir)
	
	fig, ax = plt.subplots()
	# Sauron
	size_sauron, size_unc_sauron, age_sauron, age_unc_sauron = np.loadtxt(sauron_file,
		unpack=True, usecols=(4,5,6,7), skiprows=2) 
	galaxies_sauron, fast_sauron = np.loadtxt(sauron_file, unpack=True, usecols=(0,8), 
		skiprows=2, dtype=str)

	sauron_gals = galaxy_list()
	for i, n in enumerate(galaxies_sauron):
		sauron_gals.create_galaxy('NGC'+n,
			KDC_size = size_sauron[i],
			KDC_size_uncert = size_unc_sauron[i],
			age = age_sauron[i],
			age_uncert = age_unc_sauron[i])
	sauron_gals.FR = fast_sauron=='F'

	ax.errorbar(sauron_gals.KDC_size[sauron_gals.FR], sauron_gals.age[sauron_gals.FR], 
		fmt='.', xerr=sauron_gals.KDC_size_uncert[sauron_gals.FR], 
		yerr=sauron_gals.age_uncert[sauron_gals.FR], color='k',
		label='Fast rotating SAURON')
	ax.errorbar(sauron_gals.KDC_size[~sauron_gals.FR], sauron_gals.age[~sauron_gals.FR], 
		fmt='.', xerr=sauron_gals.KDC_size_uncert[~sauron_gals.FR], 
		yerr=sauron_gals.age_uncert[~sauron_gals.FR],
		color='lightgrey', label='Slow rotating SAURON')

	# MUSE
	age_muse, age_unc_muse, OIII_eqw_muse, OIII_eqw_unc_muse= np.loadtxt(muse_core_file, 
		unpack=True, usecols=(1,2,7,8), skiprows=2)
	gals_muse1 = np.loadtxt(muse_core_file, unpack=True, usecols=(0,), 
		skiprows=2, dtype=str)
	gals_muse2, size_muse = np.loadtxt(muse_classify_file, unpack=True, usecols=(0,5), 
		dtype=str, skiprows=1)
	for g in muse_gals:
		if g.name in gals_muse1:
			i = np.where(gals_muse1 == g.name)[0][0]
			g.age = age_muse[i]
			g.age_uncert = age_unc_muse[i]
			g.OIIIew = OIII_eqw_muse[i]
			g.OIIIew_uncert = OIII_eqw_unc_muse[i]
		else:
			g.age = np.nan
			g.age_uncert = np.nan
			g.OIIIew = np.nan
			g.OIIIew_uncert = np.nan

		if g.name in gals_muse2:
			i = np.where(gals_muse2 == g.name)[0][0]
			if size_muse[i] !='-':
				g.KDC_size = angle_to_pc(g.name, float(size_muse[i]))
				g.KDC_size_uncert = angle_to_pc(g.name, 0.5)
			else: 
				g.KDC_size = np.nan
				g.KDC_size_uncert = np.nan
		else: 
			g.KDC_size = np.nan
			g.KDC_size_uncert = np.nan

	# VIMOS
	age_vimos, age_unc_vimos, OIII_eqw_vimos, OIII_eqw_unc_vimos = np.loadtxt(
		vimos_core_file, unpack=True, usecols=(1,2,7,8), skiprows=2)
	gals_vimos1 = np.loadtxt(vimos_core_file, unpack=True, usecols=(0,), 
		skiprows=2, dtype=str)
	gals_vimos2, size_vimos = np.loadtxt(vimos_classify_file, unpack=True, usecols=(0,5), 
		dtype=str, skiprows=1)
	for g in vimos_gals:
		if g.name in gals_vimos1:
			i = np.where(gals_vimos1 == g.name)[0][0]
			g.age = age_vimos[i]
			g.age_uncert = age_unc_vimos[i]
			g.OIIIew = OIII_eqw_vimos[i]
			g.OIIIew_uncert = OIII_eqw_unc_vimos[i]
		else:
			g.age = np.nan
			g.age_uncert = np.nan
			g.OIIIew = np.nan
			g.OIIIew_uncert = np.nan

		if g.name in gals_vimos2:
			i = np.where(gals_vimos2 == g.name)[0][0]
			if size_vimos[i] !='-':
				g.KDC_size = angle_to_pc(g.name, float(size_vimos[i]))
				g.KDC_size_uncert = angle_to_pc(g.name, 0.5)
			else: 
				g.KDC_size = np.nan
				g.KDC_size_uncert = np.nan
		else: 
			g.KDC_size = np.nan
			g.KDC_size_uncert = np.nan
	
	for gals in [vimos_gals, muse_gals]:
		ax.errorbar(gals.KDC_size, gals.age, fmt='.', xerr=gals.KDC_size_uncert,
			yerr=gals.age_uncert, color=gals.color, label=gals.label)
	
	# Join MUSE and VIMOS
	label = True
	for i_muse, g in enumerate(muse_gals):
		if g.name in vimos_gals.names:
			i_vimos = np.where(vimos_gals.name==g.name)[0][0]
			if label: # add just one label to legend
				ax.plot([muse_gals.KDC_size[i_muse], vimos_gals.KDC_size[i_vimos]], 
					[muse_gals.age[i_muse],vimos_gals.age[i_vimos]], 
					'k--', zorder=1, label='same galaxy in MUSE and VIMOS')
				label = False
			else:
				ax.plot([muse_gals.KDC_size[i_muse], vimos_gals.KDC_size[i_vimos]], 
					[muse_gals.age[i_muse],vimos_gals.age[i_vimos]], 
					'k--', zorder=1)



	ax.legend(facecolor='w')
	ax.set_yscale('log')#, nonposy='clip', subsy=[1,2,3,4,5,6,7,8,9])
	ax.set_xlabel('KDC size (pc)')
	ax.set_ylabel('Age of central 1 arcsec (Gyrs)')
	ax.set_title('Age and size of KDCs')


	fig.savefig('%s/Data/muse/analysis/KDC_size_age.png' % (cc.base_dir))	
	plt.close()
## ----------=========== Core OIII vs radio power ============----------
	print '[OIII] vs radio power'
	fig, ax = plt.subplots()

	for gals in [vimos_gals, muse_gals]:
		ax.errorbar(np.log10(gals.OIIIew), gals.radio, c=gals.color, fmt='x',
			label=gals.label, xerr=abs(np.log(gals.OIIIew - gals.OIIIew_uncert) - 
			np.log(gals.OIIIew)))

	first = True
	for i_muse, g in enumerate(galaxies_muse):
		if g in galaxies_vimos:
			i_vimos = np.where(galaxies_vimos==g)[0][0]
			if first: # add just one label to legend
				ax.plot(np.log10([muse_gals.OIIIew[i_muse], vimos_gals.OIIIew[i_vimos]]), 
					[muse_gals.radio[i_muse],vimos_gals.radio[i_vimos]], 'k--', 
					zorder=1, label='same galaxy in MUSE and VIMOS')
				first = False
			else:
				ax.plot(np.log10([muse_gals.OIIIew[i_muse], vimos_gals.OIIIew[i_vimos]]), 
					[muse_gals.radio[i_muse],vimos_gals.radio[i_vimos]], 'k--', zorder=1)

	
	atlas3d_file = '%s/Data/atlas3d/XXXI_tableA6.dat' % (cc.base_dir)
	galaxies_atlas2, OIII_eqw_atlas = np.loadtxt(atlas3d_file, unpack=True, 
		usecols=(0,8), skiprows=2, dtype=str)#(str,float), missing_values='-')
	for g in atlas_gals:
		if g.name in galaxies_atlas2:
			i = np.where(galaxies_atlas2 == g.name)[0][0]
			if OIII_eqw_atlas[i] != '-':
				g.OIIIew = float(OIII_eqw_atlas[i])
			else: g.OIIIew = np.nan
		else: g.OIIIew = np.nan

	ax.scatter(atlas_gals.OIIIew, atlas_gals.radio, marker='x', c='k', label='Atlas3D')
	ax.scatter(atlas_gals.OIIIew[atlas_gals.selected_27], 
		atlas_gals.radio[atlas_gals.selected_27], marker='o', edgecolor='r',
		facecolors='none', s=150, lw=1, label='Atlas3D subsample from 2.7GHz')

	ax.axvline(np.log(0.8), color='k', linestyle=':', label='AGN limit')

	ax.legend(facecolor='w')
	ax.set_xlabel(r'log(EW [OIII]/$\mathrm{\AA}$)')
	ax.set_ylabel(r'$\log(P_{1.4 \mathrm{G Hz}} / \mathrm{W \, Hz^{-1}})$')

	fig.savefig('%s/Data/muse/analysis/OIIIew_radio.png' % (cc.base_dir))	
	plt.close()

## ----------========== Radio spectral index ===========----------
	Prefig(size=(16,12*2), transparent=False)

	galaxy_properties2_file = '%s/Data/galaxies_properties2.txt' % (cc.base_dir)
	q24, q24_err, spec_index, spec_index_err = np.loadtxt(galaxy_properties2_file,
		unpack=True, usecols=(1,2,3,4), skiprows=1)
	galaxies = np.loadtxt(galaxy_properties2_file, usecols=(0,), skiprows=1, dtype=str)

	for g in vimos_gals:
		if g.name in galaxies:
			i = np.where(galaxies == g.name)[0][0]
			g.q24 = q24[i]
			g.q24_err = q24_err[i]
			g.spec_index = spec_index[i]
			g.spec_index_err = spec_index_err[i]

	for g in muse_gals:
		if g.name in galaxies:
			i = np.where(galaxies == g.name)[0][0]
			g.q24 = q24[i]
			g.q24_err = q24_err[i]
			g.spec_index = spec_index[i]
			g.spec_index_err = spec_index_err[i]



	fig,ax =plt.subplots(2, sharex=True)
	ax[0].scatter(vimos_gals.q24[vimos_gals.FR], vimos_gals.radio[vimos_gals.FR], 
		marker='x', c='k')
	ax[1].scatter(vimos_gals.q24[vimos_gals.FR], vimos_gals.spec_index[vimos_gals.FR], 
		marker='x', c='k', label='Fast Rotators')

	ax[0].scatter(vimos_gals.q24[~vimos_gals.FR], vimos_gals.radio[~vimos_gals.FR], 
		marker='^', c='k')
	ax[1].scatter(vimos_gals.q24[~vimos_gals.FR], vimos_gals.spec_index[~vimos_gals.FR], 
		marker='^', c='k', label = 'Slow Rotators')

	# Use MUSE too to include ngc1316
	ax[0].scatter(muse_gals.q24[muse_gals.FR], muse_gals.radio[muse_gals.FR], 
		marker='x', c='k')
	ax[1].scatter(muse_gals.q24[muse_gals.FR], muse_gals.spec_index[muse_gals.FR], 
		marker='x', c='k')

	ax[0].scatter(muse_gals.q24[~muse_gals.FR], muse_gals.radio[~muse_gals.FR], 
		marker='^', c='k')
	ax[1].scatter(muse_gals.q24[~muse_gals.FR], muse_gals.spec_index[~muse_gals.FR], 
		marker='^', c='k')

	for a in ax:
		a.axvspan(0.8,1.2, alpha=0.4)
	ax[0].axhline(24, ls=':', c='k')
	ax[1].axhline(-0.5, ls=':', c='k')
	lim = ax[0].get_xlim()
	ax[0].text(lim[0]*0.97, 24.05, 'High-powered')
	ax[0].text(lim[0]*0.97, 23.85, 'Low-powered')

	ax[1].text(lim[0]*0.97, -0.47, 'Flat')
	ax[1].text(lim[0]*0.97, -0.65, 'Steep')

	ax[0].set_ylabel(r'log(P$_\mathrm{1.4GHz}$)')
	ax[1].set_ylabel('Radio spectral index')
	ax[1].set_xlabel(r'q$_{24}$')
	ax[0].set_title('Radio properties')

	ax[1].legend(facecolor='w', loc=4)
	ax[0].tick_params(top=True, bottom=True, direction='in')
	ax[1].tick_params(top=True, direction='in')
	fig.subplots_adjust(hspace=0)
	fig.savefig('%s/Data/muse/analysis/radio_spectral_index.png' % (cc.base_dir))

	plt.close('all')
	Prefig(size=(16,12), transparent=False)

if __name__=='__main__':
	compare_atlas3d()