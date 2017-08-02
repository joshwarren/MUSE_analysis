# Routine to select Atlas galaxies based on the spectral index derived from WENSS, 
# VLSS and FIRST surveys. 

from checkcomp import checkcomp
cc = checkcomp()
# if 'home' not in cc.device:
# 	import matplotlib # 20160202 JP to stop lack-of X-windows error
# 	matplotlib.use('Agg') 
import numpy as np 
# import matplotlib.pyplot as plt
# from prefig import Prefig
# Prefig(transparent=False)
# from astropy import units as u
from astropy.io import fits
from astropy.coordinates import SkyCoord
from astropy.coordinates import match_coordinates_sky
from compare_atlas3d import galaxy_list

f = fits.open('%s/Data/atlas3d/FIRST_NVSS.fits' %(cc.base_dir))
d = f[1].data

atlas3d_gals = galaxy_list()
atlas3d_file = '%s/Data/atlas3d/I_table3.dat' % (cc.base_dir)
RA_atlas, dec_atlas, distance_atlas, M_k_atlas = np.loadtxt(atlas3d_file, 
	unpack=True, usecols=(1,2,7,8))
galaxies_atlas = np.loadtxt(atlas3d_file, unpack=True, usecols=(0), dtype=str)

atlas_gals = galaxy_list()
for i, n in enumerate(galaxies_atlas):
	atlas_gals.create_galaxy(n, M_k=M_k_atlas[i], distance=distance_atlas[i],
		coords=SkyCoord(str(RA_atlas[i])+' '+str(dec_atlas[i]), unit='deg'))

# WENSS
m = (d.FIRST_FINT != -99) * (d.WENSS_FLUX != -99) # Use mask to speed up reading 
												  #	to SkyCoords.
m *= np.abs(d.DEC-29) < 35 # Sky region for Atlas3d
FIRST_NVSS_coords = SkyCoord(zip(d.RA[m], d.DEC[m]), unit='deg')
for g in atlas_gals:
	idx, ang, _ = match_coordinates_sky(g.coords, FIRST_NVSS_coords)
	if ang.arcsec < 30:
		g.first = d.FIRST_FINT[m][idx]
		g.wenss = d.WENSS_FLUX[m][idx]
		g.spectral_index = (np.log(g.first)-np.log(g.wenss))/(
			np.log(1.4*10**9)-np.log(325*10**6)) # NB: Used positive sign convention
		g.wenss_projected_2_7 = g.first * (2.7/1.4)**g.spectral_index
		g.wenss_selected = g.wenss_projected_2_7 > 0.25
	else:
		g.first = np.nan
		f.wenss = np.nan
		g.spectral_index = np.nan
		g.wenss_projected_2_7 = np.nan
		g.wenss_selected = False


# VLSS
m = (d.FIRST_FINT != -99) * (d.VLSS_FLUX != -99) 
m *= np.abs(d.DEC-29) < 35 # Sky region for Atlas3d
FIRST_NVSS_coords = SkyCoord(zip(d.RA[m], d.DEC[m]), unit='deg')
for g in atlas_gals:
	idx, ang, _ = match_coordinates_sky(g.coords, FIRST_NVSS_coords)
	if ang.arcsec < 30:
		g.first = d.FIRST_FINT[m][idx]
		g.vlss = d.VLSS_FLUX[m][idx]
		g.spectral_index = (np.log(g.first)-np.log(g.vlss))/(
			np.log(1.4*10**9)-np.log(74*10**6)) # NB: Used positive sign convention

		g.vlss_projected_2_7 = g.first * (2.7/1.4)**g.spectral_index
		g.vlss_selected = g.vlss_projected_2_7 > 0.25
	else:
		g.first = np.nan
		f.wenss = np.nan
		g.spectral_index = np.nan
		g.vlss_projected_2_7 = np.nan
		g.vlss_selected = False


selected = atlas_gals.vlss_selected + atlas_gals.wenss_selected
with open('%s/Data/atlas3d/selected_galaxies.txt' % (cc.base_dir), 'w') as file:
	file.write('Galaxy   WENSS_flux   WENSS_flag  VLSS_flux  VLSS_flag\n')
	for i, g in enumerate(np.array(atlas_gals)[selected]):
		file.write(g.name+'  '+str(g.wenss_projected_2_7)+'  '+
			str(int(g.wenss_selected))+'  '+str(g.vlss_projected_2_7)+'  '+
			str(int(g.vlss_selected))+'\n')
