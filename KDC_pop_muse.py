## ==================================================================
## 		Stellar population within inner 1 arcsec
## ==================================================================
# The entire pipeline - making use of error2_muse and pop_muse routines.


from checkcomp import checkcomp
cc = checkcomp()
import numpy as np 
from astropy.io import fits
from errors2_muse import get_dataCubeDirectory, run_ppxf, set_params, apply_range
from pop_muse import population
from KDC_pop import get_areaInAperture
from Bin import trapz_uncert
from prefig import Prefig 
Prefig()

def get_specFromAperture(galaxy, app_size=1.0, inside=True, res=0.2):
	f = fits.open(get_dataCubeDirectory(galaxy))
	s = f[1].data.shape

	x = np.arange(s[1]).repeat(s[2]).reshape(s[1],s[2])
	y = np.tile(np.arange(s[2]),s[1]).reshape(s[1],s[2])

	galaxy_gals = np.loadtxt('%s/Data/muse/analysis/galaxies.txt' % (cc.base_dir),
		unpack=True, skiprows=1, usecols=(0,), dtype=str)
	x_cent, y_cent = np.loadtxt('%s/Data/muse/analysis/galaxies.txt' % (cc.base_dir),
		unpack=True, skiprows=1, usecols=(1,2), dtype=int)
	i_gal = np.where(galaxy_gals==galaxy)[0][0]
	x_cent = x_cent[i_gal]
	y_cent = y_cent[i_gal]

	area = get_areaInAperture(s[1], s[2], x_cent, y_cent, app_size/res)

	if not inside:
		area = 1 - area

	# Deal with NaN
	d = f[1].data
	d[np.isnan(d)] = 0
	n = f[2].data
	n[np.isnan(n)] = 0
	return np.einsum('ijk,jk->i', d, area), np.einsum('ijk,jk->i', n, area), \
		np.arange(s[0])*f[1].header['CD3_3'] + f[1].header['CRVAL3']


def KDC_pop(galaxy):
	params = set_params(opt='pop')
	params.reps = 10

	spec, noise, lam = get_specFromAperture(galaxy, app_size=1.0)
	CD = lam[1] - lam[0]
	spec, lam, cut = apply_range(spec, window=201, repeats=3, 
		lam=lam, return_cuts=True, set_range=params.set_range, n_sigma=2)
	noise = noise[cut]
	lamRange = np.array([lam[0],lam[-1]])

	pp = run_ppxf(galaxy, spec, noise, lamRange, CD, params, produce_plot=False)

	pop = population(pp=pp, galaxy=galaxy)
	pop.plot_probability_distribution(label=' of core region')

	data_file = "%s/Data/muse/analysis/galaxies_core.txt" % (cc.base_dir)
	age_gals, age_unc_gals, met_gals, met_unc_gals, alp_gals, alp_unc_gals, \
		OIII_eqw_gals, OIII_eqw_uncer_gals, \
		age_gals_outer, age_unc_gals_outer, met_gals_outer, met_unc_gals_outer, \
		alp_gals_outer, alp_unc_gals_outer = np.loadtxt(data_file, unpack=True, 
		skiprows=2, usecols=(1,2,3,4,5,6,7,8,9,10,11,12,13,14))
	galaxy_gals = np.loadtxt(data_file, skiprows=2, usecols=(0,),dtype=str)
	i_gal = np.where(galaxy_gals==galaxy)[0][0]

	age_gals[i_gal] = pop.age
	age_unc_gals[i_gal] = pop.unc_age
	met_gals[i_gal] = pop.metallicity
	met_unc_gals[i_gal] = pop.unc_met
	alp_gals[i_gal] = pop.alpha
	alp_unc_gals[i_gal] = pop.unc_alp
 
	# Save plot from pop before clearing from memory
	f = pop.fig
	ax = pop.ax
	OIII_pos = np.argmin(np.abs(pp.lam - 5007))
	peak_width = 20
	OIII_pos += np.argmax(pop.e_line_spec[OIII_pos - peak_width:OIII_pos + peak_width]
		) - peak_width
	flux = np.trapz(pop.e_line_spec[OIII_pos - peak_width:OIII_pos + 
		peak_width], x=pp.lam[OIII_pos - peak_width:OIII_pos + peak_width])

	OIII_eqw_gals[i_gal] = flux/pop.continuum[OIII_pos]
	
	i_OIII = np.where('[OIII]5007d' in [e for e in pp.templatesToUse 
		if not e.isdigit()])[0][0]
	flux_uncert = trapz_uncert(pp.MCgas_uncert_spec[i_OIII, 
		OIII_pos - peak_width:OIII_pos + peak_width], 
		x=pp.lam[OIII_pos - peak_width:OIII_pos + peak_width])
	cont_uncert = np.sqrt(np.sum((pp.noise**2 + pp.MCgas_uncert_spec**2)[i_OIII, OIII_pos]))

	OIII_eqw_uncer_gals[i_gal] = np.sqrt(OIII_eqw_gals[i_gal]**2 * 
		((flux_uncert/flux)**2 + (cont_uncert/pop.continuum[OIII_pos])**2))
	del pop

	# Outside apperture
	spec, noise, lam = get_specFromAperture(galaxy, app_size=1.0, inside=False)
	CD = lam[1] - lam[0]
	spec, lam, cut = apply_range(spec, window=201, repeats=3, 
		lam=lam, return_cuts=True, set_range=params.set_range, n_sigma=2)
	noise = noise[cut]
	lamRange = np.array([lam[0],lam[-1]])

	

	pp_outside = run_ppxf(galaxy, spec, noise, lamRange, CD, params, produce_plot=False)
	pop_outside = population(pp=pp_outside, galaxy=galaxy)
	pop_outside.plot_probability_distribution(f=f, ax_array=ax, label=' of outer region')

	pop_outside.fig.suptitle('%s Probability Distribution within inner 1 arcsec' % (
		galaxy.upper()), y=0.985)
	h, l = pop_outside.ax[0,0].get_legend_handles_labels()
	pop_outside.ax[1,1].legend(h, l, loc=1)
	pop_outside.fig.savefig('%s/Data/muse/analysis/%s/pop_1arcsec.png' % (cc.base_dir, 
		galaxy))

	age_gals_outer[i_gal] = pop_outside.age
	age_unc_gals_outer[i_gal] = pop_outside.unc_age
	met_gals_outer[i_gal] = pop_outside.metallicity
	met_unc_gals_outer[i_gal] = pop_outside.unc_met
	alp_gals_outer[i_gal] = pop_outside.alpha
	alp_unc_gals_outer[i_gal] = pop_outside.unc_alp
	del pop_outside

	temp = "{0:10}{1:6}{2:6}{3:6}{4:6}{5:6}{6:6}{7:9}{8:10}{9:6}{10:6}{11:6}{12:6}{13:6}{14:6}\n"
	with open(data_file, 'w') as f:
		f.write('          Core (inner 1arcsec)                                   Outer \n')
		f.write(temp.format('Galaxy', 'Age', 'error', 'Metal', 'error', 'Alpha', 'error', 
			'OIII_eqw', 'error', 'Age', 'error', 'Metal', 'error', 'Alpha', 'error'))
		for i in range(len(galaxy_gals)):
			f.write(temp.format(galaxy_gals[i], str(round(age_gals[i],2)), 
				str(round(age_unc_gals[i],2)), str(round(met_gals[i],2)), 
				str(round(met_unc_gals[i],2)), str(round(alp_gals[i],2)), 
				str(round(alp_unc_gals[i],2)), str(round(OIII_eqw_gals[i],4)),
				str(round(OIII_eqw_uncer_gals[i],4)),
				str(round(age_gals_outer[i],2)), str(round(age_unc_gals_outer[i],2)), 
				str(round(met_gals_outer[i],2)), str(round(met_unc_gals_outer[i],2)),
				str(round(alp_gals_outer[i],2)), str(round(alp_unc_gals_outer[i],2))))





if __name__ == '__main__':
	for gal in ['ic1459', 'ic4296', 'ngc1316', 'ngc1399']:
		print gal
		KDC_pop(gal)
	# KDC_pop('ic4296')

