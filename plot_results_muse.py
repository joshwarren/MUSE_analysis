## ==================================================================
## 		Plot the outputted dynamics maps
## ==================================================================
## warrenj 20150330 Process to plot the results of pPXF and GANDALF 
## routines.
## warrenj 20150727 Changing to a python script
## warrenj 20150917 Altered to plot and save all 8 plots.
## warrenj 20151216 Added section to plot residuals.
## warrenj 20160111 Add section to plot histgrams of the fields.
## warrenj 20160405 Added keyword CO to overlay CO maps from ALMA if avaible.
## This supersedes plot_results_CO.py
## warrenj 20160726 Now plots in a more object orientated way and creates a
## grid of plots too. This supersedes plot_results2.py

## *************************** KEYWORDS ************************* ##
# galaxy 		Name of the galaxy being plotted: used to find 
#				correct files and to print onto the plot.
# discard	0	Interger giving the number of rows and columns 
#				to be removed from the plot to remove edge 
#				effects.
# norm		"lwv"	Normalisation methods for velocity fields:
#				lwv: luminosity weighted mean of the whole 
#				field is set to 0.
#				lum: velocity of the brightest spaxel is set 
#				to 0.
#				sig: Noralised to the mean velocity of 5 bins with the
#				highest velocity dispersion.
#				lws: Normalised to the mean velocity of 5 bins with the highest
#				luminosity (flux) weighted velocity dispersion
# plots 	False   Boolean to show plots as routine runs.
# nointerp 	False 	Boolean to use interpolation between bins in 
#				plots or not.
# residual 	False	Method to measure the residuals:
#			mean: use the mean of the residuals in each 
#				bin.
#			median: use the median of the residuals in 
#				each bin.
#			max: use the maximum of the residuals in 
#				each bin.
#			False: do not calculate and produce plot of 
#				residuals.
# overplot   {}	Dictionary containing name of thing to overplot and its color.
# D 		None Option to pass in the Data object instead of loading it.
## ************************************************************** ##
import numpy as np # for array handling
import glob # for searching for files
from astropy.io import fits # reads fits files (is from astropy)
from astropy import units as u
from astropy.coordinates import SkyCoord
from matplotlib.collections import LineCollection
from checkcomp import checkcomp
cc = checkcomp()
if 'home' not in cc.device:
	import matplotlib # 20160202 JP to stop lack-of X-windows error
	matplotlib.use('Agg') # 20160202 JP to stop lack-of X-windows error
import matplotlib.pyplot as plt # used for plotting
from plot_velfield_nointerp import plot_velfield_nointerp
from plot_histogram import plot_histogram
import os
from sauron_colormap2 import sauron2 as sauron
import cPickle as pickle
from errors2_muse import get_dataCubeDirectory



# Give axes a saveTo property
plt.axes.saveTo = property(lambda self:str())
# Give axes an x and y on figure grid property
plt.axes.figx = property(lambda self:int())
plt.axes.figy = property(lambda self:int())
# give axes a property to hold a colorbar axes
plt.axes.cax = property(lambda self:plt.axes())
# give axes a property to hold 2 additional axes for showing other axis
plt.axes.ax2 = property(lambda self:plt.axes())
plt.axes.ax3 = property(lambda self:plt.axes())

vin_dir = '%s/Data/muse/analysis' % (cc.base_dir)
vin_dir_cube = '%s/Data/muse' % (cc.base_dir)
out_dir = '%s/Data/muse/analysis' % (cc.base_dir)

#-----------------------------------------------------------------------------
class mapping(object):
	def __init__(self):
		self.SNR = True
		self.image = True
		self.equivalent_width = True
		self.amp_noise = True
		self.kinematics = True
		self.plot_resid = True
		self.line_ratios = True
#-----------------------------------------------------------------------------

#-----------------------------------------------------------------------------
def set_lims(v, positive=False, symmetric=False):
	if all(~np.isfinite(v)):
		return 0, 0

	v = v[np.isfinite(v)]

	for i in range(2):
		av = np.median(v)
		std = np.std(v)

		include = (v >= av - 3*std) * (v <= av + 3*std)
		v = v[include]

	vmin, vmax = min(v), max(v)

	if symmetric:
		vmax = np.mean([vmax, abs(vmin)])
		vmin = -vmax

	if positive:
		vmin = max(vmin, 0)
		vmax = max(vmax, 0)

	return vmin, vmax
#-----------------------------------------------------------------------------

#-----------------------------------------------------------------------------
def add_R_e(ax, galaxy, discard=0, pa=0):
	from classify import get_R_e
	from  matplotlib.patches import Ellipse
	R_e = get_R_e(galaxy)
	
	data_file =  "%s/galaxies.txt" % (vin_dir)
	x_cent_gals, y_cent_gals = np.loadtxt(data_file, unpack=True, skiprows=1, 
		usecols=(1,2), dtype=int)
	galaxy_gals = np.loadtxt(data_file, skiprows=1, usecols=(0,),dtype=str, unpack=True)
	i_gal = np.where(galaxy_gals==galaxy)[0][0]
	x_cent_pix = x_cent_gals[i_gal]
	y_cent_pix = y_cent_gals[i_gal]

	xlims = ax.get_xlim()
	ylims = ax.get_ylim()

	x_cent = xlims[0] + (xlims[1] - xlims[0])/(40-discard*2)*x_cent_pix
	y_cent = ylims[0] + (ylims[1] - ylims[0])/(40-discard*2)*y_cent_pix



	# data_file =  '%s/galaxies2.txt' % (vin_dir)
	data_file =  '%s/Data/vimos/analysis/galaxies2.txt' % (cc.base_dir)
	ellip_gals, pa_gals = np.loadtxt(data_file, unpack=True, skiprows=1, 
		usecols=(2,3), dtype=float)
	galaxy_gals = np.loadtxt(data_file, skiprows=1, usecols=(0,),dtype=str, unpack=True)
	i_gal = np.where(galaxy_gals==galaxy)[0][0]
	ellip = ellip_gals[i_gal]
	pa += pa_gals[i_gal]

	if ax.RaDec:
		patch = Ellipse([x_cent, y_cent], R_e*(1-ellip)/60/60, R_e/60/60, angle=pa, 
			fill=False)
	else:
		patch = Ellipse([x_cent, y_cent], R_e*(1-ellip), R_e, angle=pa, fill=False)
	ax.add_patch(patch)
#-----------------------------------------------------------------------------

#-----------------------------------------------------------------------------
def add_(overplot, color, ax, galaxy, header, close=False):
	image_dir=getattr(get_dataCubeDirectory(galaxy), overplot)
	# Arcsec coords
	if os.path.exists(image_dir) and not ax.RaDec:
		f = fits.open(image_dir)[0]

		x = np.arange(f.header['NAXIS1'])*abs(f.header['CDELT1'])*60*60
		y = np.arange(f.header['NAXIS2'])*abs(f.header['CDELT2'])*60*60

		x -= max(x)/2
		y -= max(y)/2

		# Coordinates of VIMOS pointing
		mcoord = SkyCoord(header['CRVAL1'], header['CRVAL2'],
			unit=(u.deg, u.deg))

		# Coordinates of ALMA pointing
		imcoord = SkyCoord(f.header['CRVAL1'], f.header['CRVAL2'],
			unit=(u.deg, u.deg))

		# Offset between the two pointings
		x -= ((mcoord.ra.degree - header['CRPIX1']*header['CD1_1']/(60*60)) -
			(imcoord.ra.degree +
			f.header['CRPIX1']*f.header['CDELT1']/(60*60)))*60*60
				
		y += ((mcoord.dec.degree - header['CRPIX2']*header['CD2_2']/(60*60)) -
			(imcoord.dec.degree +
			f.header['CRPIX2']*f.header['CDELT2']/(60*60)))*60*60

	# RA and dec coords
	elif os.path.exists(image_dir):
		f = fits.open(image_dir)[0]

		imcoord = SkyCoord(f.header['CRVAL1'], f.header['CRVAL2'],
			unit=(u.deg, u.deg))

		# Centered by eye
		if galaxy == 'ngc1316':
			f.header['CRPIX1'] = 272
			f.header['CRPIX2'] = 457

		x = (np.arange(f.header['NAXIS1']) - f.header['CRPIX1']) *\
			f.header['CDELT1'] + imcoord.ra.degree
		y = (np.arange(f.header['NAXIS2'])-f.header['CRPIX2']) *\
			f.header['CDELT2'] + imcoord.dec.degree

	# Plot and save
	if os.path.exists(image_dir):
		#remove random extra dimenisons.
		s = f.data.shape
		image = np.sum(f.data, axis=(0,1))
		# image = np.sum(f.data, axis=np.where(np.array(s)==1)[0])
		xlim = ax.get_xlim()
		ylim = ax.get_ylim()

		# Discard noise from outer parts of the galaxy -  for radio
		if overplot == 'radio':
			lim = np.nanmean(image) + np.nanstd(image)
			image[image < lim] = lim
		image = np.log(image)

		cs = ax.contour(x, y, image, colors=color, linestyles='solid', linewidth=1)
		cs.collections[0].set_label(overplot)

		ax.set_xlim(xlim)
		ax.set_ylim(ylim)

		plt.legend(facecolor='w')

		saveTo = os.path.dirname(ax.saveTo)+"/Overplot/" + \
			os.path.basename(ax.saveTo)
		if not os.path.exists(os.path.dirname(saveTo)):
			os.makedirs(os.path.dirname(saveTo))
		plt.savefig(saveTo, bbox_inches="tight")

		if close:
			plt.close()
		else:
			# Make lines thinner for pdf by finding the line objects
			for o in ax.get_children():
				if type(o) is LineCollection:
					o.set_linewidth(0.3)
#-----------------------------------------------------------------------------

#-----------------------------------------------------------------------------
def plot_results(galaxy, discard=0, norm="lwv", plots=False, residual=False, 
	overplot={}, show_bin_num=False, D=None, mapping=None, opt='kin'):	

	pa = {'ic1459':0, 'ic4296':0, 'ngc1316':0, 'ngc1399':0}

	res = 0.2 # arcsec (spatial MUSE resolution)
	pa = pa[galaxy] # PA from reduction

	data_file =  "%s/galaxies.txt" % (vin_dir)
	# different data types need to be read separetly
	file_headings = np.loadtxt(data_file, dtype=str)[0]
	col = np.where(file_headings=='SN_%s' % (opt))[0][0]

	SN_target_gals = np.loadtxt(data_file, unpack=True, skiprows=1, 
		usecols=(col))
	galaxy_gals = np.loadtxt(data_file, skiprows=1, usecols=(0,),dtype=str)
	i_gal = np.where(galaxy_gals==galaxy)[0][0]
	SN_target=SN_target_gals[i_gal]-10

	data_file =  "%s/Data/vimos/analysis/galaxies.txt" % (cc.base_dir)
	z_gals = np.loadtxt(data_file, unpack=True, skiprows=1, 
		usecols=(1))
	galaxy_gals = np.loadtxt(data_file, skiprows=1, usecols=(0,),dtype=str)
	i_gal = np.where(galaxy_gals==galaxy)[0][0]
	z = z_gals[i_gal]

	# Remove option if overplot file does not exist.
	if overplot:
		for o in overplot.keys():
			if not getattr(get_dataCubeDirectory(galaxy),o):
				del overplot[o]

	dataCubeDirectory = get_dataCubeDirectory(galaxy)
	output = "%s/%s/%s" % (out_dir, galaxy, opt)
	out_plots = "%s/plots" % (output)
	out_nointerp = "%s/notinterpolated" % (out_plots)
	vin_dir_gasMC = "%s/%s/%s/MC" % (vin_dir, galaxy, opt) # for chi2
	out_pickle = '%s/pickled' % (output)

	cubeFile = fits.open(dataCubeDirectory)
	header = cubeFile[1].header
	cubeFile.close()
# ------------== Reading pickle file and create plot  ===----------

	# Load pickle file from pickler.py
	if D is None:
		pickleFile = open("%s/dataObj.pkl" % (out_pickle), 'rb')
		D = pickle.load(pickleFile)
		pickleFile.close()

	for bin in D.bin:
		for e in bin.e_line.itervalues():
			e.__threshold__ = 3

	if D.norm_method != norm:
		D.norm_method = norm
		D.find_restFrame()

	# Adjust by hand
	if galaxy == 'ic1459' and norm == 'lws':
		D.vel_norm -= 15
	if galaxy == 'ic4296' and norm == 'lws':
		D.vel_norm += 20
	if galaxy == 'ngc1399' and norm =='lws':
		D.vel_norm += 35

	# Create figure and array for axes
	n_rows = 2+2*len(D.e_components) + int(np.ceil(len(D.e_components)*
		(len(D.e_components)-1)/6.0))
	f = plt.figure(frameon=False)
	ax_array = []


	if mapping.SNR or mapping is None:	
		saveTo = "%s/SNR.png" % (out_nointerp)
		ax1 = plot_velfield_nointerp(D.x, D.y, D.bin_num, D.xBar, D.yBar, D.SNRatio, 
			colorbar=True, nodots=True, title='SNR', save=saveTo, close=not overplot=={}, 
			res=res)
# ------------=============== Plot image ================----------
	if mapping.image or mapping is None:
		print "    Image"
		
		title = "Total Flux"
		CBLabel = r"Flux (erg s$^{-1}$ cm$^{-2}$)"

		ax = f.add_subplot(111, aspect='equal')
		saveTo = "%s/total_image.png" % (out_nointerp)
		ax.saveTo = saveTo
		ax.figx, ax.figy = 0, 0

		fmin, fmax = set_lims(D.flux, positive=True)

		ax = plot_velfield_nointerp(D.x, D.y, D.bin_num, D.xBar, D.yBar, D.flux, 
			vmin=fmin, vmax=fmax, nodots=True, show_bin_num=show_bin_num, colorbar=True, 
			label=CBLabel, title=title, cmap='gist_yarg', ax=ax, res=res,
			flux_unbinned=D.unbinned_flux, header=header)
		if overplot:
			for o, c in overplot.iteritems():
				add_(o, c, ax, galaxy, header, close=True)
		ax_array.append(ax)
		f.delaxes(ax)
		f.delaxes(ax.cax)
		if hasattr(ax,'ax2'): f.delaxes(ax.ax2)
		if hasattr(ax,'ax3'): f.delaxes(ax.ax3)
		
		if plots:
			plt.show()
# ------------========= Plot intensity (& EW) ===========----------
	if mapping.equivalent_width or mapping is None:
		print "    gas map(s) and equivalent widths"
		for c in D.e_components:
			print "        " + c

			if 'OIII' in c:
				c_title = '[OIII]'
			elif 'Hbeta' in c:
				c_title = r'H$_\beta$'
			elif 'Hgamma' in c:
				c_title = r'H$_\gamma$'
			else:
				c_title = c

			f_title = "%s Flux" % (c_title)
			fh_title = "%s Flux Histogram" % (c_title)
			# from header
			fCBtitle = r"Flux (erg s$^{-1}$ cm$^{-2}$)"
			f_min, f_max = set_lims(D.e_line[c].flux, positive=True)

			saveTo = "%s/%s_flux_hist.png" % (out_plots, c)
			plot_histogram(D.e_line[c].flux, galaxy=galaxy.upper(), redshift=z,
				vmin=f_min,vmax=f_max, weights=D.n_spaxels_in_bin, title=fh_title,
				xaxis=fCBtitle, save=saveTo)
			
			ax = f.add_subplot(111, aspect='equal')
			saveTo = "%s/%s_img.png" % (out_nointerp, c)
			ax.saveTo = saveTo
			
			ax = plot_velfield_nointerp(D.x, D.y, D.bin_num, D.xBar, D.yBar, 
				D.e_line[c].flux, vmin=f_min, vmax=f_max, colorbar=True, nodots=True, 
				label=fCBtitle, title=f_title, cmap = 'gist_yarg', ax=ax, res=res,
				flux_unbinned=D.unbinned_flux)
			ax_array.append(ax)
			f.delaxes(ax)
			f.delaxes(ax.cax)
			if hasattr(ax,'ax2'): f.delaxes(ax.ax2)
			if hasattr(ax,'ax3'): f.delaxes(ax.ax3)
				
			if plots: plt.show()
			

			eq_title = "%s Equivalent Width" % (c_title)
			eqh_title = "%s Equivalent Width Histogram" % (c_title)
			eqCBtitle = r"Equivalent Width ($\AA$)"

			eq_min, eq_max = set_lims(D.e_line[c].equiv_width)#, positive=True)

			saveTo = "%s/%s_eqWidth_hist.png" % (out_plots, c)
			plot_histogram(D.e_line[c].equiv_width, galaxy=galaxy.upper(), redshift=z,
				vmin=eq_min,vmax=eq_max, weights=D.n_spaxels_in_bin, title=eqh_title,
				xaxis=eqCBtitle, save=saveTo)
			
			ax = f.add_subplot(111, aspect='equal')
			saveTo = "%s/%s_equiv_width.png" % (out_nointerp, c)
			ax.saveTo = saveTo


			ax = plot_velfield_nointerp(D.x, D.y, D.bin_num, D.xBar, D.yBar, 
				D.e_line[c].equiv_width, vmin=eq_min, vmax=eq_max, colorbar=True, 
				nodots=True, label=eqCBtitle, title=eq_title, ax=ax, res=res,
				flux_unbinned=D.unbinned_flux, signal_noise=D.SNRatio,
				signal_noise_target=SN_target)
			ax_array.append(ax)
			f.delaxes(ax)
			f.delaxes(ax.cax)
			if hasattr(ax,'ax2'): f.delaxes(ax.ax2)
			if hasattr(ax,'ax3'): f.delaxes(ax.ax3)
# ------------============ Amplitude/Noise ==============----------
	if mapping.amp_noise or mapping is None:
		for c in D.e_components:
			amp_title = '%s Amplitude to Noise ratio' % (c_title)
			amp_min, amp_max = set_lims(D.e_line[c].amp_noise, positive=True)
			saveTo = "%s/%s_amp_nosie.png" % (out_nointerp, c)

			ax1 = plot_velfield_nointerp(D.x, D.y, D.bin_num, D.xBar, D.yBar, 
				D.e_line[c].amp_noise, vmin=amp_min, vmax=amp_max, colorbar=True, 
				nodots=True, title=amp_title, save=saveTo, close=not overplot=={}, 
				res=res, flux_unbinned=D.unbinned_flux)
			if overplot:
				ax1.saveTo = saveTo
				for o, c in overplot.iteritems():
					add_(o, c, ax1, galaxy, header, close=True)
# ------------=========== Setting titles etc ============----------
	if mapping.kinematics or mapping is None:
		print '    Kinematics'
		# for c in ['stellar']: # For debugging
		for c in D.independent_components:
			print '        %s' % (c)

			im_type = c
			pl = c
			if im_type == "gas":
				im_type=""
				pl = 'Hbeta'
			elif im_type == "SF":
				im_type=" (Star Forming)"
				pl = '[OIII]5007d'
			elif im_type == "Shocks":
				im_type=" (Shocking)"
				pl = 'Hbeta'
			elif 'Hbeta' in im_type:
				im_type=" ("+r'H$_\beta$'+")"
			elif 'Hgamma' in im_type:
				im_type=" ("+r'H$_\gamma$'+")"
			elif 'OIII' in im_type:
				im_type=" (OIII)"
			else:
				im_type=" (" + im_type + ")"

			for k in D.components[pl].plot.keys():


				symmetric=False
				positive=False
					
				CBLabel = None
				if k == "vel":
					title = 'Velocity'
					CBLabel = "V (km s$^{-1}$)"
					symmetric=True

				if  k == "sigma":
					title = 'Velocity Dispersion'
					CBLabel = r'$\mathrm{\sigma}$ (km s$^{-1}$)'
					positive = True

				if k == "h3":
					title = 'h3'
					symmetric = True

				if k == "h4":
					title = 'h4'


				if c == "stellar":
					utitle = "Stellar Uncertainty " + title + " Map"
					htitle = "Stellar " + title + " Histogram"
					uhtitle = "Stellar Uncertainty " + title + " Histogram"
					title = "Stellar " + title + " Map"
				else:
					utitle = "Ionised" + im_type + " Gas Uncertainty " + title + " Map"
					htitle = "Ionised" + im_type + " Gas " + title + " Histogram"
					uhtitle = "Ionised" + im_type + " Gas Uncertainty " + title + \
						" Histogram"
					title = "Ionised" + im_type + " Gas\n" + title + " Map"
# ------------============ Setting v range ==============----------
				vmin, vmax = set_lims(D.components[pl].plot[k], positive=positive, 
					symmetric=symmetric)
				v_uncert_min, v_uncert_max = set_lims(D.components[pl].plot[k].uncert, 
					positive=True)
# # ------------============== Plot Histogram =============----------
				# Field histogram
				# saveTo = "%s/%s_hist_%s.png" % (out_plots, plot_title)
				# plot_histogram(D.components[c].plot[k], galaxy=galaxy.upper(), 
				# 	redshift=z, vmin=vmin,vmax=vmax, weights=D.n_spaxels_in_bin, 
				# 	title=htitle, xaxis=CBLabel, save=saveTo)
				# # Uncertainty histogram
				# saveTo = "%s/%s_hist_%s.png" % (out_plots, plot_title+'_uncert', 
				# 	wav_range)
				# plot_histogram(D.components[c].plot[k].uncert, galaxy=galaxy.upper(), 
				# 	redshift=z, vmin=v_uncert_min,vmax=v_uncert_max, 
				# 	weights=D.n_spaxels_in_bin, title=uhtitle, xaxis=CBLabel, save=saveTo)

				# if plots:
				# 	plt.show()
# ------------==== Plot velfield - no interperlation ====----------
				# Field plot
				ax = f.add_subplot(111, aspect='equal')
				saveTo = ("%s/%s_%s_field.png" % (out_nointerp, c, k))
				ax.saveTo = saveTo
				ax = plot_velfield_nointerp(D.x, D.y, D.bin_num, D.xBar, D.yBar,
					D.components[pl].plot[k], vmin=vmin, vmax=vmax, 
					flux_unbinned=D.unbinned_flux, #flux_type='notmag',
					nodots=True, show_bin_num=show_bin_num, colorbar=True, 
					label=CBLabel,galaxy = galaxy.upper(), redshift = z,
					title=title, ax=ax, res=res, signal_noise=D.SNRatio,
					signal_noise_target=SN_target, header=header)
				if overplot:
					for o, c in overplot.iteritems():
						add_(o, c, ax, galaxy, header, close=True)
				# add_R_e(ax, galaxy, pa=pa)
				if plots:
					plt.show()
				ax_array.append(ax)
				f.delaxes(ax)
				f.delaxes(ax.cax)
				if hasattr(ax,'ax2'): f.delaxes(ax.ax2)
				if hasattr(ax,'ax3'): f.delaxes(ax.ax3)
				
				# Uncertainty plot
				# saveTo = "%s/%s_%s_uncert_field.png" % (out_nointerp, c, k)
				# ax1 = plot_velfield_nointerp(D.x, D.y, D.bin_num, D.xBar, D.yBar,
				# 	D.components[pl].plot[k].uncert, vmin=v_uncert_min, 
				# 	vmax=v_uncert_max, #flux_type='notmag', 
				# 	nodots=True, show_bin_num=show_bin_num, colorbar=True, 
				# 	label=CBLabel, galaxy = galaxy.upper(), redshift = z, 
				# 	title=utitle, save=saveTo, close=not overplot=={}, res=res)
				# 	#, header=header)
				# if overplot:
				# 	ax1.saveTo = saveTo
				# 	for o, c in overplot.iteritems():
				# 		add_(o, c, ax1, galaxy, header, close=True)
				if plots:
					plt.show()
# ------------============= Plot residuals ==============----------
	if residual and (mapping.plot_resid or mapping is None):
		print "    " + residual + " residuals"

		average_residuals = np.zeros(D.number_of_bins)
		for i, bin in enumerate(D.bin):
			residuals = np.abs(bin.spectrum - bin.bestfit)
			# remove edge pixels
			residuals = np.delete(residuals, [np.arange(5), 
				len(residuals)+np.arange(-5,0)], axis=0)

			if residual=="mean":
				average_residuals[i] = np.mean(residuals)
			elif residual=="median":
				average_residuals[i] = np.median(residuals)
			elif residual=="max":
				average_residuals[i] = np.max(np.abs(residuals))
				
		minres, maxres = set_lims(average_residuals, positive=True) #mean_centered=True,
		
		CBLabel = "Residuals"
		title = str.capitalize(residual) + \
		" Residuals of Bestfit to Normalised Spectrum"
		saveTo = "%s/%s_residual.png" % (out_nointerp, residual)

		ax1 = plot_velfield_nointerp(D.x, D.y, D.bin_num, D.xBar, D.yBar,
			average_residuals, vmin=minres, vmax=maxres, flux_type='notmag',
			nodots=True, show_bin_num=show_bin_num, colorbar=True, 
			label=CBLabel, #flux_unbinned=D.unbinned_flux, 
			galaxy = galaxy.upper(), redshift = z, title=title, 
			save=saveTo, close=not overplot=={}, res=res)#, header=header)
		if plots:
			plt.show()
		if overplot:
			ax1.saveTo = saveTo
			for o, c in overplot.iteritems():
				add_(o, c, ax1, galaxy, header, close=True)
# # ------------=============== Plot Chi2/DOF =============----------
	# print "    chi2"

	# chi2 = np.zeros(D.number_of_bins)
	# for i in range(D.number_of_bins):
	# 	chi2[i] = np.loadtxt("%s/chi2/%d.dat" % (vin_dir_gasMC, i))

	# minchi2, maxchi2 = set_lims(chi2, positive = True)
	
	# CBLabel = "Chi2/DOF"
	# title = "Chi2/DOF of the bestfit"
	# saveTo = "%s/chi2_%s.png" % (out_nointerp)

	# ax1 = plot_velfield_nointerp(D.x, D.y, D.bin_num, D.xBar, D.yBar, chi2, 
	# 	vmin=minchi2, vmax=maxchi2, flux_type='notmag',
	# 	nodots=True, show_bin_num=show_bin_num, colorbar=True, 
	# 	label=CBLabel, flux_unbinned=D.unbinned_flux, 
	# 	galaxy = galaxy.upper(), redshift = z, title=title, 
	# 	save=saveTo, close=not overplot=={}, header=header)#, cmap=cm.blue)
	# if plots:
	# 	plt.show()
	# if overplot:
	# 	ax1.saveTo = saveTo
	# 	for o, c in overplot.iteritems():
	# 		add_(o, c, ax1, galaxy, header, close=True)
# ------------============ Line ratio maps ==============----------
	# if any('OIII' in o for o in D.list_components) and line_ratios:
	if len(D.list_components) > 2 and (mapping.line_ratios or mapping is None):
		print "    line ratios"

		t_num = (len(D.e_components)-1)*len(D.e_components)/2
		for n in range(t_num):
			i = 0
			m = t_num
			while m > n:
				i += 1
				m -= i

			cA = D.e_components[len(D.e_components)-i-1]
			cB = D.e_components[len(D.e_components)-i+n-m]

			line_ratio = np.log10(D.e_line[cB].flux/D.e_line[cA].flux)
			if 'OIII' in cA:
				cA_title = '[OIII]'
			elif 'Hbeta' in cA:
				cA_title = r'H$_\beta$'
			elif 'Hdelta' in cA:
				cA_title = r'H$_\delta$'
			elif 'Hgamma' in cA:
				cA_title = r'H$_\gamma$'
			else:
				cA_title = cA

			if 'OIII' in cB:
				cB_title = '[OIII]'
			elif 'Hbeta' in cB:
				cB_title = r'H$_\beta$'
			elif 'Hdelta' in cB:
				cB_title = r'H$_\delta$'
			elif 'Hgamma' in cB:
				cB_title = r'H$_\gamma$'
			else:
				cB_title = cB
				
			lr_title = "%s/%s Line Ratio" % (cB_title, cA_title)
			lrCBtitle = r"log$_{10}$ (%s/%s)" %(cB_title,cA_title)

			lr_min, lr_max = set_lims(line_ratio)


			ax = f.add_subplot(111, aspect='equal')
			saveTo = "%s/lineratio/%s_%s_line_ratio.png" % (out_nointerp, cB, cA)
			ax.saveTo = saveTo
			ax.figx, ax.figy = n%3, n_rows-int(np.ceil(t_num/3)) + int(np.ceil(n/3))


			ax = plot_velfield_nointerp(D.x, D.y, D.bin_num, D.xBar, D.yBar,
				line_ratio, vmin=lr_min, vmax=lr_max, colorbar=True,
				nodots=True, title=lr_title, label=lrCBtitle, ax=ax,
				show_bin_num=show_bin_num, galaxy = galaxy.upper(), redshift = z, 
				res=res)#,	header=header)

			ax_array.append(ax)
			f.delaxes(ax)
			f.delaxes(ax.cax)
			if hasattr(ax,'ax2'): f.delaxes(ax.ax2)
			if hasattr(ax,'ax3'): f.delaxes(ax.ax3)
# ------------============= Plot and save ===============----------
	print "    Plotting and saving"

	for i, a in enumerate(ax_array):
		f.add_axes(a)
		#a.axis('tight')
		f.add_axes(a.cax)
		if hasattr(a,'ax2'): f.add_axes(a.ax2)
		if hasattr(a,'ax3'): f.add_axes(a.ax3)
		if not os.path.exists(os.path.dirname(a.saveTo)):
			os.makedirs(os.path.dirname(a.saveTo))
		print a.get_title()
		plt.savefig(a.saveTo)#, bbox_inches="tight")

		# if overplot:
		# 	for o, c in overplot.iteritems():
		# 		add_(o, c, a, galaxy, header, close=True)

		f.delaxes(a)
		f.delaxes(a.cax)
		if hasattr(a,'ax2'): f.delaxes(a.ax2)
		if hasattr(a,'ax3'): f.delaxes(a.ax3)

	return D



##############################################################################

# Use of plot_results.py

if __name__ == '__main__':

	galaxies = ['ic1459', 'ic4296', 'ngc1316', 'ngc1399']
	galaxy = galaxies[3]

	print galaxy

	plot_results(galaxy, plots=False, residual = "median", show_bin_num=True)