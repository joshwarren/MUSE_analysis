# ==================================================================
# Rebinning using Voronoi tessellation
# ==================================================================
# warrenj 20150515 Process to rebin spaxels together in order to
# create a minimum S/N ratio. targetSN is approx 10-15 for just v and
# sigma, and 40-50 for h3 and h4. 
# warrenj 20160913 Ported to python

from checkcomp import checkcomp
cc=checkcomp()
if cc.remote:
	import matplotlib # 20160202 JP to stop lack-of X-windows error
	matplotlib.use('Agg') # 20160202 JP to stop lack-of X-windows error
import numpy as np
import glob
from astropy.io import fits
import ppxf_util as util
from voronoi_2d_binning import voronoi_2d_binning
from errors2_muse import remove_anomalies, get_dataCubeDirectory
import os


# ----------===============================================---------
# ----------======== Check overwrite of target SN =========---------
# ----------===============================================---------
def check_overwrite(new, old, auto_override=False):
	if not auto_override and new != old:
		A = raw_input('Are you sure you want to overwrite the old target ' + 
			'of %d with a new target of %d? (Y/N) ' % (old, new))
		if A == "N" or A == "n": new = old
	return new



def binning_spaxels(galaxy, targetSN=None, opt='kin', auto_override=False, debug=False,
	set_range=None):
	print '     Voronoi Binning'
# ----------===============================================---------
# ----------============ Default parameters ===============---------
# ----------===============================================---------
	dir = "%s/Data/muse" % (cc.base_dir)
	data_file = "%s/analysis/galaxies.txt" %(dir)
	# Check if file has anything in it - it does need to exsist.
	try:
		galaxy_gals = np.loadtxt(data_file, skiprows=1, usecols=(0,),dtype=str)
		x_gals, y_gals, SN_kin_gals, SN_pop_gals = np.loadtxt(data_file, skiprows=1, 
			usecols=(1,2,3,4), unpack=True, dtype='int,int,float,float')
	except StopIteration:
		galaxy_gals = np.array([])
		x_gals = np.array([])
		y_gals = np.array([])
		SN_kin_gals = np.array([])
		SN_pop_gals = np.array([])

	if opt=='kin':
		SN_used_gals = SN_kin_gals
	elif opt=='pop':
		SN_used_gals = SN_pop_gals

	i_gal = np.where(galaxy_gals == galaxy)[0]
	if len(i_gal) == 0:
		i_gal = -1
		galaxy_gals = np.append(galaxy_gals, galaxy)

	if targetSN is None and i_gal != -1:
		targetSN=SN_used_gals[i_gal]
	elif targetSN is not None and i_gal  != -1: 
		targetSN = check_overwrite(targetSN, SN_used_gals[i_gal], auto_override)
		SN_used_gals[i_gal] = targetSN
	elif targetSN is not None and i_gal == -1:
		SN_used_gals = np.append(SN_used_gals, targetSN)
	else:
		targetSN = 30.0
		SN_used_gals = np.append(SN_used_gals, targetSN)

	if i_gal == -1:
		x_gals = np.append(x_gals, 0)
		y_gals = np.append(y_gals, 0)
		if opt == 'kin':
			SN_pop_gals = np.append(SN_pop_gals, 0)
		elif opt == 'pop':
			SN_kin_gals = np.append(SN_kin_gals, 0)

# ----------================= Save SN_used ===============---------
	if opt=='kin':
		SN_kin_gals = SN_used_gals
	elif opt=='pop':
		SN_pop_gals = SN_used_gals 

	temp = "{0:12}{1:4}{2:4}{3:8}{4:8}\n"
	with open(data_file, 'w') as f:
		f.write(temp.format("Galaxy", "x", "y", "Kin SN", "Pop SN"))
		for i in range(len(galaxy_gals)):
			f.write(temp.format(galaxy_gals[i], str(int(x_gals[i])), 
				str(int(y_gals[i])), str(round(SN_kin_gals[i],2)),
				str(round(SN_pop_gals[i],2))))

# ----------================ Find S/N ================------------
# Final wildcard notes that depending on the method used the quadrants
#may or may not have been flux calibrated. 
	dataCubeDirectory = get_dataCubeDirectory(galaxy)

	f = fits.open(dataCubeDirectory)
	galaxy_data, header = f[1].data, f[1].header
	galaxy_noise = f[2].data

	## write key parameters from header - can then be altered in future	
	CRVAL_spec = header['CRVAL3']
	CDELT_spec = header['CD3_3']

	s = galaxy_data.shape

	x = np.zeros(s[1]*s[2])
	y = np.zeros(s[1]*s[2])

	if set_range is not None:
		set_range[0] = max(set_range[0], CRVAL_spec)
		set_range_pix = (set_range - CRVAL_spec)/CDELT_spec
	else:
		set_range_pix = np.array([0,s[0]])
	set_range_pix = set_range_pix.astype(int)

# collapsing the spectrum for each spaxel. 
	if debug:
		signal = np.array(galaxy_data[int(np.mean(set_range_pix)),:,:].flatten())
		# noise = np.sqrt(galaxy_noise[s[0]/2,:,:])#.flatten())
		noise = np.sqrt(np.abs(galaxy_data[int(np.mean(set_range_pix)),:,:])).flatten()
	else:
		signal = np.zeros((s[1],s[2]))
		noise = np.zeros((s[1],s[2]))
		blocks = 10
		bl_delt1 = int(np.ceil(s[1]/float(blocks)))
		bl_delt2 = int(np.ceil(s[2]/float(blocks)))
		for i in xrange(blocks):
			for j in xrange(blocks):
				signal[bl_delt1*i:bl_delt1*(i+1),bl_delt2*j:bl_delt2*(j+1)] = \
					np.nanmedian(galaxy_data[set_range_pix[0]:set_range_pix[1], 
					bl_delt1*i:bl_delt1*(i+1), bl_delt2*j:bl_delt2*(j+1)], axis=0)
				noise[bl_delt1*i:bl_delt1*(i+1),bl_delt2*j:bl_delt2*(j+1)] = \
					np.nanmedian(np.abs(galaxy_noise[set_range_pix[0]:set_range_pix[1], 
					bl_delt1*i:bl_delt1*(i+1), bl_delt2*j:bl_delt2*(j+1)]), axis=0)

		signal = signal.flatten()
		noise = noise.flatten()

		bad_pix = ~(signal > 0) + ~(noise > 0)
		signal[bad_pix] = 0
		noise[bad_pix] = 0
		noise +=0.000001


	galaxy_data = []
	del galaxy_data
	galaxy_noise = []
	del galaxy_noise


	for i in range(s[1]):
		for j in range(s[2]):
			# Assign x and y
			x[i*s[2]+j] = i
			y[i*s[2]+j] = j

	mask = (np.isfinite(signal)) * (np.isfinite(noise))

	nobin = signal/noise > targetSN*2

	signal = signal[mask + ~nobin]
	noise = noise[mask + ~nobin]
	x = x[mask + ~nobin]
	y = y[mask + ~nobin]
	n_spaxels = np.sum(mask) # include the not-for-binning-bins

	# signal = signal[mask]
	# noise = noise[mask]
	# x = x[mask]
	# y = y[mask]
	# n_spaxels = np.sum(mask)

	if not os.path.exists("%s/analysis/%s" % (dir,galaxy)):
		os.makedirs("%s/analysis/%s" % (dir, galaxy))

	# if not debug:
	binNum, xNode, yNode, xBar, yBar, sn, nPixels, scale = voronoi_2d_binning(
		x, y, signal, noise, targetSN, quiet=True, plot=False,
		saveTo='%s/analysis/%s/binning_%s.png' %(dir,galaxy, opt))
	matplotlib.pyplot.close('all')
	# else:
	# 	binNum = np.arange(len(x))
	# 	xBar = x
	# 	yBar = y

	xBar = np.append(xBar, x[nobin])
	yBar = np.append(yBar, y[nobin])
	binNum = np.append(binNum, np.arange(np.sum(nobin))+max(binNum)+1)

	order = np.argsort(binNum)
	xBin = np.zeros(n_spaxels)
	yBin = np.zeros(n_spaxels)

	# spaxel number
	i = 0
	for bin in range(max(binNum)+1):
		while i < n_spaxels and bin == binNum[order[i]]:
			xBin[order[i]] = xBar[bin]
			yBin[order[i]] = yBar[bin]
			# move onto next spaxel
			i = i + 1

# ------------================ Saving Results ===============---------------			

	temp = "{0:5}{1:5}{2:8}{3:9}{4:9}\n"
	temp2 = "{0:12}{1:12}\n"

	with open("%s/analysis/%s/voronoi_2d_binning_output_%s.txt" % (dir,galaxy,opt), 
		'w') as f:
		f.write(temp.format('X"', 'Y"', 'BIN_NUM', 'XBIN', 'YBIN'))
		for i in range(len(xBin)):
			f.write(temp.format(str(int(x[i])), str(int(y[i])), str(int(binNum[i])), 
				str(round(xBin[i],5)), str(round(yBin[i],5))))


	with open("%s/analysis/%s/voronoi_2d_binning_output2_%s.txt" % (dir,galaxy,opt), 
		'w') as f:
		f.write(temp2.format('XBAR','YBAR'))
		for i in range(len(xBar)):
			f.write(temp2.format(str(round(xBar[i],5)), str(round(yBar[i],5)))) 

	print 'Number of bins: ', max(binNum)+1