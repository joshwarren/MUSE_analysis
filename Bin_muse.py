## ==================================================================
## 		Emission lines
## ==================================================================
## warrenj 20160725 An object to hold emission lines. 

import numpy as np
import ppxf_util as util
from absorption import absorption
from glob import glob
import os
from checkcomp import checkcomp
cc = checkcomp()
from errors2_muse import get_dataCubeDirectory


c = 299792.458 # speed of light in km/s

vin_dir = '%s/Data/muse/analysis' % (cc.base_dir)
out_dir = '%s/Data/muse/analysis' % (cc.base_dir)

lines = {'Hdelta':4101.76, 'Hgamma':4340.47, 'Hbeta':4861.33, '[OIII]5007d':5004.0,
	'[NI]d':5200.0, 'Halpha':6562.8, '[SII]6716':6716.0, '[SII]6731':6731.0, 
	'[OI]6300d':6300.0, '[NII]6583d':6583.0}
dynamics = np.array(['vel', 'sigma', 'h3', 'h4'])

class Data(object):
# Attributes:
# number_of_bins: int
# bin: array of bin objects (see below) - holds most of the data - is a 
#	'setter' class
# e_line: dictionary of emission_data objects (see below) - is a 'getter' class
# e_components: string list of emission lines
# unbinned_flux: float array of unbinned flux in 2D array
# flux: float array of total integrated flux (observed - not fitted)
# xBar: float array of central point's x-coord for each bin
# yBar: as xBar
# n_spaxels_in_bin: int array of number of spaxels in each bin
# norm_method: str, normalisation methods for velocity fields:
#		lwv: 	(default) luminosity weighted mean of the whole 
#				field is set to 0.
#		lum: 	velocity of the brightest spaxel is set 
#				to 0.
#		sig: 	Noralised to the mean velocity of 5 bins with the
#				highest velocity dispersion.
# vel_norm: float normalisation factor for finding the rest frame of the galaxy.
# center_bin: int, bin number of the brightest bin.
# common_range: array of min and max of common wavelength range
# galaxy_spectrum: array of spectrum of whole galaxy
# galaxy_continuum: array of spectrum of whole galaxy with emission lines 
#	removed.
# SNRatio: (array) Signal to Noise Ratio of each bin. 
# gas_opt: 	0	(default) No emission lines
#			1	All emission lines set with LOSVD fixed to each other
#			2	All Balmers lines set with the same LOSVD, and all others set with fixed
#				LOSVD ([OIII] and [NI])
#			3	All emission lines set with independant LOSVD
# independent_components: list components with inderpendant LOSVDs.
#
# Methods:
# add_e_line (line name, line wavelength): adds a new emission line object 
#	to e_line.
# set_spaxels_in_bins (spaxel x-coord, spaxel y-coord, bin membership of spaxel): 
#	sets which bin contains which spaxel.
# absorption_line (absorption line): returns absorption line indice level
# 	from Lick like methods.
	def __init__(self, galaxy, opt):
		self.galaxy = galaxy
		self.opt = opt

		# Files
		self.tessellation_File = "%s/%s/voronoi_2d_binning_output_%s.txt" % (vin_dir, 
			galaxy, opt)
		self.tessellation_File2 = "%s/%s/voronoi_2d_binning_output2_%s.txt" %(vin_dir, 
			galaxy, opt)
		self.dataCubeDirectory = get_dataCubeDirectory(galaxy)
		if self.opt == 'kin':
			self.vin_dir_gasMC = '%s/%s/gas_MC' % (vin_dir, galaxy)
		else:
			self.vin_dir_gasMC = '%s/%s/pop_MC' % (vin_dir, galaxy)

		self.number_of_bins = int(max(self.bin_num) + 1) # zero-base
		self.norm_method = 'lwv'
		self.vel_norm = 0.0
		self.common_range = self.find_common_range()

		self.find_restFrame()


		# Check tessellation file is older than pPXF outputs (checks 
		# self.vin_dir_gasMC/0.dat only).
		if os.path.getmtime(self.tessellation_File) > os.path.getmtime('%s/0.dat' % (
			self.vin_dir_gasMC)): 
			if os.path.exists('%s/%i.dat' % (self.vin_dir_gasMC,max(self.bin_num))) and \
				not os.path.exists('%s/%i.dat' % (self.vin_dir_gasMC, 
				max(self.bin_num)+1)):
				# Issue warning, but do not stop.
				import warnings
				warnings.warn('WANING: The tesselation file '+\
					'voronoi_2d_binning_output_%s.txt may to have been changed.' % (opt))
			else:
				# Stop and raise exception
				raise UserWarning('WANING: The tesselation file '+\
					'voronoi_2d_binning_output_%s.txt has been overwritten.' % (opt))

	@property
	def x(self):
		return np.loadtxt(self.tessellation_File, unpack=True, skiprows = 1, 
			usecols=(0,), dtype=int)

	@property
	def y(self):
		return np.loadtxt(self.tessellation_File, unpack=True, skiprows = 1, 
			usecols=(1,), dtype=int)

	@property
	def bin_num(self):
		return np.loadtxt(self.tessellation_File, unpack=True, skiprows = 1, 
			usecols=(2,), dtype=int)

	@property
	def xBar(self):
		return np.loadtxt(tessellation_File2, unpack=True, skiprows = 1, usecols=(0,))

	@property
	def yBar(self):
		return np.loadtxt(tessellation_File2, unpack=True, skiprows = 1, usecols=(1,))

	@property
	def unbinned_flux(self):
		from astropy.io import fits
		galaxy_data = fits.getdata(self.dataCubeDirectory, 0)
		return np.nansum(galaxy_data, axis=0)

	@property
	def bin(self):
		return [Bin(i,self) for i in range(self.number_of_bins)]

	@property
	def center_bin(self):
		return np.nanargmax(self.flux)

	@property
	def gas(self):
		componants = [d for d in os.listdir(self.vin_dir_gasMC + "/gas") if \
			os.path.isdir(os.path.join(self.vin_dir_gasMC + "/gas", d))]

		if len(componants) == 0: return 0
		elif 'gas' in componants: return 1
		elif 'Shocks' in componants and 'SF' in componants: return 2
		else: return 3

	@property
	def independent_components(self):
		if self.gas == 0: 
			return ['stellar']
		elif self.gas == 1: 
			return ['stellar', 'gas']
		elif self.gas == 2: 
			return ['stellar', 'SF', 'Shocks']
		elif self.gas == 3:
			return self.list_components 

	@property
	def n_spaxels_in_bin(self):
		return np.array([bin.n_spaxels_in_bin for bin in self.bin])

	@property
	def flux(self):
		return np.array([bin.flux for bin in self.bin])

	@property
	def galaxy_spectrum(self):
		s = np.zeros(len(self.bin[0].spectrum))
		for bin in self.bin:
			s += bin.spectrum
		return s

	@property
	def galaxy_continuum(self):
		s = np.zeros(len(self.bin[0].continuum))
		for bin in self.bin:
			s += bin.continuum
		return s

	@property
	def galaxy_varience(self):
		s = np.zeros(len(self.bin[0].noise))
		for bin in self.bin:
			s += bin.noise**2
		return np.sqrt(s)

	@property
	def SNRatio(self):
		return np.array([np.nanmedian(bin.spectrum)/np.nanmedian(bin.noise) 
			for bin in self.bin])

	@property
	def _components(self):
		a = {line:emission_data(self, line, wav) for line, wav in lines.iteritems()}
		a['stellar'] = stellar_data(self)
		return a

	@property
	def e_components(self):
		return list(self.e_line.keys())
	@property
	def list_components(self):
		return list(self.components.keys())

	@property
	def e_line(self):
		return {k:v for k,v in self._components.iteritems() if k!='stellar' and not
			all(v.mask)}

	@property
	def components(self):
		return {k:v for k,v in self._components.iteritems() if k=='stellar' or not
			all(v.mask)}



	def find_common_range(self):
		ranges = np.array([bin.lamLimits for bin in self.bin])
		return np.array([np.max(ranges[:,0]), np.min(ranges[:,1])])

	# Calculate the rest frame for velocities.
	def find_restFrame(self):
		if self.norm_method == "lum":
			self.vel_norm = 0.0
			self.vel_norm = self.components['stellar'].\
				plot['vel'][self.center_bin]
		if self.norm_method == "lwv":
			self.vel_norm = 0.0
			lwv = self.components['stellar'].plot['vel'].unbinned*self.unbinned_flux
			self.vel_norm = np.nanmean(lwv)*np.sum(self.n_spaxels_in_bin)/ \
				np.nansum(self.unbinned_flux)
		if self.norm_method == "sig":
			self.vel_norm = 0.0
			s_sort = sorted(np.unique(self.components['stellar'].plot['sigma']))
			c = np.where(self.components['stellar'].plot['sigma'] > s_sort[-6])
			self.vel_norm = np.mean(D.components['stellar'].plot['velocity'][c])

	def absorption_line(self, absorption_line, uncert=False):
		ab = np.zeros(self.number_of_bins)
		ab_uncert = np.zeros(self.number_of_bins)
		for i, bin in enumerate(self.bin):
			if uncert:
				ab[i], ab_uncert[i] = bin.absorption_line(absorption_line, 
					uncert=uncert)
			else:
				ab[i] = bin.absorption_line(absorption_line, uncert=uncert)
		if uncert:
			return ab, ab_uncert
		else:
			return ab

	

	

	



# taken from http://docs.scipy.org/doc/numpy/user/basics.subclassing.html
# class myArray(np.ndarray):
# 	def __new__(cls, input_array, uncert=None):
# 		# Input array is an already formed ndarray instance
# 		# We first cast to be our class type
# 		obj = np.asarray(input_array).view(cls)
# 		# add the new attribute to the created instance
# 		obj.uncert = uncert
# 		# Finally, we must return the newly created object:
# 		return obj

# 	def __array_finalize__(self, obj):
# 		# see InfoArray.__array_finalize__ for comments
# 		if obj is None: return
# 		self.uncert = getattr(obj, 'uncert', None)
# 	#pass



class _data(object):
# To be used with __init__ method within other classes (stellar_data and 
#	emission_data). NB: these are all 'getter' functions.
# Attributes:
# plot: (dictionary) listing attributes below
# vel (_uncert): (array)
# sigma (_uncert): (array)
# h3 (_uncert): (array)
# h4 (_uncert): (array)
	
	def __init__(self):
		pass
		# self.moment = 0

	@property
	def plot(self):
		if self.moment == 1:
			return {'vel':self.vel}
		elif self.moment == 2:
			return {'vel':self.vel, 'sigma':self.sigma}
		elif self.moment == 3:
			return {'vel':self.vel, 'sigma':self.sigma, 'h3':self.h3}
		elif self.moment == 4:
			return {'vel':self.vel, 'sigma':self.sigma, 'h3':self.h3, 'h4':self.h4}
			
	# # If value is unmeasured.
	# def unset(self, attr):
	# 	for i, bin in enumerate(self.__parent__.bin):
	# 		del bin.components[self.name].__dict__[attr]

	# def setkin(self, attr, value):
	# 	# Ensures only acts on kinematics
	# 	if attr in ['vel','sigma','h3','h4']:
	# 		self.moment = max(self.moment, 
	# 			np.where(np.array(['vel','sigma','h3','h4'])==attr)[0][0]+1)
	# 		for i, bin in enumerate(self.__parent__.bin):
	# 			try:
	# 				sav = bin.components[self.name].__dict__[attr].uncert
	# 				bin.components[self.name].__dict__[attr] = myFloat(value[i])
	# 				bin.components[self.name].__dict__[attr].uncert = sav
	# 			except KeyError:
	# 				pass

	# def setkin_uncert(self, attr, value):
	# 	# Ensures only acts on kinematics
	# 	if attr in ['vel','sigma','h3','h4']:
	# 		for i, bin in enumerate(self.__parent__.bin):
	# 			try:
	# 				bin.components[self.name].__dict__[attr].uncert = value[i]
	# 			except KeyError:
	# 				pass

	# __getattr__ only called when attribute is not found by other means. 
	def __getattr__(self, attr):
		if attr in dynamics:
			m = self.mask_dynamics
			# Normalise to rest frame of stars
			if attr == 'vel':
				kinematics = np.array([bin.components[self.name].__getattr__(attr) 
					- self.__parent__.vel_norm if not m[i] else np.nan 
					for i, bin in enumerate(self.__parent__.bin)])
			else:
				kinematics = np.array([bin.components[self.name].__getattr__(attr) 
					if not m[i] else np.nan 
					for i, bin in enumerate(self.__parent__.bin)])
			# kinematics.uncert = myArray(
			# 	[bin.components[self.name].__dict__[attr].uncert 
			# 	if not m[i] else np.nan 
			# 	for i, bin in enumerate(self.__parent__.bin)])

			# unbinned = np.zeros(self.__parent__.unbinned_flux.shape)
			# uncert_unbinned = np.zeros(self.__parent__.unbinned_flux.shape)

			# for spaxel in range(len(self.__parent__.x)):
			# 	unbinned[self.__parent__.x[spaxel],self.__parent__.y[spaxel]] = \
			# 		kinematics[self.__parent__.bin_num[spaxel]]
			# 	uncert_unbinned[self.__parent__.x[spaxel],self.__parent__.y[spaxel]] = \
			# 		kinematics.uncert[self.__parent__.bin_num[spaxel]]
			# kinematics.unbinned = unbinned
			# kinematics.uncert.unbinned = uncert_unbinned

			# def __getattr__(kinematics, attr):
			# 	if uncert in attr:
			# 		uncert.unbinned = property(lambda self:uncert_unbinned)
			# 		return uncert
			return kinematics
		else:
			return object.__getattribute__(self,attr)


	


class stellar_data(_data):
# Attributes:
# name: (str) 'stellar'
# mask(_dynamics): (boolean array) all set to False
# Inherited attributes from _data object for reading ppxf fitted stellar kinematics
	def __init__(self, parent):
		self.__parent__ = parent
		self.name = 'stellar'
		self.moment = 4
		_data.__init__(self)
	@property
	def mask(self):
		return np.array([False]*self.__parent__.number_of_bins)
	@property
	def mask_dynamics(self):
		return np.array([False]*self.__parent__.number_of_bins)



class emission_data(_data):
# Attributes:
# flux: array of integrated fitted spectrum
# name: emission line name (will also be directory key in Data object)
# wav: emission line quoted wavelength
# equiv_width: array of equivelent width for each bin
# mask: boolean array of if the emission line is masked in bin
# mask_dynamics: boolean array of if all emission lines fixed to the same LOSVD are 
#	masked in bin
# amp_noise: array of Amplitude to Noise ratio for each bin
#
# Inherited attributes from _data object for reading ppxf fitted kinematics
	def __init__(self, parent, name, wav):
		self.__parent__ = parent
		self._name = name
		self._wav = wav
		self.moment = 2
		_data.__init__(self)

	@property	
	def flux(self):
		p = []
		for bin in self.__parent__.bin:
			try:
				p.append(bin.e_line[self._name].flux)
			except KeyError:
				p.append(np.nan)
		return np.array(p)

	@property
	def name(self):
		return self._name

	@property
	def wav(self):
		return self._wav

	@property
	def equiv_width(self):
		return np.array(self.flux/[bin.continuum[np.argmin(np.abs(bin.lam-self.wav))] 
			for bin in self.__parent__.bin])

	@property
	def mask(self):
		p = np.ones(self.__parent__.number_of_bins).astype(bool) # all masked
		# if self.__parent__.gas == 3:
		for bin in self.__parent__.bin:
			try:
				p[bin.bin_number] = bin.e_line[self._name].mask
			except KeyError:
				p[bin.bin_number] = True
		return np.array(p)
		# elif self.__parent__.gas == 2:
		# 	pass
		# elif self.__parent__.gas == 1:
		# 	# Check all emission lines for a detection (i.e. mask = False)
		# 	for bin in self.__parent__.bin:
		# 		for c in bin.e_line.keys():
		# 			if not bin.e_line[c].mask:
		# 				p[bin.bin_number] = bin.e_line[c].mask
		# return p

	@property
	def mask_dynamics(self):
		p = np.ones(self.__parent__.number_of_bins).astype(bool) # all masked
		if self.__parent__.gas == 3:
			p = self.mask
		elif self.__parent__.gas == 2:
			pass
		elif self.__parent__.gas == 1:
			# Check all emission lines for a detection (i.e. mask = False)
			# components = self.__parent__.e_components
			# for bin in self.__parent__.bin:
			# 	for c in components:
			# 		if not bin.e_line[c].mask:
			# 			p[bin.bin_number] = bin.e_line[c].mask
			for k, i in self.__parent__.e_line.iteritems():
				p *= i.mask
		return p

	@property
	def amp_noise(self):
		p = []
		for bin in self.__parent__.bin:
			try:
				p.append(bin.e_line[self._name].AmpNoi)
			except KeyError:
				p.append(np.nan)
		return np.array(p)




class Bin(object):
# Attributes:
# bin_number: (int)
# e_line: (dictionary) of emission_line objects
# temp_weight: (dictionary) of template weights
# lam: (array) x axis of spectrum
# loglam: (array) log of lam
# limLimits: (array) limits of lam
# FWHM_gal: (float) Full-width-half-maximum used for calculating emission lines
# bestfit: (array) bestfit from ppxf
# spectrum: (array) observed spectrum (with cuts from errors2.py)
# noise: (array) observed noise (with cuts from errors2.py)
# contiuum: (array) observed spectrum with emission lines removed
# xspaxels: (array) x-coords of spaxels in bin
# yspaxels: as xspaxels
# xBar: (float) x-coord of center of bin
# yBar: as xBar
# n_spaxels_in_bin: (int) number of spaxels in this bin
# stellar: (_bin_data) object for storing ppxf fitted stellar kinematics
# apweight: (array) weighting of *addative* polynomial used in pPXF fit
# mpweight: (array) weighting of *multiplicative* polynomial used in pPXF fit
# unconvolved_spectrum: (array) The unconvolved spectrum created using the template
#	weights and ap and mp weights.
# unconvolved_lam: (array) wavelength array corresponding to unconvolved_spectrum 
#	(see above)
#
# Methods:
# **DEPRECATED** set_emission_lines (FWHM of observations): requires self._lam is set. 
#	Uses ppxf utilites to determine which lines are within the observed spectrum.
# set_templates (template name, ppxf fitted weighting): Sets the template weights
#	particularly in the emission_line objects. Also now loads the unconvolved 
#	spectrum instead of method below.
# **DEPRECATED** Unconvolved spectrum (wavelength of templates, templates used in 
#	dictionary): 
# 	(array) Returns the unconvolved templates. NB: this are not in the same 
#	wavelength bins as spectrum, bestfit, lam etc, but are binned as the stellar 
#	templates.
# absorption_line (line name): returns absorption line strength (of LICK style 
#	indicies). Use uncert (boolean ) keyword to return uncertainty as well 
#	(defualt False). 

	def __init__(self, bin_number, parent):
		self.__parent__ = parent
		self.bin_number = int(bin_number)



	@property
	def xBar(self):
		return self.__parent__.xBar[self.bin_number]

	@property
	def yBar(self):
		return self.__parent__.yBar[self.bin_number]

	@property
	def xspaxels(self):
		return self.__parent__.x[np.where(self.__parent__.bin_num == bin_number)[0]]

	@property
	def yspaxels(self):
		return self.__parent__.y[np.where(self.__parent__.bin_num == bin_number)[0]]


	@property
	def spectrum(self):
		return np.loadtxt("%s/input/%d.dat" % (self.__parent__.vin_dir_gasMC,
			self.bin_number), unpack=True)

	@property
	def noise(self):
		return np.loadtxt("%s/noise_input/%d.dat" % (self.__parent__.vin_dir_gasMC,
			self.bin_number), unpack=True)

	@property
	def lam(self):
		return np.loadtxt("%s/lambda/%d.dat" % (self.__parent__.vin_dir_gasMC,
			self.bin_number), unpack=True)

	@property
	def bestfit(self):
		return np.loadtxt("%s/bestfit/%d.dat" % (self.__parent__.vin_dir_gasMC,
			self.bin_number), unpack=True)

	@property
	def apweight(self):
		if os.path.isfile("%s/apweight/0.dat" % (self.__parent__.vin_dir_gasMC)):
			return np.loadtxt("%s/apweight/%d.dat" % (self.__parent__.vin_dir_gasMC,
				self.bin_number), unpack=True)
		else:
			return np.array([])

	@property
	def mpweight(self):
		if os.path.isfile("%s/mpweight/0.dat" % (self.__parent__.vin_dir_gasMC)):
			return np.loadtxt("%s/mpweight/%d.dat" % (self.__parent__.vin_dir_gasMC,
				self.bin_number), unpack=True)
		else:
			return np.array([])

	@property
	def n_spaxels_in_bin(self):
		return len(self.xspaxels)

	@property
	def components(self):
		c = {line:emission_line(self, line, wav) for line, wav in lines.iteritems()}
		c['stellar'] = _bin_data(self, 'stellar')
		return c

	@property
	def temp_weight(self):
		temp_name, temp_weight = np.loadtxt("%s/temp_weights/%d.dat" % (
			self.__parent__.vin_dir_gasMC, self.bin_number), unpack=True, 
			dtype='str')
		temp_weight = temp_weight.astype(float)
		return {temp_name[i]:temp_weight for i in xrange(len(temp_name))}

	@property
	def unconvolved_spectrum(self):
		## Solving an error in the Uni system with loading large numbers of files from 
		## the home directory - they have to come from the Data partition instead.
		if cc.getDevice() == 'uni':
			files = glob('%s/Data/idl_libraries/ppxf/MILES_library/' % (cc.base_dir) +
				'm0[0-9][0-9][0-9]V')
		else:
			files = glob("%s/models/miles_library/m0[0-9][0-9][0-9]V" % (cc.home_dir))
		wav = np.loadtxt(files[0], usecols=(0,), unpack=True)
		# ***** NB: spec is binned as the stellar templates are binned, NOT as 
		# self.spectrum, bestfit etc i.e. wav != self.lam  *****
		a = [min(np.where(wav>=self.lam[0])[0]), max(np.where(wav<=self.lam[-1])[0])]
		unc_spec = np.zeros(a[1]-a[0])

		for i,n in enumerate(name):
			self.temp_weight[n] = weight[i]
			if n.isdigit() and weight[i] !=0:
				template =  np.loadtxt(files[int(n)], usecols=(1,), unpack=True)
				unc_spec += template[a[0]:a[1]]*weight[i]
		if len(self.mpweight)!=0:
			unc_spec *= np.polynomial.legendre.legval(np.linspace(-1,1,
				len(unc_spec)), np.append(1, self.mpweight))
		if len(self.apweight)!=0:
			unc_spec += np.polynomial.legendre.legval(np.linspace(-1,1,
				len(unc_spec)), self.apweight)
		return unc_spec

	@property
	def unconvolved_lam(self):
		if cc.getDevice() == 'uni':
			files = glob('%s/Data/idl_libraries/ppxf/MILES_library/' % (cc.base_dir) +
				'm0[0-9][0-9][0-9]V')
		else:
			files = glob("%s/models/miles_library/m0[0-9][0-9][0-9]V" % (cc.home_dir))
		wav = np.loadtxt(files[0], usecols=(0,), unpack=True)
		# ***** NB: spec is binned as the stellar templates are binned, NOT as 
		# self.spectrum, bestfit etc i.e. wav != self.lam  *****
		a = [min(np.where(wav>=self.lam[0])[0]), max(np.where(wav<=self.lam[-1])[0])]
		return wav[a[0]:a[1]]

	@property
	def e_line(self):
		return {k:v for k,v in self.components.iteritems() if k!='stellar'}

	@property
	def stellar(self):
		return self.components['stellar']

	@property
	def loglam(self):
		return np.log(self.lam)
	# @loglam.setter
	# def loglam(self, loglam):
	# 	self._lam = np.exp(loglam)

	@property
	def lamLimits(self):
		lam = self.lam
		return np.array([lam[0], lam[-1]])

	@property
	def continuum(self):
		# NB: Masks not used
		return np.array(self.spectrum - np.nansum([line.spectrum_nomask 
			for key, line in self.e_line.iteritems()],axis=0))


	@property
	def flux(self):
		# NB: only calculated for the common wavelength range.
		a = [min(np.where(self.lam > self.__parent__.common_range[0])[0]),
			max(np.where(self.lam < self.__parent__.common_range[1])[0])]
		return np.trapz(self.spectrum[a[0]:a[1]], x=self.lam[a[0]:a[1]]
			)/self.n_spaxels_in_bin

	def absorption_line(self, absorption_line, uncert=False):
		convolved = self.bestfit - np.nansum([line.spectrum_nomask for key, 
			line in self.e_line.iteritems()], axis=0)
		lam = self.lam/(1+self.components['stellar'].vel/c)
		if uncert:
			return absorption(absorption_line, lam, self.continuum, noise=self.noise,
				unc_lam=self.unconvolved_lam, unc_spec=self.unconvolved_spectrum, 
				conv_spec=convolved)
		else:
			return absorption(absorption_line, lam, self.continuum, 
				unc_lam=self.unconvolved_lam, unc_spec=self.unconvolved_spectrum, 
				conv_spec=convolved)



class _bin_data(object):
# Attributes:
# vel (.uncert): float
# sigma (.uncert): float
# h3 (.uncert): float
# h4 (.uncert): float
	def __init__(self, parent, name):
		self.__parent__ = parent
		self.name = name


	def __getattr__(self, attr):
		if attr in dynamics:
			glamdring_file = '%s/%s.dat' % (self.__parent__.__parent__.vin_dir_gasMC,
				self.__parent__.bin_number)
			t = np.loadtxt(glamdring_file, unpack=True, usecols=(0,), dtype=str)
			i = np.where(t == self.name)[0][0]
			print 'attr ', attr
			print 'dynamics ', dynamics
			j = np.where(dynamics == attr)[0] + 1
			print 'i ', i
			print 'j ', j


			with open(glamdring_file, 'r') as g:
				rows = g.read().splitlines()
				row = np.array(rows[i].split('   '))

			try:
				return float(row[j])
			except IndexError:
				pass
		else:
			return object.__getattribute__(self,attr)



class emission_line(_bin_data):
# Attributes:
# name: str, name of emission line
# _weight: float, weighting from ppxf fit
# wav: wavelength of emission
# flux: float, integrated fitted spectrum.
# spectrum: array if not masked of fitted spectrum
# spectrum_nomask: as spectrum overides masking
# AmpNoi: float, amplitude of fitted spectrum vs noise from Bin object
# mask: boolean: True if emission line is masked in this bin
# Inherited attributes from _bin_data object for storing ppxf fitted kinematics
	def __init__(self, parent, name, wav):
		self.__parent__ = parent
		self.name = str(name)
		self.wav = wav
		self.__threshold__ = 4.0
		if self.__parent__.__parent__.gas == 1:
			t = 'gas'
		elif self.__parent__.__parent__.gas == 2:
			if 'H' in self.name:
				t = 'SF'
			else:
				t = 'Shocks'
		elif self.__parent__.__parent__.gas == 3:
			t = self.name
		_bin_data.__init__(self, self.__parent__, t)

		
	@property
	def flux(self):
		if not self.mask:
			return np.trapz(self.spectrum, x=self.__parent__.lam)
		else:
			return np.nan

	# @flux.setter
	# def flux(self,no_mask=False):
	# 	if not self.mask or no_mask:
	# 		return np.trapz(self.spectrum, x=self.__parent__.lam)
	# 	else:
	# 		return np.nan
	# def spectrum(self, no_mask=False):
	# 	if not self.mask or no_mask:
	# 		return self._spectrum*self.weight
	# 	else:
	# 		return self._spectrum*np.nan

	@property
	def weight(self):
		return self.__parent__.temp_weight[self.name]

	@property
	def _spectrum(self):
		matrix = np.loadtxt("%s/bestfit/matrix/%d.dat" % (
			self.__parent__.__parent__.vin_dir_gasMC, 
			self.__parent__.bin_number),dtype=str)
		return matrix[np.where(matrix[:,0] == self.name)[0], 1:].astype(float)

	@property
	def spectrum(self):
		if not self.mask:
			return self._spectrum*self.weight
		else:
			return self._spectrum*np.nan

	@property
	def spectrum_nomask(self):
		return self._spectrum*self.weight

	
	@property
	def AmpNoi(self):
		return max(self.spectrum_nomask)/self.__parent__.noise[np.argmin(np.abs(
			self.__parent__.lam-self.wav))]

	@property
	def mask(self):
		return self.AmpNoi < self.__threshold__ #or self.AmpNoi > 10**3 
		#and np.isnan(self.vel)
