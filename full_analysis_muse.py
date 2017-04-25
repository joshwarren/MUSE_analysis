# ==================================================================
# Run full IDL analysis
# ==================================================================
# warrenj 20150918 Routine to call a neccessary wrapper scripts for
# binning, finding best inital guesses, finding templates to use, and
# actually running the pPXF and Gandalf code. 
# By default all routine will analyse NGC3557

import numpy as np
from binning_spaxels_muse import binning_spaxels
from find_template_muse import find_template
import traceback, sys


def full_analysis(galaxy=None, opt='kin'):

	galaxies = np.array([
		'ic1459', 
		'ic4296', 
		'ngc1316',
		'ngc1399'])

	gal=6
	if galaxy is None:
		galaxy = galaxies[gal]
	print galaxy

	targetSN = 60
	set_range = np.array([2000,5500])#7350])

	binning_spaxels(galaxy, targetSN=targetSN, opt=opt, auto_override=True, 
		set_range=set_range) #, debug=True)

	# find_template(galaxy, set_range=set_range)

if __name__=="__main__":
	galaxies = [
		'ic1459',
		'ic4296', 
		'ngc1316',
		'ngc1399'
		]
	# for g in galaxies: full_analysis(galaxy=g, opt='pop')
	try:
		full_analysis(galaxy='ngc1316', opt='kin')
	except Exception as e:
		print 'ngc1316 failed'
		print e
		traceback.print_exc()