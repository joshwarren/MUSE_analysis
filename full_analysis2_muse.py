## ==================================================================
## Run full python analysis
## ==================================================================
## warrenj 20150923 Routine to call a neccessary wrapper scripts for 
## plotting results, and calculating the alignment angles. 

from checkcomp import checkcomp
cc = checkcomp()
if 'home' not in cc.device:
	import matplotlib # 20160202 JP to stop lack-of X-windows error
	matplotlib.use('Agg') # 20160202 JP to stop lack-of X-windows error
from pickler_muse import pickler
from sav_for_kinemetry_muse import sav_for_kinemetry
from plot_results_muse import plot_results, mapping
from kinematics_muse import kinematics
from rotation_curve_muse import rotation_curve
from plot_absorption_muse import plot_absorption
import matplotlib.pyplot as plt # used for plotting
from plot_stellar_pop_muse import plot_stellar_pop
from use_kinemetry_muse import use_kinemetry
from classify_muse import classify
import traceback, sys
from BPT import BPT
from compare_atlas3d import compare_atlas3d

galaxies = [
			'ic1459', 
			'ic4296',
			'ngc1316',
			'ngc1399'
			]
# galaxies = ['ic1459']
# galaxies = ['ic4296']
# galaxies = ['ngc1316']
galaxies = ['ngc1399']

m=mapping()
# m.SNR = False
# m.image = False
# m.equivalent_width = False
# m.amp_noise = False
# m.kinematics = False
# m.plot_resid = False
# m.line_ratios = False

discard = 0
norm = 'lws' #'lwv'
MC_dir=''#'_low_res'

# Arrays for error catching
gal_err, err, trace =[], [], []
for galaxy in galaxies:
	D = None
	print galaxy
	try:
		# D = pickler(galaxy, discard=discard, norm=norm, opt='kin'+MC_dir)
		# D = sav_for_kinemetry(galaxy, opt='kin'+MC_dir, D=D)		
		# D = plot_results(galaxy, discard=discard, overplot = {'radio':'r', 'xray':'c'}, 
		# 	residual="median", norm=norm, D=D, show_bin_num=True, mapping=m, 
		# 	opt='kin'+MC_dir)
		# plt.close("all")
		# D = kinematics(galaxy, discard=discard, D=D, opt='kin') # Only run 
		# 														# for opt='kin'
		# D = rotation_curve(galaxy, D=D, opt='kin'+MC_dir) 
		# BPT(galaxy, D=D, opt='kin'+MC_dir)
		# plt.close("all")

		# Requires the IDL kinemetry routine to have been run. 
		classify(galaxy)
		use_kinemetry(galaxy)

		D = None
		# D = pickler(galaxy, discard=discard, norm=norm, opt='pop'+MC_dir)
		# D = plot_absorption(galaxy, D=D, opt='pop'+MC_dir, uncert=False)
		# plot_stellar_pop(galaxy, opt='pop'+MC_dir, D=D)
	except Exception as e:
		gal_err.append(galaxy)
		err.append(e)
		trace.append(sys.exc_info())
		# traceback.print_exc()
		 
# # v_vd_ellip(wav_range=wav_range)
compare_atlas3d()

# Display errors
for i in range(len(gal_err)):
	print ''
	print gal_err[i], ' FAILED'
	print err[i]
	exc_type, exc_value, exc_traceback = trace[i]
	traceback.print_exception(exc_type, exc_value, exc_traceback, file=sys.stdout)
	print ''
