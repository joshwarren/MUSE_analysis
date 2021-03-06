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
from sav_for_kinemetry import sav_for_kinemetry
from plot_results_muse import plot_results, mapping
from find_centre import find_centre
from kinematics import kinematics
from rotation_curve_muse import rotation_curve
from plot_absorption_muse import plot_absorption
import matplotlib.pyplot as plt # used for plotting
from plot_stellar_pop_muse import plot_stellar_pop
from use_kinemetry_muse import use_kinemetry
from classify_muse import classify
import traceback, sys
from BPT_muse import BPT
from compare_atlas3d import compare_atlas3d
from fit_disk_binned import fit_disk
from find_limits import find_limits

galaxies = [
			'ic1459', 
			'ic4296',
			'ngc1316',
			'ngc1399'
			]
# galaxies = ['ic1459']
# galaxies = ['ic4296']
# galaxies = ['ngc1316']
# galaxies = ['ngc1399']

m=mapping()
# m.SNR = False
# m.image = False
# m.equivalent_width = False
# m.amp_noise = False
# m.kinematics = False
# m.plot_resid = False
# m.line_ratios = False

discard = 0
norm = 'fit_disk' # 'lws' #'lwv'
MC_dir='' #'_low_res'

# Arrays for error catching
gal_err, err, trace =[], [], []
for galaxy in galaxies:
	D = None
	print galaxy
	try:
		# D = pickler(galaxy, discard=discard, norm=norm, 
		# 	opt='kin'+MC_dir)
		# D = find_centre(galaxy, discard=discard, opt='kin'+MC_dir,
		# 	D=D, instrument='muse')
		# D = sav_for_kinemetry(galaxy, opt='kin'+MC_dir, D=D, 
		# 	instrument='muse')
		# D = plot_results(galaxy, discard=discard, 
		# 	overplot = {'radio':'r', 'CO':'c'}, residual="median",
		# 	norm=norm, D=D, show_bin_num=True, mapping=m, 
		# 	opt='kin'+MC_dir)
		# plt.close("all")
		# find_limits(galaxy, opt='kin'+MC_dir, norm=norm, D=D, 
		# 	instrument='muse')
		# D = rotation_curve(galaxy, D=D, opt='kin'+MC_dir) 
		# BPT(galaxy, D=D, opt='kin'+MC_dir)
		# plt.close("all")
		# fit_disk(galaxy, D=D, opt='kin'+MC_dir)
		# plt.close("all")


		# Requires the IDL kinemetry routine to have been run. 
		# classify(galaxy)
		# use_kinemetry(galaxy)

		# D = kinematics(galaxy, discard=discard, D=D, opt='kin',
		# 	instrument='muse')
		# plt.close("all")

		D = None
		D = pickler(galaxy, discard=discard, norm=norm, 
			opt='pop'+MC_dir)
		# find_limits(galaxy, opt='pop'+MC_dir, norm=norm, D=D, 
		# 	instrument='muse')
		# D = plot_results(galaxy, discard=discard, 
		# 	overplot = {'radio':'r', 'CO':'c'}, residual="median", 
		# 	norm=norm, D=D, show_bin_num=True, mapping=m, 
		# 	opt='pop'+MC_dir)
		# plt.close()
		# D = plot_absorption(galaxy, D=D, opt='pop'+MC_dir, 
		# 	uncert=True, overplot = {'radio':'r', 'CO':'c'})
		# plt.close()
		# D = plot_stellar_pop(galaxy, opt='pop'+MC_dir, 
		# 	method='mostlikely', overplot = {'radio':'r', 'CO':'c'}, 
		# 	D=D, gradient='only')
		# plt.close()
		# D = BPT(galaxy, D=D, opt='pop'+MC_dir, norm=norm)
		# plt.close()
		# D = fit_disk(galaxy, D=D, opt='pop'+MC_dir, instrument='muse')
	except Exception as e:
		gal_err.append(galaxy)
		err.append(e)
		trace.append(sys.exc_info())
		# traceback.print_exc()
		 
# # v_vd_ellip(wav_range=wav_range)
# compare_atlas3d()

# Display errors
for i in range(len(gal_err)):
	print ''
	print gal_err[i], ' FAILED'
	print err[i]
	exc_type, exc_value, exc_traceback = trace[i]
	traceback.print_exception(exc_type, exc_value, exc_traceback, 
		file=sys.stdout)
	print ''
