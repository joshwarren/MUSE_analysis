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
from plot_results_muse import plot_results
from kinematics_muse import kinematics
# from GH_plots import GH_plots
# from plot_absorption import plot_absorption
import matplotlib.pyplot as plt # used for plotting
# from stellar_pop import stellar_pop
# from use_kinemetry import use_kinemetry
# from classify import classify
import traceback, sys

galaxies = [
			'ic1459', 
			'ic4296',
			# 'ngc1316',
			# 'ngc1399',
			]
galaxies = ['ic1459']
# galaxies = ['ic4296']
# galaxies = ['ngc1316']
# galaxies = ['ngc1399']


discard = 0
norm = 'lws' #'lwv'

# Arrays for error catching
gal_err=[]
err = []
trace =[]
for galaxy in galaxies:
	D = None
	print galaxy
	try:
		D = pickler(galaxy, discard=discard, norm=norm, kinemetry=False)
		D = plot_results(galaxy, discard=discard, CO = False, residual="median", 
			norm=norm, D=D, show_bin_num=True)
		# plt.close("all")
		# GH_plots(galaxy, wav_range=wav_range)
		# plt.close("all")
		kinematics(galaxy, discard=discard, D=D)
		# plt.close("all")

		# Requires the IDL kinemetry routine to have been run. 
		# use_kinemetry(galaxy)
		# classify(galaxy)
		
		# D = pickler(galaxy, discard=discard, wav_range=wav_range, norm=norm, opt='pop')
		# D = plot_absorption(galaxy, wav_range=wav_range, vLimit=vLimit, D=D)#, uncert=False)
		# D = stellar_pop(galaxy, wav_range=wav_range, vLimit=vLimit, D=D)
	except Exception as e:
		gal_err.append(galaxy)
		err.append(e)
		trace.append(sys.exc_info())
		# traceback.print_exc()
		 
#v_vd_ellip(wav_range=wav_range)

# Display errors
for i in range(len(gal_err)):
	print ''
	print gal_err[i], ' FAILED'
	print err[i]
	exc_type, exc_value, exc_traceback = trace[i]
	traceback.print_exception(exc_type, exc_value, exc_traceback, file=sys.stdout)
	print ''