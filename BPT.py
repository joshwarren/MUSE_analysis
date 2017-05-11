from checkcomp import checkcomp
cc = checkcomp()
if 'home' not in cc.device:
	import matplotlib # 20160202 JP to stop lack-of X-windows error
	matplotlib.use('Agg') # 20160202 JP to stop lack-of X-windows error
import cPickle as pickle
import matplotlib.pyplot as plt
import numpy as np
from plot_velfield_nointerp import plot_velfield_nointerp




def BPT(galaxy, D=None, opt='kin'):
	print '   BPT'

	analysis_dir = "%s/Data/muse/analysis" % (cc.base_dir)
	# galaxiesFile = "%s/galaxies.txt" % (analysis_dir)

	output = '%s/%s/%s' % (analysis_dir, galaxy, opt)


	if D is None:
		pickleFile = open('%s/pickled/dataObj.pkl' % (output))
		D = pickle.load(pickleFile)
		pickleFile.close()

	for bin in D.bin:
		for e in bin.e_line.itervalues():
			e.__threshold__ = 0
# ------------=============== BPT diagram =================----------
	fig, ax = plt.subplots()
	x = np.log(D.e_line['[OIII]5007d'].flux/D.e_line['Halpha'].flux)
	y = np.log(D.e_line['[SII]6716'].flux/D.e_line['Hbeta'].flux)

	Seyfert2 = (0.72/(x - 0.32) + 1.30 < y) * (1.89 * x + 0.76 < y)
	LINER = ((0.72/(x - 0.32) + 1.30 < y) * (1.89 * x + 0.76 > y)) + (x > 0.32)
	SF = (0.72/(x - 0.32) + 1.30 > y) - (x > 0.32)
	# o = (1.89 * x + 0.76 > y) * (x > 0.32)

	ax.scatter(x[Seyfert2], y[Seyfert2], c='r')
	ax.scatter(x[LINER], y[LINER], c='g')
	ax.scatter(x[SF], y[SF], c='b')
	# ax.scatter(x[o],y[o],c='k')

	x_line1 = np.arange(-2.0, 0.1, 0.001)
	y_line1 = 0.72/(x_line1 - 0.32) + 1.30
	ax.plot(x_line1, y_line1,'k')


	y_line2 = 1.89 * x_line1 + 0.76

	m = y_line2 > y_line1
	a = np.min(x_line1[m])
	x_line2 = np.arange(a, 0.60, 0.001)
	y_line2 = 1.89 * x_line2 + 0.76

	ax.plot(x_line2, y_line2, 'k')

	ax.set_xlim([-1.2, 0.7])
	ax.set_ylim([-1.2, 1.5])

	fig.savefig('%s/plots/BPT.png' % (output))
# ------------================= BPT map ===================----------
	m = np.ones(D.number_of_bins)*np.nan
	m[Seyfert2] = +1
	m[LINER] = 0
	m[SF] = -1

	saveTo = '%s/plots/AGN_map.png' % (output)
	plot_velfield_nointerp(D.x, D.y, D.bin_num, D.xBar, D.yBar, m, vmin=-1.3, vmax=1.3,
		nodots=True, title='Map of AGN', save=saveTo, res=0.2)


	return D


if __name__=='__main__':
	BPT('ic1459')