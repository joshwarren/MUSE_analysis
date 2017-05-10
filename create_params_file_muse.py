import numpy as np
from checkcomp import checkcomp
cc = checkcomp()

opt = 'kin_MC_low_res' 		# kin or pop

galaxies = ['ic1459', 'ic4296', 'ngc1316', 'ngc1399']
gals=[3]

output_file = "%s/MUSE/analysis/params.txt" % (cc.home_dir)
f = open(output_file, 'w')
for gal in gals:
# for gal in range(4):
	galaxy = galaxies[gal]

	tessellation_File = '%s/Data/muse/analysis/%s/%s/setup/voronoi_2d_binning_output.txt' % (
		cc.base_dir, galaxy, opt)

	bin_num = np.loadtxt(tessellation_File, usecols=(2,), skiprows=1, unpack=True)
	n_bins = int(max(bin_num)+1)


	for i in range(n_bins):
		if "kin" in opt:
			f.write("python errors2_muse.py " + str(gal) + " " + str(i) + "\n")
		elif "pop" in opt:
			f.write("python errors3_muse.py " + str(gal) + " " + str(i) + "\n")
f.write("push.sh 'Glamdring MUSE run finished'")

f.close()
print "Done"
