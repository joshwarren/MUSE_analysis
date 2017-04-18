import numpy as np
from checkcomp import checkcomp
cc = checkcomp()

opt = 'kin' 		# kin or abs

galaxies = ['ic1459', 'ic4296', 'ngc1316', 'ngc1399']
gals=[1]

output_file = "params.txt"
f = open(output_file, 'w')
for gal in gals:
# for gal in range(4):
	galaxy = galaxies[gal]

	tessellation_File = '%s/Data/muse/analysis/%s/voronoi_2d_binning_output_%s.txt' % (
		cc.base_dir, galaxy, opt)

	bin_num = np.loadtxt(tessellation_File, usecols=(2,), skiprows=1, unpack=True)
	n_bins = int(max(bin_num)+1)


	for i in range(n_bins):
		if opt == "kin":
			f.write("python errors2_muse.py " + str(gal) + " " + str(i) + "\n")
		elif opt == "abs":
			f.write("python errors3_muse.py " + str(gal) + " " + str(i) + "\n")
f.write("push 'Glamdring MUSE run finished'")

f.close()
print "Done"
