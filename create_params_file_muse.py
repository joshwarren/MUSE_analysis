import numpy as np
from checkcomp import checkcomp
cc = checkcomp()

opt = 'pop_no_Na' 	# must contain kin or pop

galaxies = ['ic1459', 'ic4296', 'ngc1316', 'ngc1399']
gals=[2]

output_file = "%s/MUSE_project/analysis/params_muse.txt" % (cc.home_dir)
f = open(output_file, 'w')
for gal in gals:
	galaxy = galaxies[gal]

	tessellation_File = '%s/Data/muse/analysis/%s/%s/' % (cc.base_dir, galaxy, opt) + \
		'setup/voronoi_2d_binning_output.txt' 

	bin_num = np.loadtxt(tessellation_File, usecols=(2,), skiprows=1, unpack=True)
	n_bins = int(max(bin_num)+1)

	for i in range(n_bins):
		f.write("python errors2_muse.py " + str(gal) + " " + opt + " " + str(i) + 
			"\n")

	f.write("push.sh 'Glamdring MUSE finished %s'\n"%(galaxy))

if 'pop' in opt:
	for gal in gals:
		galaxy = galaxies[gal]

		tessellation_File = '%s/Data/muse/analysis/%s/%s/' % (cc.base_dir, galaxy, opt) + \
			'setup/voronoi_2d_binning_output.txt' 

		bin_num = np.loadtxt(tessellation_File, usecols=(2,), skiprows=1, unpack=True)
		n_bins = int(max(bin_num)+1)

		for i in range(n_bins):
			# if 'kin' in opt:
			f.write("python pop.py muse " + str(gal) + " " + opt + " " + str(i) + 
				"\n")
f.write("push.sh 'Glamdring MUSE run finished'")

f.close()
print "Done"
