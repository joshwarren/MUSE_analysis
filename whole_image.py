# Routine to create a 'bin' which contains the entire image.
import numpy as np
import os
from astropy.io import fits
from errors2_muse import errors2, get_dataCubeDirectory
from checkcomp import checkcomp
cc = checkcomp()

def whole_image(galaxy):
	galaxies = np.array(['ic1459', 'ic4296', 'ngc1316', 'ngc1399'])
	i_gal = np.where(galaxies == galaxy)[0][0]

	f = fits.open(get_dataCubeDirectory(galaxy))
	s = f[1].data.shape

	x = np.tile(np.arange(s[1]),s[2])
	y = np.arange(s[1]).repeat(s[2])
	binNum = np.zeros(s[1]*s[2])
	
	xBar = np.array([s[1]/2])
	yBar = np.array([s[2]/2])
	xBin = np.ones(s[1] * s[2]) * xBar[0]
	yBin = np.ones(s[1] * s[2]) * yBar[0]

	temp = "{0:5}{1:5}{2:8}{3:10}{4:10}\n"
	temp2 = "{0:12}{1:12}\n"


	saveTo = "%s/Data/muse/analysis" % (cc.base_dir)
	with open("%s/galaxies.txt" % (saveTo), 'r') as f:
		galaxies = f.readlines()
	if "SN_kin_whole_gal" not in galaxies[0]:
		with open("%s/galaxies.txt" % (saveTo), 'w') as f:
			f.write(galaxies[0]+"   SN_kin_whole_gal \n")
			for i in xrange(1, len(galaxies)):
				f.write(galaxies[i] + "   1000 \n")



	saveTo = "%s/Data/muse/analysis/%s/kin_whole_gal/setup" % (cc.base_dir, galaxy)

	if not os.path.exists(saveTo):
		os.makedirs(saveTo)

	with open("%s/voronoi_2d_binning_output.txt" % (saveTo), 'w') as f:
		f.write(temp.format('X"', 'Y"', 'BIN_NUM', 'XBIN', 'YBIN'))
		for i in range(len(xBin)):
			f.write(temp.format(str(int(x[i])), str(int(y[i])), str(int(binNum[i])), 
				str(round(xBin[i],5)), str(round(yBin[i],5))))


	with open("%s/voronoi_2d_binning_output2.txt" % (saveTo), 'w') as f:
		f.write(temp2.format('XBAR','YBAR'))
		for i in range(len(xBar)):
			f.write(temp2.format(str(round(xBar[i],5)), str(round(yBar[i],5))))

	errors2(i_gal, 'kin_whole_gal', 0)



if __name__=='__main__':
	for g in ['ic1459', 'ic4296', 'ngc1316', 'ngc1399']:
		whole_image(g)