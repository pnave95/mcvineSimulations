# the purpose of this program is to determine the average squared error of a simulation

import numpy as np
import h5py as hf



# read in file
ifile = "sim-0.0-results-flat.h5"
f = hf.File(ifile, "r")

EarrayHDF = f.get('E')
WarrayHDF = f.get('weights')
Earray = np.array(EarrayHDF)
Warray = np.array(WarrayHDF)

N = len(Warray)
sumErr2 = 0.0
sumWeights = 0.0
for i in range(N):
	#if Warray[i] > 0.05:
	err2 = Earray[i][3]*Warray[i]
	sumErr2 += err2
	sumWeights += Warray[i]

meanErr2 = sumErr2 / sumWeights
print "mean error squared = " + str(meanErr2)
