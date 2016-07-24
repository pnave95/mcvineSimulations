# The purpose of this program is to estimate the "center of mass" (center of intensity) for the time of flight from beam monitor data

# status: INCOMPLETE 
# NOTE:  not sure what the meaning of these beam monitor data objects is

import h5py as hf
import numpy as np
import sys
from matplotlib import pyplot as plt

filename = sys.argv[1]		# take name of hdf5 file to read
f = hf.File(filename, "r")	# open file in read-only mode

# recursively print all groups and datasets in file
def printname(name):
	print name

#f.visit(printname)

tofHDF = f.get('I(tof)/data')
tof = np.array(tofHDF)

#print len(tof)

# compute mean tof
mean = 0.0
N = len(tof)
mean = sum(tof) / N
tof_min = min(tof)
tof_max = max(tof)
print "min tof = " + str(tof_min)
print "max tof = " + str(tof_max)
print "mean tof = " + str(mean)

# try plotting tof
plt.plot(tof)
plt.savefig("tof.png")