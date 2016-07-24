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

energyHDF = f.get('ienergy/data')
energy = np.array(energyHDF)

print "len(energy) = " + str(len(energy))

# compute mean tof
mean = 0.0
N = len(energy)
mean = sum(energy) / N
E_min = min(energy)
E_max = max(energy)
print "min tof = " + str(E_min)
print "max tof = " + str(E_max)
print "mean tof = " + str(mean)

# try plotting tof
#plt.plot(tof)
#plt.savefig("tof.png")