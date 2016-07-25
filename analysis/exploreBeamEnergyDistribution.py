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
print "min Ei = " + str(E_min)
print "max Ei = " + str(E_max)
print "mean Ei = " + str(mean)

#print(energy)

# chop off all energies below or above certain thresholds
guess = 150.0
low = 130.0
high = 170.0
num = 0
for i in range(N):
	if energy[i] < low or energy[i] > high:
		num += 1

print "num neutrons in range = " + str(num)

E = np.zeros(num)
#print len(E)
count = 0
E_mean = 0.0
for i in range(N):
	if energy[i] < low or energy[i] > high:
		E[count] = energy[i]
		count += 1
		E_mean += energy[i]
#print count
E_mean = E_mean / num
print "E_mean = " + str(E_mean)



# try plotting energy
plt.plot(E)
plt.savefig("E_inRange.png")