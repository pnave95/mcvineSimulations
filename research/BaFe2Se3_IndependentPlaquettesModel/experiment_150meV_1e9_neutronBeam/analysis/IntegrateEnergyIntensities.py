# This program will integrate over all Q values and then plot the energy vs intensity of the powder scattering results

from matplotlib import pyplot as plt
import numpy as np
#from histogram.hdf import load
#import histogram as H
import h5py as hf
import sys


#ifile = sys.argv[1]		# read hdf5 file as first argument (iqe.h5)
#iqe = hh.load('/home/patrick/Master_Lenovo/Frustrated Magnetism Research/summer 2016 research/mcvineSimulations_repo/mcvineSimulations/research/BaFe2Se3_IndependentPlaquettesModel/experiment_150MeV_1e9_neutronBeam/scattering/iqe.h5')
#iqe = load('iqe.h5')

filename = sys.argv[1]		# take name of hdf5 file to read
f = hf.File(filename, "r")	# open file in read-only mode

# recursively print all groups and datasets in file
def printname(name):
	print name

#f.visit(printname)

# Try to obtain data set named "iqe/data" from "iqe.h5"
rawData = f.get('iqe/data')
energyBinCenters = f.get('iqe/grid/energy/bin centers')

# read 2-D intensity data into numpy array
#rawData = iqe.I
data = np.nan_to_num(rawData)

# transpose data to put Q on x-axis, energy on y-axis
dataT = np.transpose(data)

# sum along rows (integrate intensity over all Q for given energy)
summed = np.sum(dataT, axis=1)
#summed = np.zeros(first_sum.shape[0])
#for i in range(summed.shape[0]):
#	summed[i] = first_sum[i]
smax = np.amax(summed)
smin = np.amin(summed)
#Irange = np.linspace(smin, smax, num=150)
bins = summed.shape[0]
Erange = np.linspace(0, len(summed), num=len(summed))

# debugger
print "Shape of integrated intensity = " + str(len(summed))

# # Debugger:  plot results
# plt.xlabel("Energy Transfer (meV)")
# plt.ylabel("Intensity")
# plt.plot(energyBinCenters, summed)
# plt.savefig("E_vs_I_from_h5py")


# save array to a csv file
outfile = "integrated_intensity.csv"
points = np.vstack((energyBinCenters, summed))
points = points.T
print "Points Dimension: " + str(points.shape)
np.savetxt(outfile, points, delimiter=',', header='energy,intensity',comments='')
