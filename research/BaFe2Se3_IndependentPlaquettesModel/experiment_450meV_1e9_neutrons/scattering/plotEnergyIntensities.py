# This program will integrate over all Q values and then plot the energy vs intensity of the powder scattering results

from matplotlib import pyplot as plt
import numpy as np
import histogram.hdf as hh
import histogram as H


ifile = "iqe.h5"
iqe = hh.load(ifile)

# read 2-D intensity data into numpy array
rawData = iqe.I
data = np.nan_to_num(rawData)

# transpose data to put Q on x-axis, energy on y-axis
dataT = np.transpose(data)

# sum along rows (integrate intensity over all Q for given energy)
summed = np.sum(dataT, axis=1)

# plot
plt.xlabel("Energy Transfer (meV)")
plt.ylabel("Intensity")
plt.plot(iqe.energy, summed)
plt.savefig("scatteredneutron_energy_vs_intensity")

