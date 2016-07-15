from matplotlib import pyplot as plt
import numpy as np
import histogram.hdf as hh
import histogram as H


# load file into histogram object
ie = hh.load("out/ienergy.h5")

plt.xlabel("neutron energy (meV)")
plt.ylabel("neutron intensity")
plt.title("10 meV beam profile")
plt.plot(ie.energy, ie.I)
plt.savefig("10meV_beamProfile")
