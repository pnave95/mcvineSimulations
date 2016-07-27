

import numpy as np
import h5py as hf
from matplotlib import pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import sys

# take input file
ifile = sys.argv[1]

# read hdf input file
f = hf.File(ifile, 'r')

# get E, H, K
E_hdf = f.get('E')
HKL_hdf = f.get('HKL')

# convert to numpy arrays
E_all = np.array(E_hdf)
HKL = np.array(HKL_hdf)
E = E_all[10000:15000,0]
E = np.array(E)
print "Shape of E = " + str(E.shape)

print "HKL.shape = " + str(HKL.shape)
H = HKL[10000:15000,0]
H = np.array(H)
print "shape of H = " + str(H.shape)
K = HKL[10000:15000,1]
K = np.array(K)
print "shape of K = " + str(K.shape)


fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
#ax = Axes3D(plt.gcf())
ax.scatter(H, K, E)
plt.show()







