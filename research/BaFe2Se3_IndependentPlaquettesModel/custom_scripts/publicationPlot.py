from matplotlib import pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import numpy as np
import histogram.hdf as hh, histogram as H

#TODO: take input data file as terminal argument
infile = "iqe.h5"

# use histogram.hdf to load hdf5 file into histogram data structure (for the histogram module)
iqe = hh.load(infile)

# get a list containing strings representing the name of each axis
axisNames = iqe.axisNameList()

# retrieve data from histogram object into manipulatable numpy array
rawdata = iqe.I

# convert all "NaN" (Not A Number) elements of 'data' into zeros
data = np.nan_to_num(rawdata)

# transpose the data set to put Q on x-axis and E on y-axis
dataT = np.transpose(data)

# get axis ticks/ranges
Qlist = list(iqe.Q)
Elist = list(iqe.energy)
Qstep = Qlist[1] - Qlist[0]
Estep = Elist[1] - Elist[0]
Qmin = Qlist[0] - 0.5*Qstep
Qmax = Qlist[-1] + 0.5*Qstep
Emin = Elist[0] - 0.5*Estep
Emax = Elist[1] + 0.5*Estep

# set maximum cutoff intensity (this makes lower intensities seem "brighter")
maxIntensity = 0.02

# set up plot
plt.figure(1)
plt.xlabel(axisNames[0])  # plot Q along x-axis
plt.ylabel(axisNames[1])  # plot energy along y-axis
plt.title("I(Q,E)")


# Qlist and Elist contain 'bin centers' for the x and y axes, respectively, to properly scale the plot
plt.pcolormesh(Qlist, Elist, dataT, vmax=maxIntensity)
plt.colorbar()
plt.plot()

# set x and y (Q and E) limits of what to actually display
Xmin = 0.0
Xmax = Qlist[-1]
Ymin = -20.0
Ymax = Elist[-1]
plt.xlim( (Xmin, Xmax) )
plt.ylim( (Ymin, Ymax) )

plt.savefig('I_Q,E_pmeshcolor')
