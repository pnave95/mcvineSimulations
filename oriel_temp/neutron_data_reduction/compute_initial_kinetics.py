'''
The purpose of this script is to compute the average incident energy of the neutron beam from beam energy-intensity data.  It uses the /out/mon1-itof-focused.h5 and mon2-itof-focused.h5 files.  These files contain intensity values (proportional to number of neutrons) for a finite number of time bins.  That is, the range of tof of the neutrons in the incident beam is partitioned and the intensity associated to each of these bins/partitions is proportional to the number of neutrons with tof (as measured up to the monitor position) which fall into that bin (in MCViNE, the probability weight of the neutron is probably also used)

It would probably be beneficial in the future to use the data before it is binned, if possible.
'''


import numpy as np
import h5py as hf
import sys

# compute center of intensity just as for center of mass (i.e. a weighted sum)
def centerOfMassIntensity(monitor):
	com = 0.0	# com = center of mass
	m = 0.0		# m = total mass
	intensity = np.array(monitor.get('I(tof)/data'))
	binCenters = np.array(monitor.get('I(tof)/grid/tof/bin centers'))
	numBins = len(binCenters)

	for i in range(numBins):
		com += intensity[i]*binCenters[i]
		m += intensity[i]

	com /= m
	return com

#  This function expects a path to the 'out' directory where the monitor energy-intensity hdf5 files are stored (or other directory if the files have been moved)
def computeInitialKinetics(pathToOutDirectory):
	
	# get monitors 1 and 2
	if pathToOutDirectory[-1] != "/":
		m1_string = pathToOutDirectory + "/mon1-itof-focused.h5"
		m2_string = pathToOutDirectory + "/mon2-itof-focused.h5"
	else:
		m1_string = pathToOutDirectory + "mon1-itof-focused.h5"
		m2_string = pathToOutDirectory + "mon2-itof-focused.h5"

	m1 = hf.File(m1_string, 'r')
	m2 = hf.File(m2_string, 'r')

	# compute "center of mass" for each beam monitor (actually center of intensity along time axis)
	com1 = centerOfMassIntensity(m1)
	com2 = centerOfMassIntensity(m2)

	# compute time, distance, and velocity between monitors
	t = com2 - com1
	d = 6.67 # meters
	v = d / t

	# mass of neutron
	m = 1.674929 * 10**-27

	# energy
	E = 0.5 * m * v**2

	# momentum
	p = m*v

	# convert energy from Joules to meV
	#E = E * 6.242 * 10**21

	# make list to hold initial energy, momentum, and velocity (in SI units)
	initial = [E, p, v]

	return initial


def convert_Joules_to_meV(E):
	E *= 6.242 * 10**21
	return E



if __name__ == '__main__':

	outDirectory = sys.argv[1]
	Ei, pi, vi = computeInitialKinetics(outDirectory)
	E = convert_Joules_to_meV(Ei)

	print "Ei = " + str(E) + " meV"

