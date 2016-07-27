# The purpose of this program is to compute the average incident energy of the neutron beam from beam energy-intensity data


import numpy as np
import h5py as hf
import sys

# compute center of intensity just as for center of mass
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

def computeEi(pathToOutDirectory):
	
	# get monitors 1 and 2
	m1_string = pathToOutDirectory + "/mon1-itof-focused.h5"
	m2_string = pathToOutDirectory + "/mon2-itof-focused.h5"

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

	# convert energy from Joules to meV
	E = E * 6.242 * 10**21

	return E


if __name__ == '__main__':

	outDirectory = sys.argv[1]
	Ei = computeEi(outDirectory)

	print "Ei = " + str(Ei)

