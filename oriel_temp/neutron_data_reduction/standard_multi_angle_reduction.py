'''
The purpose of this script is to perform all standard data reduction for all specified angles:
	- tof, spherical pixel coordinates, weight
	- pixel cartesian instrument coordinates
	- Ei
	- Qx, Qy, Qz, |vQ|, E, omega
	- H, K, L
'''
import numpy as np




'''
Description:
	Arguments:
		scattering_dir (string):  directory holding all "work_#/sim-#.nxs" simulation results

		out_dir (string):  "out" beam simulation directory

		minAngle, maxAngle (float):  min and max sample rotation angles measured

		incrementAngle (float):  sample angle difference between sample measurements

		orientation (list of sublists):  [u, v] where u and v are both 3-element lists representing the initial orientation of the sample

		lattice_param_vectors (list of np arrays):  holds Bravais lattice vectors

		debug:  debug mode

	Returns:
'''
def standard_reduction(scattering_dir, out_dir, results_dir, minAngle, maxAngle, incrementAngle, orientation, lattice_param_vectors, debug=0):

	if debug == 1:
		print ""
		print "(function:  standard_multi_angle_reduction.standard_reduction  )"

	# import reduction scripts
	import reduce_to_tof_and_spherical
	import compute_instrument_coordinates
	import compute_initial_kinetics
	import compute_QEomega
	import compute_HKL

	# get basic tof and pixel coordinate / id information
	reduce1 = reduce_to_tof_and_spherical.reduce_all_angles(results_dir, minAngle, maxAngle, incrementAngle, scattering_dir, debug)

	# compute Cartesian coordinates of events' pixels
	reduce2 = compute_instrument_coordinates.record_instrument_cartesian_all_angles(results_dir, minAngle, maxAngle, incrementAngle, debug)

	# compute Q,E,omega
	reduce3 = compute_QEomega.reduce_all(results_dir, out_dir, minAngle, maxAngle, incrementAngle, debug)

	# compute HKL
	u = orientation[0]
	v = orientation[1]
	reduce4 = compute_HKL.reduce_all(results_dir, lattice_param_vectors, u, v, minAngle, maxAngle, incrementAngle, debug)

	return 0




if __name__ == "__main__":

	import sys
	import time

	# checkpoint
	start = time.time()

	scattering_dir = sys.argv[1]
	outDir = sys.argv[2]
	recording_dir = sys.argv[3]

	minAngle = -90.0
	maxAngle = 90.0
	increment = 3.0
	debug = 1

	a = 5.0
	b = 5.0
	c = 5.0
	a1 = np.array([a, 0, 0])
	a2 = np.array([0, b, 0])
	a3 = np.array([0, 0, c])
	lattice_param_vectors = [a1, a2, a3]

	u = [1, 0, 2]
	v = [1, 0, 0]

	orientation = [u, v]

	# test
	success = standard_reduction(scattering_dir, outDir, recording_dir, minAngle, maxAngle, increment, orientation, lattice_param_vectors, debug)

	# checkpoint
	end = time.time()

	print "completion time:  " + str(end - start) + " seconds"