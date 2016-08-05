'''
The purpose of this script is to compute vQ, |Q|, omega, and E
'''

import numpy as np
import h5py as hf


'''
Description:
	Arguments:
		spherical (2-dim numpy array):  spherical coordinates for each event's recording pixel

'''
def reduce_numpy_array(spherical_arr, tof_arr, initials):

	# determine the number of data points
	N = len(tof_arr)

	# mass of neutron (kg)
	m = 1.674929*(10**-27)

	# Planck's constant (J s)
	hBar = 1.0545718*(10**-34)
	# hBar_scaled is multiplied by 10**10 so that it will give a result in inverse angstroms when Q is divided by it
	hBar_scaled = 1.0545718*(10**-24)



	# identify the initiall energy, momentum, and velocity (SI units)
	Ei = initials[0]
	pi = initials[1]
	vi = initials[2]

	# Q_arr = Qx, Qy, Qz, |Q|, E, omega
	Q_arr = np.zeros( (N,6) )

	# compute time for travel from moderator to sample
	L_mod2sample = 13.60
	t_ms = L_mod2sample / vi

	# iterate through each data point
	for i in range(N):

		# compute tof from sample to pixel
		t_sp = tof_arr[i][0] - t_ms

		# compute final speed
		d = spherical_arr[i][2]
		vf = d / t_sp

		# compute final scalar momentum
		pf = m*vf

		# compute instrument Cartesiann unit vectors
		polar = spherical_arr[i][0]
		azimuthal = spherical_arr[i][1]
		zHat = np.cos(polar)
		xyHat = np.sin(polar)
		xHat = xyHat * np.cos(azimuthal)
		yHat = xyHat * np.sin(azimuthal)	

		# get vector momentum change (SI units)
		delta_px = xHat * pf
		delta_py = yHat * pf
		delta_pz = zHat * pf - pi	

		# convert to wavevector change (units of inverse angstroms)
		Qx = delta_px / hBar_scaled
		Qy = delta_py / hBar_scaled
		Qz = delta_pz / hBar_scaled

		# compute change in energy (Joules)
		Ef = 0.5*m*vf**2
		E = Ei - Ef		# this is backwards because neutron scientists use "energy transfer"

		# compute energy in meV
		E_meV = E*6.2415093419*(10**21)

		# compute omega, in complementary units to Q
		omega = E / hBar_scaled

		# compute scalar Q
		Q = np.sqrt( Qx**2 + Qy**2 + Qz**2 )

		# record results in the numpy array
		Q_arr[i][0] = Qx
		Q_arr[i][1] = Qy
		Q_arr[i][2] = Qz
		Q_arr[i][3] = Q
		Q_arr[i][4] = E_meV
		Q_arr[i][5] = omega


	# return the numpy array of Q,E,omega data
	return Q_arr



def record_QEomega(QEomega, angle, Ei, debug=0):

	outfile = "angle_" + str(angle) + "_.h5"

	# create an output file
	out = hf.File(outfile, 'w')

	# create a group
	dataGroup = out.create_group("data")

	# record numpy array to file
	dataQEomega = out.create_dataset("data/QEomega", data=QEomega)
	dataQEomega.attrs['columns'] = "Qx, Qy, Qz, |vQ|, E (meV), omega"
	dataQEomega.attrs['units'] = "all Q quantities are in inverse angstroms.  omega is in the corresponding 10^-10 inverse seconds.  E is in meV to conform to the usual standard in neutron scattering."

	# record incident energy
	dataEi = out.create_dataset("data/Ei", data=np.array([Ei]))
	dataEi.attrs['columns'] = "incident energy" 

	return outfile


'''
Description:
	Arguments:
		outDir (string):  the path of the "out" directory produced by a MCViNE beam simulation
'''
def reduce_all(recording_dir, outDir, minAngle, maxAngle, incrementAngle, debug=0):

	if debug == 1:
		print ""
		print "(function:  compute_QEomega.reduce_all  )"

	import os
	import glob

	# record current directory
	current_dir = os.getcwd()

	# first change to recording_dir to create a subdirectory to specifically hold the instrument cartesian coordinate data
	os.chdir(recording_dir)

	# handle case of recording_dir with or without trailing slash
	if recording_dir[-1] != "/":
		recording_dir[-1] += "/"

	# append specific subdir
	data_dir = recording_dir + "QEomega-data"

	# make the directory if it does not exist;  change into the directory
	if os.path.isdir(data_dir) == False:
		os.mkdir(data_dir)
		os.chdir(data_dir)
	else:
		os.chdir(data_dir)
		map(os.remove, glob.glob("*h5"))	# remove old files if they exist


	# compute initial energy, momentum, velocity
	import compute_initial_kinetics
	Epv = compute_initial_kinetics.computeInitialKinetics(outDir, debug)
	Ei = compute_initial_kinetics.convert_Joules_to_meV(Epv[0])

	# now iterate through all angles (and thus through all files)
	angle = minAngle
	while angle <= maxAngle:
		filename = "angle_" + str(angle) + "_.h5"
		filepath = recording_dir + "raw-spherical-data/" + filename

		if os.path.isfile(filepath):
			# load file data (spherical coordinates and tof)
			f = hf.File(filepath, 'r')
			spherical = np.array( f.get('data/spherical_coordinates') )
			tof = np.array( f.get('data/tof') )


			# compute Q,E,omega coordinates
			data = reduce_numpy_array(spherical, tof, Epv)

			# record computed data to a new hdf5 file
			newfile = record_QEomega(data, angle, Ei, debug=debug)
		else:
			if debug==1:
				print "Error!  Non-existent file path: " + filepath

		if debug == 1:
			print "angle " + str(angle) + " has been processed"

		angle += incrementAngle


	# return to original directory
	os.chdir(current_dir)

	return 0




if __name__ == "__main__":

	import sys
	import time

	# checkpoint
	start = time.time()

	recording_dir = sys.argv[1]
	outDir = sys.argv[2]

	minAngle = -90.0
	maxAngle = 90.0
	increment = 3.0
	debug = 1

	reduce_all(recording_dir, outDir, minAngle, maxAngle, increment, debug=debug)

	# checkpoint
	end = time.time()

	print "completion time:  " + str(end - start) + " seconds"