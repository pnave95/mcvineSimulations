'''
The purpose of this script is to compute the the pixel location of each event in Instrument Cartesian coordinates
'''

# Status:  UNTESTED

import numpy as np
import h5py as hf


'''
Description:
	Arguments:
		spherical (numpy array):  a 2-dimensional numpy array which holds the spherical coordinates of the pixel for each recorded event
		debug (int):  debug mode (defaults to no debug mode)
	Returns:
		cartesian (numpy array):  a 2-dimensional numpy array holding the x,y,z (instrument Cartesian coordinates) of the pixel associated to each detected event listed in 'spherical'
'''
def get_instrument_cartesian_from_array(spherical, debug=0):


	# sanity check
	if spherical.shape[1] != 3:
		print "Error!  Array of spherical coordinates was expected to have 3 columns, but instead has " + str(spherical.shape[1])

	# get number of data points
	N = spherical.shape[0]

	# create new array to hold Cartesian coordinates
	cartesian = np.zeros( (N, 3) )

	# iterate through all data points, computing x,y,z for each
	for i in range(N):
		polar = spherical[i][0]
		azimuthal = spherical[i][1]
		r = spherical[i][2]

		z = r*np.cos(polar)
		xy = r*np.sin(polar)
		x = xy*np.cos(azimuthal)
		y = xy*np.sin(azimuthal)

		cartesian[i][0] = x
		cartesian[i][1] = y
		cartesian[i][2] = z

	# return the Cartesian results (in a numpy array)
	return cartesian


'''
Description:
	Arguments:
		coordinates (2-dim numpy array):  N x 3 numpy array holding x,y,z values (instrument Cartesian coordinates) of the pixel for each event
		angle (float):  the sample angle which this data was recorded at
		debug (int):  debug mode code
	Returns:
		outfile (string):  the name of the hdf5 file to which the data was recorded
'''
def record_instrument_cartesian_for_single_angle(coordinates, angle, debug=0):

	outfile = "angle_" + str(angle) + "_.h5"

	# create an output file
	out = hf.File(outfile, 'w')

	# create a group
	dataGroup = out.create_group("data")

	# record numpy array to file
	dataCoord = out.create_dataset("data/instrument_cartesian", data=coordinates)
	dataCoord.attrs['columns'] = "x, y, z (instrument Cartesian coordinates of pixel"

	return outfile






def record_instrument_cartesian_all_angles(recording_dir, minAngle, maxAngle, incrementAngle, debug=0):
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
	cartesian_dir = recording_dir + "cartesian-data"

	# make the directory if it does not exist;  change into the directory
	if os.path.isdir(cartesian_dir) == False:
		os.mkdir(cartesian_dir)
		os.chdir(cartesian_dir)
	else:
		os.chdir(cartesian_dir)
		map(os.remove, glob.glob("*h5"))	# remove old files if they exist


	# now iterate through all angles (and thus through all files)
	angle = minAngle
	while angle <= maxAngle:
		filename = "angle_" + str(angle) + "_.h5"
		filepath = recording_dir + "raw-spherical-data/" + filename

		if os.path.isfile(filepath):
			# load file data (spherical coordinates)
			f = hf.File(filepath, 'r')
			spherical = np.array( f.get('data/spherical_coordinates') )

			# compute cartesian coordinates
			cartesian = get_instrument_cartesian_from_array(spherical)

			# record cartesian coordinates to a new hdf5 file
			newfile = record_instrument_cartesian_for_single_angle(cartesian, angle)
		else:
			if debug==1:
				print "Error!  Non-existent file path: " + filepath

		if debug == 1:
			print "angle " + str(angle) + " has been processed"

		angle += incrementAngle


	# return to original directory
	os.chdir(current_dir)




if __name__ == '__main__':
	import sys
	import time

	# starting checkpoint
	start = time.time()

	recording_dir = sys.argv[1]

	minAngle = -90.0
	maxAngle = 90.0
	increment = 3.0

	record_instrument_cartesian_all_angles(recording_dir, minAngle, maxAngle, increment, debug=1)

	end = time.time()

	print "completion time:  " + str(end - start) + " seconds"