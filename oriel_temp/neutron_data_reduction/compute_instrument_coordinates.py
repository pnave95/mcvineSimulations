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