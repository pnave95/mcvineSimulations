'''
The purpose of this script is to compute H, K, L for each data point
'''

import numpy as np
import h5py as hf


'''
Description:
	Arguments:
		vQ (2-dim numpy array):  Qx, Qy, Qz

		lattice_param_vectors (list):  [a1, a2, a3], where ai is a lattice vector (np array), in the crystal's own local Cartesian coordinate system, which depends on the lattice paramters a, b, c

		uList, vList (3-entry lists):  vectors in 3-dimensional HKL basis.  u is a vector in HKL which is parallel to the instrument z-axis (beam axis), and v is any vector in the xz plane which is linearly independent of u

		angle (float):  sample rotation angle
'''
def reduce_numpy_array(vQ, lattice_param_vectors, uList, vList, angle):

	# define the Bravais lattice
	a1 = lattice_param_vectors[0]
	a2 = lattice_param_vectors[1]
	a3 = lattice_param_vectors[2]

	# compute the reciprocal lattice vectors in terms of the Bravais vectors
	b1 = 2.0*np.pi*np.cross(a2, a3) / np.dot(a1, np.cross(a2, a3))
	b2 = 2.0*np.pi*np.cross(a3, a1) / np.dot(a2, np.cross(a3, a1))
	b3 = 2.0*np.pi*np.cross(a1, a2) / np.dot(a3, np.cross(a1, a2))

	#print "b1 = "
	#print b1
	#print "b2 = "
	#print b2
	#print "b3 = "
	#print b3

	# define u,v in terms of H,K,L (coefficients of reciprocal lattice vectors)
	u = np.matrix(uList)
	uArray = np.array(uList)
	v = np.matrix(vList)
	vArray = np.array(vList)

	# convert u,v into same space as a1,a2,a3 (crystal Cartesian coordinates)
	U = uArray[0]*b1 + uArray[1]*b2 + uArray[2]*b3
	V = vArray[0]*b1 + vArray[1]*b2 + vArray[2]*b3

	# Now U and V are represented in terms of the crystal's Cartesian coordinates
	#print "U = "
	#print U
	#print "V = "
	#print V

	# temporarily, let's ignore the rotation of the angle (e.i. assume that angle = 0.0);  then what we must do next is to obtain unit vectors ex, ey, ez in terms of the crystal's Cartesian coordinate system
	ez_ = U / np.linalg.norm(U)
	ex1_ = V
	ey_ = np.cross(ez_, ex1_); ey_ /= np.linalg.norm(ey_)
	ex_ = np.cross(ey_, ez_)
	#print "ex_ = "
	#print ex_
	#print "ey_ = "
	#print ey_
	#print "ez_ = "
	#print ez_

	# Now, we will try to account for sample rotation;  to do this, we think about the following:
	#   A positive rotation angle corresponds to rotating the crystal from the positive z axis towards the positive x axis;  this is essentially the same as if we had rotated the instrument coordinates in the opposite direction (by the negative of that angle); what we will do is to rotate vectorQ by that negative angle:
	# Q = RQ
	phi = float(angle)*2*np.pi/360.0	# this is sample rotation angle
	phi *= -1.0							# now we take the reverse angle to rotate the instrument coordinates
	# Make rotation matrix
	R1 = [np.cos(phi), 0.0, np.sin(phi)]
	R2 = [0.0, 1.0, 0.0]
	R3 = [-np.sin(phi), 0.0, np.cos(phi)]
	R = np.matrix([R1, R2, R3])
	R = R.T
	#print "R = "
	#print R


	# Now we would convert the Qx,y,z into a form with basis in the crystal's Cartesian coordinate system
	# P is the matrix which will convert Crystal Cartesian coordinates to x,y,z instrument coordinates (not accounting for angle)
	P = np.matrix([ex_, ey_, ez_]); P = P.T
	#print "P = "
	#print P 
	# Now to convert from instrument coordinates to crystal coordinates, we need the inverse of P
	Pinverse = np.linalg.inv(P)
	#print "Pinverse = "
	#print Pinverse

	# So, first we want to rotate Qx,y,z to account for sample rotation (by applying R), then we apply Pinverse to convert those rotated instrument coordinates into crystal coordinates;  so Q = Pinverse*R*Q = MQ
	M = Pinverse*R

	# Q = MQ should now convert Q in terms of x,y,z into Q in terms of Crystal's Cartesian x,y,z coordinates
	# After doing that, we will need express Q as a linear combination of b1,b2,b3;  the coefficients will then be H,K,L
	C = np.matrix([b1, b2, b3]); C = C.T  # this matrix, when multiplied by an "a" basis will give a set of "b" basis vectors in terms of "a" basis vectors;  inverting it should then give "a" basis vectors in terms of "b" basis vectors
	Cinv = np.linalg.inv(C)

	# Now, we want to use Cinv to convert vectors in an a basis into vectors in a b basis
	# so, Q -> CinvMQ = SQ
	S = Cinv*M



	# Now, transpose to get a 3 x N array of data instrument coordinate Q values
	vQ = vQ.T
	# make it into a matrix
	vQ = np.matrix(vQ)
	# Next, we act on this matrix with S;  then convert back to numpy array
	HKL = np.array(S * vQ)
	# Now, transpose back to original N x 3 shape
	HKL = HKL.T

	# return the resulting numpy array
	return HKL


def record_HKL(HKL, angle, debug=0):

	outfile = "angle_" + str(angle) + "_.h5"

	# create an output file
	out = hf.File(outfile, 'w')

	# create a group
	dataGroup = out.create_group("data")

	# record numpy array to file
	dataQEomega = out.create_dataset("data/HKL", data=HKL)
	dataQEomega.attrs['columns'] = "H, K, L"
	dataQEomega.attrs['units'] = "inverse angstroms"

	return outfile



def reduce_all(recording_dir, lattice_param_vectors, uList, vList, minAngle, maxAngle, incrementAngle, debug=0):

	if debug == 1:
		print ""
		print "(function:  compute_HKL.reduce_all  )"

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
	data_dir = recording_dir + "HKL-data"

	# make the directory if it does not exist;  change into the directory
	if os.path.isdir(data_dir) == False:
		os.mkdir(data_dir)
		os.chdir(data_dir)
	else:
		os.chdir(data_dir)
		map(os.remove, glob.glob("*h5"))	# remove old files if they exist


	# now iterate through all angles (and thus through all files)
	angle = minAngle
	while angle <= maxAngle:
		filename = "angle_" + str(angle) + "_.h5"
		filepath = recording_dir + "QEomega-data/" + filename

		if os.path.isfile(filepath):
			# load file data (vector Q)
			f = hf.File(filepath, 'r')
			QEomega = np.array( f.get('data/QEomega') )
			vQ = QEomega[0:,0:3]


			# compute H,K,L
			data = reduce_numpy_array(vQ, lattice_param_vectors, uList, vList, angle)

			# record computed data to a new hdf5 file
			newfile = record_HKL(data, angle, debug=debug)
		else:
			if debug==1:
				print "Error!  Non-existent file path: " + filepath

		if debug == 1:
			print "angle " + str(angle) + " has been processed"

		angle += incrementAngle


	# return to original directory
	os.chdir(current_dir)

	return 0






if __name__ == '__main__':

	import sys

	resultsdir = sys.argv[1]

	a = 5.0
	b = 5.0
	c = 5.0
	a1 = np.array([a, 0, 0])
	a2 = np.array([0, b, 0])
	a3 = np.array([0, 0, c])
	lattice = [a1, a2, a3]

	u = [1, 0, 2]
	v = [1, 0, 0]

	minAngle = -90.0
	maxAngle = 90.0
	increment = 3.0

	debug = 1

	reduce_all(resultsdir, lattice, u, v, minAngle, maxAngle, increment, debug=debug)
