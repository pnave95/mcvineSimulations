'''
This script contains first approximations for beam distribution analysis
'''
import numpy as np
import h5py as hf
from matplotlib import pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

'''
Description:
	Arguments:
		params (list):  a list of 4 numpy arrays; params[0],[1] are the bin centers of the two functions to be convolved while params[2],[3] are the probability masses at each bin center
		final_length (int):  length of the numpy arrays (bin centers and probabilities) of the convolved distribution which is returned
	Returns:
		Z:  a 2-element list where the first element is a numpy array of the bin centers and the second is a numpy array of probability intensities (for the convolution of the two probability distributions which were passed as arguments in params)
'''
def convolution_1D(params, final_length, debug=0):
	if debug == 1:
		print ""
		print "(function:  convolution_1D  )"

	bins1 = params[0]
	bins2 = params[1]
	Prob1 = params[2]
	Prob2 = params[3]

	if debug == 1:
		print "np.sum(Prob1) = " + str(np.sum(Prob1))
		print "np.sum(Prob2) = " + str(np.sum(Prob2))

	# determine length of arrays
	L1 = len(bins1)
	L2 = len(bins2)
	# sanity check
	if L1 != len(Prob1):
		print "Error!  should have len(bins1) = len(Prob1)"
	if L2 != len(Prob2):
		print "Error!  should have len(bins2) = len(Prob2)"

	# Z = X + Y, so P(z) = \int_-inf^+inf P_X(X) P_Y(Z-X) dX
	Z_bins = np.zeros(L1 * L2)
	Z_prob = np.zeros(L1 * L2)

	for i in range(L1):
		for j in range(L2):
			Z_bins[i*L1 +j] = bins1[i] + bins2[j]
			Z_prob[i*L1 +j] = Prob1[i] * Prob2[j]

	if debug == 1:
		print "np.sum(Z_prob) = " + str(np.sum(Z_prob))
		print "np.sum(Z_bins) = " + str(np.sum(Z_bins))

	# now sort Z_bins into ascending order
	# Z_length = L1 + L2
	# for i in range(Z_length):
	# 	for j in range(Z_length - i - 1):
	# 		if Z_bins[j+i+1] < Z_bins[i]:
	# 			# then swap so that smaller element is first
	# 			Z_bins_low_temp =
	Z_bins_list = Z_bins.tolist()
	Z_prob_list = Z_prob.tolist() 
	Z = zip(Z_bins_list, Z_prob_list)
	Z_bins_list_sorted = [x for (x,y) in sorted(Z, key=lambda pair: pair[0])]
	Z_prob_list_sorted = [y for (x,y) in sorted(Z, key=lambda pair: pair[0])]

	if debug == 1:
		print "Z length should be len(X) * len(Y) = L1 + L2 = " + str(L1 * L2)
		print "should = len(Z_prob_list_sorted) = len(Z_bins_list_sorted) = " + str(len(Z_prob_list_sorted)) + " = " + str(len(Z_bins_list_sorted))

	# now, combine bins to reduce the total number of bins to the same number of bins as L1
	# TODO:  make this more flexible to take # bins desired as an argument; also, do preprocessing to combine any identical bins that may exist
	reduced_Z_bins = np.zeros(L1)
	reduced_Z_prob = np.zeros(L1)
	for i in range(L1):
		total_prob = 0.0
		com = 0.0			# center of "mass" (probability)
		for j in range(L2):
			index = i*L2 + j
			com += Z_bins_list_sorted[index] * Z_prob_list_sorted[index]
			total_prob += Z_prob_list_sorted[index]
		if total_prob == 0.0:
			com = Z_bins_list_sorted[i*L2 + L2 - 1] - Z_bins_list_sorted[i*L2]
		else:
			com /= total_prob
		reduced_Z_bins[i] = com
		reduced_Z_prob[i] = total_prob

	if debug == 1:
		print "np.sum(reduced_Z_prob) = " + str(np.sum(reduced_Z_prob))
		print "len(reduced_Z_prob) = " + str(len(reduced_Z_prob))
		print "should = len(bins1) = " + str(len(bins1))
		print ""

	Z = [reduced_Z_bins, reduced_Z_prob]

	return Z




'''
Description:
	Arguments:
		pathToOutDirectory (string):  absolute path to the simulation directory "out" which contains the monitor 1 and 2 beam data
		final_length (int):  the length of the final convolved t12 distribution numpy array which is returned (actually 2 arrays are returned)
			* the default = -1 signals to use the length of the first beam monitor arrays
		debug (int):  debug mode
	Returns:
		t12:  a 2-element list.  The first element is a numpy array of bin centers (values of t_12) and the second element is a numpy array of probabilities at each bin center location
'''
def retrieve_t12_distribution(pathToOutDirectory, final_length=-1, debug=0):
	if debug == 1:
		print ""
		print "(function:  incident_beam_distributions.retrieve_t12_distribution  )"

	# get monitors 1 and 2
	if pathToOutDirectory[-1] != "/":
		m1_string = pathToOutDirectory + "/mon1-itof-focused.h5"
		m2_string = pathToOutDirectory + "/mon2-itof-focused.h5"
	else:
		m1_string = pathToOutDirectory + "mon1-itof-focused.h5"
		m2_string = pathToOutDirectory + "mon2-itof-focused.h5"

	m1 = hf.File(m1_string, 'r')
	m2 = hf.File(m2_string, 'r')

	# retrieve bin centers
	m1_bins = np.array(m1.get('I(tof)/grid/tof/bin centers'))
	m2_bins = np.array(m2.get('I(tof)/grid/tof/bin centers'))

	# retrieve intensities
	m1_intensities = np.array(m1.get('I(tof)/data'))
	m2_intensities = np.array(m2.get('I(tof)/data'))

	# if final_length= -1, then use the first beam monitor data length
	if final_length == -1:
		final_length = len(m1_bins)

	# debugger:
	# if debug==1:
	# 	plt.figure("tof m1")
	# 	plt.plot(m1_bins, m1_intensities)
	# 	plt.show()


	# rescale the intensities into the range [0,1]
	m1_max = np.amax(m1_intensities)
	m1_min = np.amin(m1_intensities)
	m1_range = m1_max - m1_min
	m1_intensities = (m1_intensities - m1_min) / m1_range

	m2_max = np.amax(m2_intensities)
	m2_min = np.amin(m2_intensities)
	m2_range = m2_max - m2_min
	m2_intensities = (m2_intensities - m2_min) / m2_range

	# normalize the intensities to sum to 1
	m1_total_intensity = np.sum(m1_intensities)
	m2_total_intensity = np.sum(m2_intensities)
	m1_intensities /= m1_total_intensity
	m2_intensities /= m2_total_intensity

	# convolute probability distributions of -t1 and t2
	prob_minus_t1 = m1_intensities
	bins_minus_t1 = - m1_bins

	prob_t2 = m2_intensities
	bins_t2 = m2_bins

	# if debug==1:
	# 	plt.figure("t2, -t1 (prob distributions)")
	# 	plt.plot(bins_minus_t1, prob_minus_t1)
	# 	plt.plot(bins_t2, prob_t2)
	# 	plt.show()

	# t12 = [bins, probabilities]
	t12 = convolution_1D([bins_minus_t1, bins_t2, prob_minus_t1, prob_t2], final_length, debug)


	# debug
	if debug == 1:
		plt.figure("t21")
		plt.plot(t12[0], t12[1])
		plt.show()

	return t12



def get_rid_of_distribution_zeros(t12, debug=0):
	t12_bins = t12[0]
	t12_probs = t12[1]

	t12_bins_reduced = []
	t12_probs_reduced = []
	for i in range(len(t12_probs)):
		if t12_probs[i] > 0.0 and t12_bins[i] != 0.0:
			t12_bins_reduced.append(t12_bins[i])
			t12_probs_reduced.append(t12_probs[i])

	t12_new = [np.array(t12_bins_reduced), np.array(t12_probs_reduced)]

	# debug
	if debug == 1:
		plt.figure("t21_new")
		plt.plot(t12_new[0], t12_new[1])
		plt.show()

	return t12_new



def get_vi_distribution(t12, debug=0):

	t12_bins = t12[0]
	t12_probs = t12[1]

	# determine length of arrays
	N = len(t12_bins)

	# create new arrays for vi
	d12 = 6.67
	vi_bins = d12 / t12_bins 
	vi_probs = t12_probs 

	# debug
	if debug == 1:
		plt.figure("vi distribution")
		plt.plot(vi_bins, vi_probs)
		plt.show()

	return [vi_bins, vi_probs]



def get_pi_distribution(vi, debug=0):
	# mass of neutron
	m = 1.674929 * 10**-27

	pi_bins = m*vi[0]
	pi_probs = vi[1]
	pi = [pi_bins, pi_probs]

	return pi


def get_Ei_distribution(vi, debug=1):
	# mass of neutron
	m = 1.674929 * 10**-27

	Ei_bins = 0.5*m*(vi[0])**2
	Ei_probs = vi[1]
	Ei = [Ei_bins, Ei_probs]

	# debug
	if debug == 1:
		plt.figure("Ei distribution")
		plt.plot(Ei_bins, Ei_probs)
		plt.show()

	return Ei


def get_vec_pi_distribution(pi, theta, phi, debug=0):
	pi_bins = pi[0]
	pi_probs = pi[1]

	pz_bins = pi_bins * np.cos(theta)
	px_bins = pi_bins * np.sin(theta) * np.cos(phi)
	py_bins = pi_bins * np.sin(theta) * np.sin(phi)

	vec_p = [px_bins, py_bins, pz_bins, pi_probs]

	if debug==1:
		from mpl_toolkits.mplot3d import Axes3D
		fig = plt.figure("distribution in pi_x")
		ax = fig.add_subplot(111, projection='3d')
		ax.scatter(px_bins, py_bins, pz_bins, c=pi_probs)
		plt.show()


	return vec_p


def get_vec_ki_resolution(pi, theta, phi, debug=0):
	pi_bins = np.copy(pi[0])
	pi_probs = np.copy(pi[1])

	# remove bins from both sides until only 50% of total probability is left
	left_removed = 0.0
	left_done = False
	index = 0
	while left_done == False:
		if pi_probs[index] < (0.25 - left_removed):
			left_removed += pi_probs[index]
			pi_probs[index] = 0.0
			index += 1
		else:
			left_done = True 

	right_removed = 0.0
	right_done = False
	index = -1
	while right_done == False:
		if pi_probs[index] < (0.25 - right_removed):
			right_removed += pi_probs[index]
			pi_probs[index] = 0.0
			index -= 1
		else:
			right_done = True

	# now remove the bins with zero probability
	pi_bins_list = []
	pi_probs_list = []
	for i in range(len(pi_bins)):
		if pi_probs[i] > 0.0:
			pi_bins_list.append(pi_bins[i])
			pi_probs_list.append(pi_probs[i])

	pi_bins = np.array(pi_bins_list)
	pi_probs = np.array(pi_probs_list)

	if debug == 1:
		print "np.amin(pi_probs) = " + str(np.amin(pi_probs))


	pz_bins = pi_bins * np.cos(theta)
	px_bins = pi_bins * np.sin(theta) * np.cos(phi)
	py_bins = pi_bins * np.sin(theta) * np.sin(phi)

	#vec_p = [px_bins, py_bins, pz_bins, pi_probs]

	# hBar_scaled is multiplied by 10**10 so that it will give a result in inverse angstroms when Q is divided by it
	hBar_scaled = 1.0545718*(10**-24)

	Qz_bins = pz_bins / hBar_scaled
	Qx_bins = px_bins / hBar_scaled
	Qy_bins = py_bins / hBar_scaled

	Qz_range = np.amax(Qz_bins) - np.amin(Qz_bins)
	Qx_range = np.amax(Qx_bins) - np.amin(Qx_bins)
	Qy_range = np.amax(Qy_bins) - np.amin(Qy_bins)

	if debug==1:
		fig = plt.figure("distribution in vec_pi")
		ax = fig.add_subplot(111, projection='3d')
		ax.scatter(Qx_bins, Qy_bins, Qz_bins, c=pi_probs)
		#plt.show()

		print "ranges:  (Qx, Qy, Qz) = (" + str(Qx_range) + ", " + str(Qy_range) + ", " + str(Qz_range) + ")"

	#ki_res = np.array([[Qx_range, 0.0, 0.0], [0.0, Qy_range, 0.0], [0.0, 0.0, Qz_range]])
	ki_res = np.array([Qx_range, Qy_range, Qz_range])

	# # rotation matrix
	# e11 = [np.cos(phi), -np.sin(phi), 0.0]
	# e21 = [np.sin(phi), np.cos(phi), 0.0]
	# e31 = [0.0, 0.0, 1.0]
	# R1 = np.array([e11, e21, e31]); R1 = R1.T 

	# e12 = [np.sin(theta), np.sin(theta), -np.cos(theta)]
	# e22 = [np.sin(theta), np.sin(theta), -np.cos(theta)]
	# e32 = [np.sin(theta), np.sin(theta), np.cos(theta)]
	# R2 = np.array([e12, e22, e32]); R2 = R2.T 

	# R = np.dot(R1, R2)
	#ki_res = np.dot(R, ki_res)

	return ki_res






if __name__ == '__main__':
	import sys

	outdir = sys.argv[1]
	debug = 1
	if debug == 1:
		from matplotlib import pyplot as plt

	t12 = retrieve_t12_distribution(outdir, debug=0)
	t12 = get_rid_of_distribution_zeros(t12, debug=1)
	vi = get_vi_distribution(t12, debug=0)
	pi = get_pi_distribution(t12, debug=0)
	Ei = get_Ei_distribution(t12, debug=1)

	#plt.figure("pi or Ei")
	#plt.plot(pi[0], pi[1])
	#plt.plot(Ei[0], Ei[1])
	#plt.show()

	theta1 = 1.0
	phi1 = 1.0
	#vec_pi = get_vec_pi_distribution(pi, 0.5, 2.2, debug)
	vec_ki_res = get_vec_ki_resolution(pi, theta1, phi1, debug=0)

	# rotation matrix
	# e11 = [np.cos(phi), -np.sin(phi), 0.0]
	# e21 = [np.sin(phi), np.cos(phi), 0.0]
	# e31 = [0.0, 0.0, 1.0]
	# R1 = np.array([e11, e21, e31]); R1 = R1.T 

	# e12 = [np.sin(theta), np.sin(theta), -np.cos(theta)]
	# e22 = [np.sin(theta), np.sin(theta), -np.cos(theta)]
	# e32 = [np.sin(theta), np.sin(theta), np.cos(theta)]
	# R2 = np.array([e12, e22, e32]); R2 = R2.T 

	# R = np.dot(R1, R2)

	# invert the covariance matrix
	#SigmaInv = np.linalg.inv(vec_ki_res)

	# Now, plot ellipse using the Qx, Qy, Qz parameters
	vec_ki_res = vec_ki_res / np.linalg.norm(vec_ki_res)
	theta_prime = np.linspace(0.0, np.pi, num=100)
	phi_prime = np.linspace(0.0, 2.0*np.pi, num=100)
	# x = np.zeros( 50*50 )
	# y = np.zeros( 50*50 )
	# z = np.zeros( 50*50 )
	# index = 0
	# for theta in theta_prime:
	# 	for phi in phi_prime:
	# 		y[index] = (np.sin(theta)) * (np.sin(phi))
	# 		x[index] = (np.sin(theta)) * (np.cos(phi))
	# 		index += 1
	# x = vec_ki_res[0] * x
	# y = y * vec_ki_res[1]


	# index = 0
	# for theta in theta_prime:
	# 	for phi in phi_prime:
	# 		z[index] = (np.cos(theta)) * vec_ki_res[2]
	# 		index += 1



	theta, phi = np.meshgrid(theta_prime, phi_prime)
	X = np.sin(theta) * np.cos(phi) * vec_ki_res[0]
	Y = np.sin(theta) * np.sin(phi) * vec_ki_res[1]
	Z = np.cos(theta) * vec_ki_res[2] + phi*0.0
	#x, y = np.meshgrid(x, y)
	#z = np.meshgrid(z)
	X = X*np.cos(phi1) - Y*np.sin(phi1)
	Y = X*np.sin(phi1) + Y*np.cos(phi1)

	X = X*np.sin(theta1)
	Y = Y*np.sin(theta1)
	Z = Z*np.cos(theta1) + (X + Y)*np.sin(theta1)


	fig = plt.figure("vQ resolution")
	ax = fig.add_subplot(111, projection='3d')
	ax.plot_surface(X, Y, Z)
	xmin = np.amin(X)
	xmax = np.amax(X)
	ymin = np.amin(Y)
	ymax = np.amax(Y)
	zmin = np.amin(Z)
	zmax = np.amax(Z)
	xyz_max = max(xmax, ymax, zmax)
	# ax.axis([xmin,xmax,ymin,ymax,zmin,zmax])
	#ax.set_xlim3d(xmin,xmax)
	#ax.set_ylim3d(ymin,ymax)
	ax.set_xlabel("dQx/Q")
	ax.set_ylabel("dQy/Q")
	ax.set_zlabel("dQz/Q")
	ax.set_xticks([xmin, 0.0, xmax])
	ax.set_yticks([ymin, 0.0, ymax])
	ax.set_zticks([zmin, 0.0, zmax])
	x_location = np.sin(theta1)*np.cos(phi1)
	y_location = np.sin(theta1)*np.sin(phi1)
	z_location = np.cos(theta1)
	title = "Q-resolution ellipsoid at normalized position \n[" + str(x_location) + ", " + str(y_location) + ", " + str(z_location) + "]"
	ax.set_title(title)
	
	plt.show()




