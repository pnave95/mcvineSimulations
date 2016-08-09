'''
The purpose of this script is to propagate intensity distributions from tof and uncertainties to final values like E_f and p_f
'''

import numpy as np
import h5py as hf


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




def multiply_distributions(X, Y, final_length=-1, debug=0):

	if debug == 1:
		print ""
		print "(function:  multiply_distributions  )"

	x_bins = X[0]
	x_probs = X[1]
	y_bins = Y[0]
	y_probs = Y[1]


	if final_length == -1:
		final_length = len(y_bins)

	# Z_probs = np.zeros( (len(x_bins), len(y_bins)) )
	# Z_bins = np.zeros( (len(x_bins), len(y_bins)) )
	# for i in range(len(x_bins)):
	# 	for j in range(len(y_bins)):
	# 		Z[i][j] = x_probs

	# determine length of arrays
	L1 = len(x_bins)
	L2 = len(y_bins)
	# sanity check
	if L1 != len(x_probs):
		print "Error!  should have len(x_bins) = len(x_probs)"
	if L2 != len(y_probs):
		print "Error!  should have len(y_bins) = len(y_probs)"

	# Z = X * Y, so P(z) = P(x) * P(y)
	Z_bins = np.zeros(L1 * L2)
	Z_prob = np.zeros(L1 * L2)

	for i in range(L1):
		for j in range(L2):
			Z_bins[i*L1 +j] = x_bins[i] * y_bins[j]
			Z_prob[i*L1 +j] = x_probs[i] * y_probs[j]

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



'''
Description:
	Arguments:
		var:  variance of the Gaussian distribution for t12
		mu:  expected value of d12 (mean)
		length:  length of the numpy array which will hold the numerical values of the Gaussian
'''
def set_initial_d12_distribution_Gaussian(length, mu=6.67, var=0.01, debug=0):
	sd = np.sqrt(var)
	upperlim = mu + 2.0*sd
	downlim = mu - 2.0*sd
	dist_range = upperlim - downlim

	bins = np.linspace(upperlim, downlim, num=length, endpoint=True)
	probabilities = (1.0 / (sd*np.sqrt(2.0*np.pi)) )*np.exp(- (bins - mu)**2 /(2.0*var) )

	d12 = [bins, probabilities]

	if debug == 1:
		plt.plot(bins, probabilities)
		plt.show()

	return d12



def get_vi_distribution(d12, t12, final_length=-1, debug=0):
	d12_bins = d12[0]
	d12_probs = d12[1]
	t12_bins = t12[0]
	t12_probs = t12[1]
	params = [d12_bins, t12_bins, d12_probs, t12_probs]

	# if final_length= -1, then use length of t12
	if final_length == -1:
		final_length = len(t12_bins)

	# vi = d12 / t12 so need inverse of t12
	inv_t12_bins = 1.0 / t12_bins

	inv_t12 = [inv_t12_bins, t12_probs]

	# multiply distributions
	vi = multiply_distributions(d12, inv_t12, final_length, debug=0)

	if debug == 1:
		plt.figure("vi")
		plt.plot(vi[0], vi[1])
		plt.show()

	return vi




if __name__ == '__main__':
	import sys

	outdir = sys.argv[1]
	debug = 1
	if debug == 1:
		from matplotlib import pyplot as plt

	t12 = retrieve_t12_distribution(outdir, debug=debug)
	d12 = set_initial_d12_distribution_Gaussian(20, debug=0)
	vi = get_vi_distribution(d12, t12, debug=debug)

