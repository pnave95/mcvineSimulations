# This program is designed to read a given sim-#.nxs file generated from a single slice (angle) of crystal data and convert from time offset, t_zero, and pixel location maps to tof (microseconds), distance (m), polar angle (rad), azimuthal angle (rad), and weight AND then to vector Q (instrument coordinates) and E and then to H,K,L,E -- and to save those values to a standard csv file

# current status: WORKING

# Note:  Q is given in inverse angstroms, E is given in meV

# Notes:  h5py available by default on SNS cluster, but not with miniconda2, so don't source lj7/.use-miniconda2 before running this program

import numpy as np
import h5py as hf
import sys
import computeEi_fromBeamMonitors
import computeAverageErrorSquared

# compute the theoretical energy for a given value of H,K,L
# this function is only for the specific case of the Mourigal 2D square lattice
def E_Q(H, K, L, Emeasured):
	# define first dispersion
	gamma_k = 0.5*(np.cos(2.0*np.pi*H) + np.cos(2.0*np.pi*K))
	twoTheta = 0.200334842323
	E_1 = 40.0*np.sqrt( (1.0+gamma_k)*(1 - gamma_k*np.cos(twoTheta) ))
	# define second dispersion
	E_2 = 0.0
	# determine which dispersion is vertically (i.e. in energy dimension) closer to the point; assume point belongs to that dispersion
	Etrue = E_1
	if abs(Emeasured - E_1) > abs(Emeasured - E_2):
		Etrue = E_2
	error = Emeasured - Etrue
	error2 = error**2
	Elist = [Etrue, error, error2]
	return Elist
	

def results(ifile, E_i):
	

	# mass of neutron (kg)
	power = 0 - 27
	m = 1.674929*(10**power)
	
	# Reduced planck's constant, also scaled to give inverse angstroms when multiplied by momentum
	power = 0 - 34
	hBar = 1.0545718*(10**power)
	power_scaled = 0 - 24
	hBar_scaled = 1.0545718*(10**power_scaled)	

	# convert from meV to Joules
	power = 0 - 22
	Ei = E_i*1.6021766*(10**power)

	# compute initial velocity and momentum
	v_i = np.sqrt(2.0*Ei / m)
	p_i = m*v_i
	k_i = p_i / hBar_scaled

	# debugger:
	#print "k_i = " + str(k_i)
	#print "v_i = " + str(v_i) + " m/s"
	#print "Ei = " + str(Ei) + " Joules"


	# read in file
	f = hf.File(ifile, "r")

	# debugger + code aid:
	def printname(name):
		print name
	#f.visit(printname)

	# compute time for travel from moderator to sample
	L_mod2sample = 13.60
	t_ms = L_mod2sample / v_i	

	
	# strip angle out of input file name:

	# input is sim-#.nxs; so ifile[-4] should be the period before nxs
	stripped_input = ifile[0:-4]
	# result should now be 4 characters (sim-) plus 3-5 more (#.# with possible -)
	stripped_input = stripped_input[-9:]
	# have now stripped most of any other path components; either / or s is first
	while stripped_input[0] != "s":
		stripped_input = stripped_input[1:]
	# now strip off "sim-"
	angle = stripped_input[4:]
	# debugger:
	#print "Angle = " + angle



	# Now open an output file
	outfile = stripped_input + "-results-v7.h5"  
	out = hf.File(outfile, 'w')  # create new hdf5 file to hold all data for the angle slice in an organized fashion
	# Create groups to store data
	data = out.create_group("data")

	# Now count total number of events
	N_total = 0
	for bank in range(115):
		b = bank + 1  # now starts at 1
		bankName = "bank" + str(b)
		# debugger:
		#print bankName + "\n"
		N = f.get("entry/" + bankName + "_events/total_counts")
		N = np.array(N)
		N = N[0]		# N is now an int # of events for this bank
		N_total += N
	N_total = int(N_total)

	# Now create numpy arrays for the data
	
	# Qarray = Qx,Qy,Qz,|vQ|  for each data point (each row)
	Qarray = np.zeros((N_total,4))
	# HKLarray = H, K, L
	HKLarray = np.zeros( (N_total, 3) )
	# Earray = E, Etrue, Eerror, EerrorSquared
	Earray = np.zeros( (N_total, 4) )
	# CoordinatesArray = polar, azimuthal, d
	CoordinatesArray = np.zeros( (N_total, 3) )
	# Instrument Cartesian Coordinates
	InstrumentCartesian = np.zeros( (N_total, 3) )
	# WeightsArray = weights
	WeightsArray = np.zeros( (N_total, 1) )	

	# Define vectors u,v here for now; later can take them as arguments
	# u is along z-direction (beam direction), and v is a linearly independent vector in the xz plane (but must be less than pi radians rotated from positive z axis otherwise our calculations will result in a left-handed coordinate system where y is vertically down instead of up)

	# Define lattice parameters
	a = 5.0
	b = 5.0
	c = 5.0  
	# define the (cubic) Bravais lattice
	a1 = np.array([a, 0, 0])
	a2 = np.array([0, b, 0])
	a3 = np.array([0, 0, c])
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
	u = np.matrix([1, 0, 2])
	uArray = np.array([1, 0, 2])
	v = np.matrix([1, 0, 0])
	vArray = np.array([1, 0, 0])

	# convert u,v into same space as a1,a2,a3 (crystal Cartesian coordinates)
	U = uArray[0]*b1 + uArray[2]*b3
	V = vArray[0]*b1 + vArray[2]*b3

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
	C = np.matrix([b1, b2, b3])  # this matrix, when multiplied by an "a" basis will give a set of "b" basis vectors in terms of "a" basis vectors;  inverting it should then give "a" basis vectors in terms of "b" basis vectors
	Cinv = np.linalg.inv(C)

	# Now, we want to use Cinv to convert vectors in an a basis into vectors in a b basis
	# so, Q -> CinvMQ = SQ
	S = Cinv*M


	# Perform Data conversion from tof,pixels to vQ,E
	
	Ncounter = 0
	# iterate through each detector bank (1-115)
	for bank in range(115):
		b = bank + 1  # now starts at 1
		bankName = "bank" + str(b)
		# debugger:
		#print bankName + "\n"
		N = f.get("entry/" + bankName + "_events/total_counts")
		N = np.array(N)
		N = N[0]		# N is now an int # of events for this bank
		#print N

		# get number of different tz's (time zero's)
		tzHDF = f.get("entry/" + bankName + "_events/event_time_zero")
		tz = np.array(tzHDF)
		# now rescale from seconds to microseconds
		tz *= 10**(6)
		NumTz = len(tz)
		#print NumTz

		# get time record for each event
		timesHDF = f.get("entry/" + bankName + "_events/event_time_offset")
		times = np.array(timesHDF)

		# get weights for each event
		weightsHDF = f.get("entry/" + bankName + "_events/event_weight")
		weights = np.array(weightsHDF)

		# get pixel_id for each event
		event_pixelHDF = f.get("entry/" + bankName + "_events/event_id")
		event_pixel = np.array(event_pixelHDF)
		
		# pixel_id to pixel tube and subpixel
		pixel_mapHDF = f.get("entry/" + bankName + "_events/pixel_id")
		pixel_map = np.array(pixel_mapHDF)
		lowest_pixel_id = pixel_map[0][0]

		# get pixel coordinates
		# system:  z along beam, y vertical, x horizontal and orth. to beam
		#   azimuthal_angle is angle from x in xy plane (rad)
		#   polar_angle is angle from z (from beam direction) (rad)
		#   distance is distance from sample (spherical coordinates)
		folder = "entry/instrument/" + bankName
		polar_angleHDF = f.get(folder + "/polar_angle")
		polar_angle = np.array(polar_angleHDF)
		azimuthal_angleHDF = f.get(folder + "/azimuthal_angle")
		azimuthal_angle = np.array(azimuthal_angleHDF)
		distanceHDF = f.get(folder + "/distance")
		distance = np.array(distanceHDF)	
		#print distance.shape

		# get starting index for each tz
		tzStartIndexHDF = f.get("entry/" + bankName + "_events/event_index")
		tzStartIndex = np.array(tzStartIndexHDF)
		nextStartIndex = N  #this is default in case there is only 1 value
		#print tzStartIndex
		if NumTz > 1:
			nextStartIndex = tzStartIndex[1]
		start = 1
		for e in range(N):  # iterate through each event for this bank	
		#for e in range(10):	
			if e >= nextStartIndex:
				start += 1
				nextStartIndex = tzStartIndex[start]
			t_zero = tz[start-1]
			#print t_zero
			tof = times[e] - t_zero
			pixel = event_pixel[e]
			shifted_pixel = pixel - lowest_pixel_id
			pixel_tube = int(np.floor( shifted_pixel / 128 ))
			subpixel = int(shifted_pixel - (pixel_tube*128))
			polar = polar_angle[pixel_tube][subpixel]
			azimuth = azimuthal_angle[pixel_tube][subpixel]
			d = distance[pixel_tube][subpixel]
			
			w = weights[e]
			
			# perform computations of kinetic variables
			
			# time from sample to pixel
			t_sp = (tof / (10**6)) - t_ms # do I need to account for emission time?


			# final velocity
			v_f = d / t_sp

			# final energy and scalar momentum
			p_f = m*v_f
			#E_f = (p_f**2)/(2.0*m)
			E_f = 0.5*m*(v_f**2)

			# energy transfer TO the material
			E = Ei - E_f

			# energy transfer in meV
			EmeV = E*6.2415093419*(10**21)

			# get Cartesian unit vector towards pixel
			z = np.cos(polar)
			xy = np.sin(polar)
			x = xy*np.cos(azimuth)
			y = xy*np.sin(azimuth)

			# get vector momentum change
			pf_x = x*p_f
			pf_y = y*p_f
			pf_z = z*p_f - p_i

			# convert to wavevector change
			Qx = pf_x / hBar_scaled
			Qy = pf_y / hBar_scaled
			Qz = pf_z / hBar_scaled

			# get Q as numpy matrix
			Qmat = np.matrix([Qx, Qy, Qz])

			# compute H,K,L with transformation matrix S
			HKL = np.array(S*Qmat.T)

			H = HKL[0]
			K = HKL[1]
			L = HKL[2]

			
			# Compute true energy and errors
			Elist = E_Q(H, K, L, EmeV)
			Etrue = Elist[0]
			error = Elist[1]
			error2 = Elist[2]


			# compute magnitude of Q
			scalarQ = np.sqrt(Qx**2 + Qy**2 + Qz**2 )

			# get instrument Cartesian coordinates
			X = float(d*x)
			Y = float(d*y)
			Z = float(d*z)

			# Now record this information to output file
			Qarray[Ncounter][0] = Qx
			Qarray[Ncounter][1] = Qy
			Qarray[Ncounter][2] = Qz
			Qarray[Ncounter][3] = scalarQ
			CoordinatesArray[Ncounter][0] = polar
			CoordinatesArray[Ncounter][1] = azimuth
			CoordinatesArray[Ncounter][2] = d
			InstrumentCartesian[Ncounter][0] = X
			InstrumentCartesian[Ncounter][1] = Y
			InstrumentCartesian[Ncounter][2] = Z
			Earray[Ncounter][0] = EmeV
			Earray[Ncounter][1] = Etrue
			Earray[Ncounter][2] = error
			Earray[Ncounter][3] = error2
			HKLarray[Ncounter][0] = H
			HKLarray[Ncounter][1] = K
			HKLarray[Ncounter][2] = L
			WeightsArray[Ncounter][0] = w


			Ncounter += 1				
			




	# Close output file (after writing data to file)
	dataQ = out.create_dataset("data/Q", data=Qarray)
	dataQ.attrs['columns'] = "Qx, Qy, Qz, |vQ|"

	dataHKL = out.create_dataset("data/HKL", data=HKLarray)
	dataHKL.attrs['columns'] = "H, K, L"

	dataE = out.create_dataset("data/E", data=Earray)
	dataE.attrs['columns'] = "E (meV), E_theoretical (meV), error, error^2"

	dataW = out.create_dataset("data/weights", data=WeightsArray)
	dataW.attrs['columns'] = "probability of neutron (weight)"

	dataCoord = out.create_dataset("data/Coordinates", data=CoordinatesArray)
	dataCoord.attrs['columns'] = "polar angle (angle from beam direction -- rad) , azimuthal angle (rad),  d (distance to pixel -- meters)"
	
	dataInstrument = out.create_dataset("data/InstrumentCartesianCoordinates", data=InstrumentCartesian)
	dataInstrument.attrs['columns'] = "pixel: x, y, z"
	# how to close the "out" hdf5 results file correctly?

	return outfile	


if __name__ == '__main__':
	# run this if this file is executed as a script 
	
	ifile = sys.argv[1] # path to sim-#.nxs file
	#Ei = float(sys.argv[2]) # incident energy in meV

	# compute Ei (meV)
	beamOutDir = sys.argv[2]
	Ei = computeEi_fromBeamMonitors.computeEi(beamOutDir)

	# perform data reduction and computation
	outfile = results(ifile,Ei)

	# calculate errors
	computeAverageErrorSquared.computeErrors(outfile)


	
