# This program is designed to read a given sim-#.nxs file generated from a single slice (angle) of crystal data and convert from time offset, t_zero, and pixel location maps to tof (microseconds), distance (m), polar angle (rad), azimuthal angle (rad), and weight -- and to save those values to a standard csv file

# current status: INCOMPLETE

# Notes:  h5py available by default on SNS cluster, but not with miniconda2, so don't source lj7/.use-miniconda2 before running this program

import numpy as np
import h5py as hf
import sys

def rawInfo(ifile, Ei):
	# get sim-#.nxs file name from first terminal argument
	#ifile = sys.argv[1]

	# read in file
	f = hf.File(ifile, "r")

	# debugger + code aid:
	def printname(name):
		print name
	#f.visit(printname)


	# strip angle out of input file name:

	# input is sim-#.nxs; so ifile[-4] should be the period before nxs
	stripped_input = ifile[0:-4]
	# result should now be 4 characters (sim-) plus 3-4 more (#.# with possible -)
	stripped_input = stripped_input[-8:]
	# have now stripped most of any other path components; either / or s is first
	if stripped_input[0] != "s":
		stripped_input = stripped_input[1:]
	# now strip off "sim-"
	angle = stripped_input[4:]
	# debugger:
	#print "Angle = " + angle



	# Now open an output file
	outfile = stripped_input + ".csv"
	out = open(outfile, 'w+')  # overwrite file if exists
	# debugger:
	#out.write("Success!")
	out.write("Qx.Qy,Qz,E,weight\n")


	# Perform Data conversion from tof,pixels to vQ,E

	# get monitor 1 distance from sample
	m1_distHDF = f.get("entry/monitor1/data")
	m1_distNP = np.array(m1_distHDF)
	m1_dist = np.abs(m1_distNP[0])

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
		tz *= 10**(-6)
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
			if e >= nextStartIndex:
				start += 1
				nextStartIndex = tzStartIndex[start]
			t_zero = tz[start-1]
			#print t_zero
			tof = times[e] - t_zero
			pixel = event_pixel[e]
			shifted_pixel = pixel - lowest_pixel_id
			pixel_tube = np.floor( shifted_pixel / 128 )
			subpixel = shifted_pixel - (pixel_tube*128)
			polar = polar_angle[pixel_tube][subpixel]
			azimuth = azimuthal_angle[pixel_tube][subpixel]
			d = distance[pixel_tube][subpixel]
			
			w = weights[e]
			
			# perform computations of kinetic variables
			v_i = 

			# Now record this information to output file
			#info = str(tof) + "," + str(d) + "," + str(polar) + "," + str(azimuth) + "," + str(w) + "\n"
			out.write(info)
		
			
			




	# Close output file
	out.close()


if __name__ == '__main__':
	# getSingleAngleRaw.py executed as script 
	# do something
	rawInfo(ifile)	
