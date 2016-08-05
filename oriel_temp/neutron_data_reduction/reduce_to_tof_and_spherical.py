'''
This program is designed to read a given sim-#.nxs file generated from a single slice (angle) of crystal data and convert from time offset, t_zero, and pixel location maps to tof (seconds), distance (m), polar angle (rad), azimuthal angle (rad), pixel_id, and weight 
'''

import numpy as np
import h5py as hf


'''
Description:
	Arguments:
		filename (string expected):  the absolute path of a "sim-#.nxs" results file from a MCViNE simulation.
			* may or may not be immediately usable for experimental data; would need to check if experimental data has a "weights" dataset
	Returns:
		file name of the newly created hdf5 file which contains the tof, spherical coordinates, pixel_id, and weight of each event 
'''
def reduce_single_angle(ifile):

	# read in file
	f = hf.File(ifile, "r")

	# debugger + code aid:
	def printname(name):
		print name
	#f.visit(printname)


	# strip angle out of input file name:

	# input is sim-#.nxs; so ifile[-4] should be the period before nxs
	stripped_input = ifile[0:-4]
	# result should now be 4 characters (sim-) plus 3-5 more ( (-)(#)#.# = 3 to 5)
	stripped_input = stripped_input[-9:]
	# have now stripped most of any other path components; either / or s is first
	while stripped_input[0] != "s":
		stripped_input = stripped_input[1:]
	# now strip off "sim-"
	angle = stripped_input[4:]
	# debugger:
	#print "Angle = " + angle


	# Now open an output file
	outfile = "angle_" + angle + "_.h5"  
	out = hf.File(outfile, 'w')  # create new hdf5 file to hold all data for the angle slice in an organized fashion

	# Create group to store data
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

	# record total number of events (summed from all banks)
	Total_events = out.create_dataset("data/total_events", data=np.array([N_total]))

	# record angle
	this_angle = out.create_dataset("data/sample_angle", data=np.array([float(angle)]))


	# Now create numpy arrays for the data

	# tof_array (seconds) = time_offset - time_zero
	tof_array = np.zeros( (N_total,1) )

	# coordinates_array = polar, azimuthal, d
	coordinates_array = np.zeros( (N_total,3) )

	# pixel_id_array = pixel_id
	pixel_id_array = np.zeros( (N_total,1), dtype=int )

	# WeightsArray = weight (probability of neutrons in simulation)
	WeightsArray = np.zeros( (N_total,1) )


	# Perform Data conversion from time offset, time zero, and pixel maps
	
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
			if e >= nextStartIndex:
				start += 1
				nextStartIndex = tzStartIndex[start]
			t_zero = tz[start-1]
			#print t_zero
			tof = times[e] - t_zero
			# convert tof to seconds
			tof /= 10**6
			pixel = event_pixel[e]
			shifted_pixel = pixel - lowest_pixel_id
			pixel_tube = int(np.floor( shifted_pixel / 128 ))
			subpixel = int(shifted_pixel - (pixel_tube*128))
			polar = polar_angle[pixel_tube][subpixel]
			azimuth = azimuthal_angle[pixel_tube][subpixel]
			d = distance[pixel_tube][subpixel]
			
			w = weights[e]

			# Now record this information to numpy arrays
			tof_array[Ncounter][0] = tof

			coordinates_array[Ncounter][0] = polar
			coordinates_array[Ncounter][1] = azimuth
			coordinates_array[Ncounter][2] = d

			pixel_id_array[Ncounter][0] = pixel 

			WeightsArray[Ncounter][0] = w


			Ncounter += 1



	# Close output file (after writing data to file)
	dataTOF = out.create_dataset("data/tof", data=tof_array)
	dataTOF.attrs['columns'] = "time-of-flight (seconds)"

	dataCoord = out.create_dataset("data/spherical_coordinates", data=coordinates_array)
	dataCoord.attrs['columns'] = "polar angle (angle from beam direction -- rad) , azimuthal angle (rad),  d (distance to pixel -- meters)"

	dataPixels = out.create_dataset("data/pixel_id", data=pixel_id_array)
	dataPixels.attrs['columns'] = "pixel_id"

	dataW = out.create_dataset("data/weights", data=WeightsArray)
	dataW.attrs['columns'] = "probability of neutron (weight)"
	
	# how to close the "out" hdf5 results file correctly?

	return outfile


'''
Description:
	Summary:
		This function records the current directory, changes into a specified "recording directory", creates a new subdirectory named "raw-spherical-data", changes into that subdirectory, then systematically processes the simulation data of the range of angles specified within the scattering directory (dirpath), saving the processed data to new hdf5 files, and finally changes back to the original starting directory
	Arguments:
		recording_dir (string):  directory where all reduced data will be stored
		
		minAngle, maxAngle (floats):  minimum and maximum angles that were measured
		
		incrementAngle (float):  angular change between successive measurements
		
		dirpath (string):  the scattering directory where all the "work_#/sim-#.nxs" files can be found
	
	Returns:  0
'''
def reduce_all_angles(recording_dir, minAngle, maxAngle, incrementAngle, dirpath, debug=0):
	#import os.path
	if debug == 1:
		print ""
		print "(function:  reduce_to_tof_and_spherical.reduce_all_angles)"

	import os
	import glob

	# record current directory
	current_dir = os.getcwd()

	# change into recording_dir
	os.chdir(recording_dir)

	# handle case of recording_dir with or without trailing slash
	if recording_dir[-1] != "/":
		recording_dir[-1] += "/"

	# create a subdirectory to hold these preliminary reduced data files (if it already exists, then erase old data files inside);  change into this new subdirectory
	data_dir  = recording_dir + "raw-spherical-data"
	if os.path.isdir(data_dir) == False:
		os.mkdir(data_dir)
		os.chdir(data_dir)
	else:
		os.chdir(data_dir)
		map(os.remove, glob.glob("*h5"))

	# now handle cases of dirpath with or without trailing slash
	if dirpath[-1] != "/":
		dirpath += "/"

	# iterate through all angles which were part of the experiment
	angle = minAngle
	while angle <= maxAngle:
		filename = "sim-" + str(angle) + ".nxs"
		filepath = dirpath + "work_" + str(angle) + "/" + filename

		# implement a check to make sure the file actually exists
		if os.path.isfile(filepath):
			outfile = reduce_single_angle(filepath)
			# debugger:
			print "file: " + filename + " has been processed"

		angle += incrementAngle

	# change back to original directory
	os.chdir(current_dir)

	return 0








if __name__ == '__main__':

	import sys
	import time

	start = time.time()

	scattering_path = sys.argv[1] # path to scattering directory
	results_path = sys.argv[2]


	minAngle = -90.0
	maxAngle = 90.0
	incrementAngle = 3.0

	debug=1

	# perform data reduction and computation
	#outfile = reduce_single_angle(ifile)

	success = reduce_all_angles(results_path, minAngle, maxAngle, incrementAngle, scattering_path, debug)

	end = time.time()

	print "time:  " + str(end - start) + " seconds"

