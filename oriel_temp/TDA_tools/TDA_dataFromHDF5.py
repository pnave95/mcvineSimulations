'''
This script contains functions for retrieving data from hdf5 files and converting that data to numpy arrays
'''

import numpy as np
import h5py as hf


'''
Description:
	Inputs:
		filename:  the file name of the hdf5 file containing the metric data (i.e. the coordinates) for the data points
		dataLocations:  a list containing one or more strings; each string represents the hdf5 absolute path to a dataset which contains at least one column of metric data
		dataColumns:  a list of sublists; the number of sublists should equal the number of strings in the 'dataLocations' list; each sublist contains nonnegative integers corresponding to the column indices of the metric data for the corresponding data set in 'dataLocations'
		minRow:  the first row of data (i.e. the index, starting at zero, of the first data point to be considered)
		maxRow:  the last row of data to be considered
		debugMode:  debug mode
	Returns:
		metricData:  a numpy array which contains all the metric data
'''
def get_data_from_single_hdf5(filename, dataLocations, dataColumns, minRow=0, maxRow=-1, debugMode=0):
	
	# read in hdf5 file
	f = hf.File(filename, 'r')

	# get all data sets as numpy arrays contained in a list
	dataSets = []
	for dset in dataLocations:
		new = np.array( f.get(dset) )
		dataSets.append(new)

	# now, get specific columns as numpy arrays and add them to a new list
	data = []
	for i in range(len(dataLocations)):
		for j in dataColumns[i]:
			data.append( dataSets[i][minRow:maxRow, j] )
			#print dataSets[i].shape


	data = np.transpose(np.array(data))


	if debugMode == 1:
		print ""
		print "(function:  get_data_from_single_hdf5) "
		print "data.shape = " + str(data.shape)
		print ""
	
	return data


'''
Description:
	WARNING:  this function currently works for unix-type file paths and needs to be adapted to also work with Windows
	Arguments:
		dirname:  (string) the path to the directory where the hdf5 data files are located
		dataLocations:  a list containing one or more strings; each string represents the hdf5 absolute path to a dataset which contains at least one column of metric data
		dataColumns:  a list of sublists; the number of sublists should equal the number of strings in the 'dataLocations' list; each sublist contains nonnegative integers corresponding to the column indices of the desired data for the corresponding dataset in 'dataLocations'
		minRow:  the first row of data (i.e. the index, starting at zero, of the first data point to be considered)
		maxRow:  the last row of data to be considered
		debugMode:  debug mode

'''
def get_all_directory_hdf5_data(dirname, file_ending, dataLocations, dataColumns, minRow=0, maxRow=-1, debugMode=0):
	if debugMode == 1:
		print ""
		print "(function:  get_all_directory_hdf5_data)"

	import fnmatch
	import os

	# create accumulating numpy array to hold all metric data
	data = np.array([])

	# create file name pattern
	pattern = '*.' + file_ending

	# iterate through all files with the specified file suffix (file_ending)
	numfiles = 0
	for file in os.listdir(dirname):
		if fnmatch.fnmatch(file, pattern):
			numfiles += 1
			filename = dirname + file
			if numfiles == 1:
				data = get_data_from_single_hdf5(filename, dataLocations, dataColumns, minRow, maxRow, 0)
			else:
				more_data = get_data_from_single_hdf5(filename, dataLocations, dataColumns, minRow, maxRow, debugMode)
				data = np.append(data, more_data, axis=0 )
			# if debugMode == 1:
			# 	print "Successfully processed file: " + filename

	if debugMode == 1:
		print "data.shape = " + str(data.shape)

	return data


if __name__ == '__main__':

	import sys

	dirname = sys.argv[1]
	file_ending = 'h5'

	# set which datasets, and which columns of those datasets, within the hdf5 file given by 'filename' should be used as coordinates for the data points
	metricDataLocations = ['data/E', 'data/Q']
	metricDataColumns = [[0], [0,1,2]]


	minRow = 0
	maxRow = 15000

	data = get_all_directory_hdf5_data(dirname, file_ending, metricDataLocations, metricDataColumns, debugMode=1)

