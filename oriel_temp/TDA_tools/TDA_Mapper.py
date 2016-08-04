# The purpose of this program is to attempt to implement the TDA Mapper algorithm

# current status: Appears to be working

# debug modes:
#  0 = none, 1 = debug, 2 = verbose debug (not yet implemented)

'''
TODO:  
- implement default debugMode=0
- change the location of the "(function:  functionName)  " debug status print to occur at the beginning of the function call instead of the end; this way debug messages may be called anywhere throughout the function in an identifiable manner
'''

import numpy as np
import h5py as hf
import sys
import fastcluster
import pygraphviz as pgv

# import modules for TDA:
import TDA_plot_graphs
import TDA_dataFromHDF5



'''
Description:
	Arguments (mandatory):
		filterImage:  (expects) a 1-dimensional non-empty array with the values of the filter function for each data point
		N: = number of intervals in the cover 
		x: = decimal percent overlap with one neighbor
	Arguments (optional):
		debugMode:  (default 0)  specify debug level
	Returns:
		a list with four elements: [filterImage_min, filterImage_max, interval_length, offset]
		* offset is the change in starting position of one element of the cover to the next
'''
def get_cover_params(filterImage, N, x, debugMode=0):

	if debugMode == 1:
		print ""
		print "(function:  get_cover_params)  "

	fmin = np.amin(filterImage)				# cover_min = filterImage_min
	fmax = np.amax(filterImage)				# cover_max = filterImage_max
	R = fmax - fmin	# filterImage_range
	L = R / ( N*(1.0 - x) + x)				# interval_length
	#L = R / (N*(1.0 - x))
	t = L*(1.0 - x)							# offset

	if debugMode == 1:
		print "R = " + str(R)
		print "[cover_min, cover_max, interval_length, offset] = " + str([fmin, fmax, L, t])
		print ""

	return [fmin, fmax, L, t]


'''
Description:
	Arguments:
		filterImage: (expects) a 1-D numpy array containing the filter function image; the cover parameters; and the coverElementNumber, (an integer 0, 1, ... ) which determines where to start from in the range of the image
	Returns:
		coveringSetElements:  a list of the indices of each data point which has an image in the range of the specified element of the cover
		*NOTE:  the indices in the sublists will be with respect to the original 'filterImage' numpy array
'''
def get_single_cover_element(filterImage, params, coverElementNumber,debugMode=0):
	
	if debugMode == 1:
		print ""
		print "(function:  get_single_cover_element)  "

	# determine how many data points are in the image
	N = len(filterImage)

	# determine bounds of this element of the cover
	cover_min = params[0]
	cover_max = params[1]
	L = params[2]
	t = params[3]
	coverElementMin = cover_min + float(coverElementNumber)*t
	coverElementMax = coverElementMin + L

	# create list to hold indices of each element of the image which is within the range of this cover element
	coveringSetElements = []

	# iterate through the filter function image
	for i in range(N):
		if filterImage[i] > coverElementMin:
			if filterImage[i] < coverElementMax:
				coveringSetElements.append(i)

	# debug / test
	if debugMode == 1:
		print "image range = [" + str(cover_min) + ", " + str(cover_max) + "]"
		print "coverElementNumber = " + str(coverElementNumber) + ": [" + str(coverElementMin) + ", " + str(coverElementMax) + "]"
		print "len(coveringSetElements) = " + str(len(coveringSetElements))
		print ""

	return coveringSetElements



'''
Description:
	inputs:
		pointCoordinates: the numpy array of metric coordinates for all data points
		coveringSetElements:  a list which contains indices of all points in a particular set which is an element of the cover

	returns:
		a numpy array which is a subset of 'pointCoordinates', but which only contains those elements which are in the coveringSetElements
'''
def get_single_inverse_image(pointCoordinates, coveringSetElements, debugMode=0):

	if debugMode == 1:
		print ""
		print "(function:  get_single_inverse_image)  "

	# get the total number of data points
	N = pointCoordinates.shape[0]		# TODO:  handle case of 1-dimensional coordinates (will this cause a problem or not?)

	inverseImage = []
	for i in coveringSetElements:
		inverseImage.append(pointCoordinates[i,:])
		

	inverseImage = np.array(inverseImage)

	if debugMode == 1:
		print "total number of data points: N = " + str(N)
		print "inverseImage.shape = " + str(inverseImage.shape)
		print ""

	return inverseImage


'''
Description:
	Inputs:
		inverseImage:  the numpy array which consists of the data point coordinates for all points within a single element of the cover
		coveringSetElements:  a list containing the original indices of each data point in the inverse image
		metricName:  which metric to use for measuring distance between data points (valid options are those supported by the fastcluster module)
		clusterMethod:  a string (which will be passed to the fastcluster module) specifying which clustering algorithm to use
		dmax = maximum distance to cluster data points before deciding what the final clusters produced are.  
				* Currently this is quite arbitrary, and there is definitely a need for improvement in how this is selected / handled
	Returns:
		clusters:  a list of sets; each set corresponds to one cluster in the inverse image, and the contents of that set are the original indices of the data points in that cluster
'''
def get_single_inverse_image_clusters(inverseImage, coveringSetElements, dmax, debugMode, metricName='euclidean', clusterMethod='single'):

	if debugMode == 1:
		print ""
		print "(function:  get_single_inverse_image_clusters)  "


	# check that at least two data points are in the inverse image
	numElements = len(coveringSetElements)
	if numElements == 0:
		return [set()]
	elif numElements == 1:
		return [{coveringSetElements[0]}]
	else:
		# perform clustering in the inverse image
		links = fastcluster.linkage_vector(inverseImage, method=clusterMethod, metric=metricName)

		# determine total number of data points in this inverse image
		N = inverseImage.shape[0]

		# count the number of clusters; to start, no points have been merged so number of Clusters =  number of data points = N
		numClusters = N

		# create an initial dictionary of clusters; the labels will be in the range [0, N-1] and these will correspond to how the initial clusters (the nodes) are labeled in the 'links' array.  The values will be sets; to start, each set will be a singleton that contains the original point index as its only element.
		clusterDict = {i: [coveringSetElements[i]] for i in range(N)}

		# The third column (index 2) of 'links' contains the merging distance for the two clusters listed in that row.  We will now loop through  'links' until that merging distance is greater than dmax, or until we reach the end of 'links', whichever comes first.
		index = 0	# this is out "counter" index
		while index < (N-1):
			if links[index][2] < dmax:
				# get the fastcluster indices for the two nodes merged at this step
				p1 = links[index][0]
				p2 = links[index][1]
				# Now create a new cluster (that is, a new dictionary entry) which has a label that is one higher than whatever is currently the highest label and which has a value that is the list produced by combining the lists which correspond to clusters p1 and p2 (remember, p1 and p2 are also labels in the dictionary)
				newLabel = N + index	# start at N, then increase
				newCluster = clusterDict[p1] + clusterDict[p2]	# this concatenates the lists which are indexed in the dictionary by p1 and p2
				clusterDict.pop(p1); clusterDict.pop(p2)	# remove these clusters
				clusterDict[newLabel] = newCluster 			# add the new cluster which was formed by merging the old 2

			index += 1 	# increment our index counter after each merge

		# Take each cluster, which is stored in a a dictionary value, and place it into a list.  This list of sublists now contains a sublist which holds the original indices of each point in the cluster which is represented by that sublist
		clusters = clusterDict.values()

		# EXPERIMENT: try converting each element of 'clusters' (elements being lists) to a set
		for i in range(len(clusters)):
			listlength = len(clusters[i])
			clusters[i] = set(clusters[i])
			setlength = len(clusters[i])
			# if lengths don't match, then something went wrong
			if listlength != setlength and debugMode == 1:
				print "ERROR!! (function:  get_single_inverse_image_clusters):  listlength != setlength"

		if debugMode == 1:
			print "links.shape = " + str(links.shape)
			print "number of data points in this inverse image = " + str(N)
			print "should = number of elements in coveringSetElements = " + str(len(coveringSetElements))
			print "Number of clusters = " + str(len(clusters))
			print "first merging distance = " + str(links[0][2])
			print "last merging distance = " + str(links[N-2][2])
			print ""

		return clusters


'''
Description:
	Arguments (mandatory):
		dataCoordinates:  numpy array of all data point metric coordinates
		filterImage:  the column of data representing the value of the filter function for each data point
		intervals:  number of intervals in cover 
		CoverParams:  a 4-element list in the form returned by the 'get_cover_params()' function
		dmax:  max cluster merging distance (purely heuristic and quite arbitrary right now)
	Arguments (optional):
		debugMode: (default 0)
		metricName:  a string representing the metric to be used when computing distance between two data points.  Acceptible metrics are those accepted by the 'fastcluster' module
		clusterMethod:  a string representing the clustering algorithm to use (by 'fastcluster' module)
	Returns:
		Allclusters:  a list of sets; each set represents a cluster in the inverse image of some element of the cover of the range of the filter function.  Elements of these sets are original indices of data points
'''
def get_all_clusters(dataCoordinates, filterImage, intervals, CoverParams, dmax, debugMode=0, metricName='euclidean', clusterMethod='single'):

	if debugMode == 1:
		print ""
		print "(function:  get_all_clusters)  "

	# create a list which will hold the sets corresponding to each cluster
	# NOTE:  implementing the logic this way (lumping clusters from all cover elements together) is inefficient but simple;  there is definitely room for improvement here
	Allclusters = []

	# iterate through the cover.  For each element of the cover, take the inverse image, cluster, and then return the clusters; add them to the list 'Allclusters'
	for c in range(intervals):
		# get list of original indices of each data point in this element of the cover
		covering_set_elements = get_single_cover_element(filterImage, CoverParams, c, 0)
		# Now get the inverse image of this element of the cover (this covering set)
		inverseImage = get_single_inverse_image(dataCoordinates, covering_set_elements, 0)
		# Next get the clusters in that inverse image
		clusters = get_single_inverse_image_clusters(inverseImage, covering_set_elements, dmax, 0)
		# now, add those clusters to 'Allclusters'
		Allclusters = Allclusters + clusters



	if debugMode == 1:
		print "total number of clusters = len(Allclusters) = " + str(len(Allclusters))
		print ""

	return Allclusters


'''
Description:
	Summary:  This function finds all pairs of distinct clusters which share at least one point; edges are then created for each such pair of clusters
'''
def get_edges(all_clusters, debugMode=0):

	if debugMode == 1:
		print ""
		print "(function:  get_edges)  "

	# determine total number of clusters
	numClusters = len(all_clusters)

	# make a list to collect edges in
	edges = []

	# find all set intersections (sets represent clusters)
	for i in range(numClusters):
		for j in range(numClusters):
			if i != j:		# we only want to compare distinct clusters
				if len( all_clusters[i] & all_clusters[j] ) != 0:
					# then we have an intersection!  Add an edge!
					edges.append([i,j])


	if debugMode == 1:
		print "total number of clusters = numClusters = " + str(numClusters)
		print "total number of edges = " + str(len(edges))
		print ""

	return edges




'''
Description:
	Arguments:
		dirname: (string)  directory where data files are located
		file_suffix: (string)  file ending (e.g. 'h5');  currently only 'h5' is acceptable
		metricDataLocations:  
	Returns:
		1-dimensional simplicial complex (i.e. a graph with edges and nodes)
		* AND a list of sets which represent the data points in each cluster (each node);  these two objects are together in a list
'''
# EXPERIMENT:  try taking 'dirname' in place of 'filename', and also take 'file_suffix'
def Mapper(dirname, file_suffix, metricDataLocations, metricDataColumns,minRow, maxRow, filterDataName, filterDataColNumber, intervals, overlap, dmax, debugMode=0):

	if debugMode == 1:
		print ""
		print "(function:  Mapper)  "

	# get metric data from hdf5;  this should later be changed so that options are available

	coordinates = TDA_dataFromHDF5.get_all_directory_hdf5_data(dirname, file_suffix, metricDataLocations, metricDataColumns, minRow, maxRow, debugMode)

	# get filter function column
	# TODO: add ability to use general filter functions and compute the corresponding column of filter function data
	
	filterCol = TDA_dataFromHDF5.get_all_directory_hdf5_data(dirname, file_suffix, filterDataName, filterDataColNumber, minRow, maxRow, debugMode)

	# get the parameters of the cover
	cover_params = get_cover_params(filterCol, intervals, overlap, debugMode)

	# get a list of sets (sets represent clusters)
	clusters = get_all_clusters(coordinates, filterCol, intervals, cover_params, dmax, debugMode)

	# make a list of nodes (just numbers from 0 to numClusters - 1)
	nodes = range(len(clusters) )

	# check for intersections between distinct clusters; intersections would be shared points.  For any two distinct clusters which share points, we add an edge
	edges = get_edges(clusters, debugMode)

	# right now we are only computing a graph, but this can later be extended to abstract simplicial complexes
	simplicial_complex = [nodes, edges]

	if debugMode == 1:
		print "total number of clusters = len(nodes) = " + str(len(nodes))
		print "total number of edges = " + str(len(edges))
		print ""

	return [simplicial_complex, clusters]



# if run as a script, perform this code
if __name__ == '__main__':
	
	import time
	start = time.time()

	# parameters for building / debugging
	debugMode = 1
	dirname = sys.argv[1]

	# set which datasets, and which columns of those datasets, within the hdf5 file given by 'filename' should be used as coordinates for the data points
	metricDataLocations = ['data/E', 'data/Q']
	metricDataColumns = [[0], [0,1,2]]

	# TODO:  implement safeguard so that minRow cannot accidentally be set to something larger than exists
	minRow = 500
	maxRow = 1000
	#maxRow = -1
	filterNames = ['data/E']	
	filterColumns = [[3]]	# this is squared energy error

	# set numIntervals and decimal percent overlap
	#intervals = 19
	intervals = 50
	#overlap = 0.5
	overlap = 0.4

	d_max = 2.7

	# get coloring data
	#colorCol = get_filter_col_from_hdf5(filename, filterDataName, filterDataColNumber)
	colorCol = TDA_dataFromHDF5.get_all_directory_hdf5_data(dirname, 'h5', filterNames, filterColumns)



	# debug / test

	augmented_simplicial_complex = Mapper(dirname, 'h5',metricDataLocations, metricDataColumns, minRow, maxRow, filterNames, filterColumns, intervals, overlap, d_max, debugMode)

	#graph = augmented_simplicial_complex[0]

	graph_filename = 'TDA_Mapper_graph.png'


	TDA_plot_graphs.plot_pygraphviz(augmented_simplicial_complex, graph_filename, colorCol)
	TDA_plot_graphs.plot_matplotlib_networkX(augmented_simplicial_complex, graph_filename, debugMode)

	end = time.time()

	print "time elapsed = " + str(end - start)
