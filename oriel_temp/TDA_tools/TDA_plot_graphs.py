'''
This program is designed for plotting graphs and calculating graph colorings
'''
import numpy as np
import pygraphviz as pgv



'''
Description:
	Summary:  this function takes a data column (which may be the same as the filter function column or different), computes the range of the data in that column;  then, it averages the values for all data points that are members of a given node, scales into an appropriate coloring range, and returns that color for that particular node. 
	Returns:
		a 3-element list containing the RGB values for the given input node
'''
def get_graph_node_colors_by_column( node_element_indices, coloringDataColumn, debugMode=0):
	# if in debug mode, print the name of this function to make it clear where the program is working
	if debugMode == 1:
		print ""
		print "(function:  get_graph_node_colors_by_column)  "

	# check if node_element_indices is empty
	if node_element_indices == set():
		if debugMode == 1:
			print "Node with no substituent elements detected!"
		return [0, 255, 0]
	else:
		# get min and max values of the coloring function
		cmin = np.amin(coloringDataColumn)
		cmax = np.amax(coloringDataColumn)
		# compute range of color function values
		R = cmax - cmin

		# extract elements of coloringDataColumn which are also elements of the node
		node_cf_values = []
		for i in node_element_indices:
			cf_value = coloringDataColumn[i]
			node_cf_values.append(cf_value)

		# now, normalize color function values to lie in [0,1]
		node_cf_values = np.array(node_cf_values)
		node_cf_values = (node_cf_values - cmin) / R

		# create normalized intensity values for the coloring function as defined on each data point
		intensity_norm = np.sum(node_cf_values) / len(node_cf_values)
		
		# create RBG values
		red_value = np.rint( intensity_norm*255 )

		blue_value = 255 - red_value

		green_value = 1.0 - (np.sqrt((intensity_norm - 0.5)**2)/0.5)
		green_value = np.rint( green_value*100 )


		if debugMode == 1:
			print "np.amin(red_value) = " + str(np.amin(red_value))
			print "np.amax(red_value) = " + str(np.amax(red_value))
			print ""

		return [red_value, green_value, blue_value]


'''
Description:
	Summary:  takes an abstract graph and a column of coloring function values and produces a graph using a pygraphviz algorithm.  The resulting graph is then saved as an image file
'''
def plot_pygraphviz(graph, graph_filename, coloringColumn, debugMode=0):
	
	# separate nodes and edges
	nodes = graph[0][0]
	edges = graph[0][1]
	clusters = graph[1]

	# create an empty graph
	G = pgv.AGraph()

	# set default attributes
	G.graph_attr['label'] = 'TDA graph'
	G.node_attr['shape'] = 'circle'
	G.node_attr['style'] = 'filled'

	# add disconnected nodes
	for i in range(len(nodes)):
		vertex = nodes[i]
		clusterPoints = clusters[i]
		#red, green, blue = TDA_plot_graphs.get_graph_node_colors_by_column(clusterPoints, coloringColumn)
		red, green, blue = get_graph_node_colors_by_column(clusterPoints, coloringColumn)
		
		node_color = "#%2x%2x%2x" %(red,green,blue)
		G.add_node(str(vertex), color=node_color)


	# add edges
	for e in edges:
		# get both boundary vertices
		v1 = str(e[0])
		v2 = str(e[1])
		# add the edge
		G.add_edge(v1, v2)


	# add positions to the nodes using a Graphviz layout algorithm
	G.layout() # default to neato

	# render graph to an image
	G.draw(graph_filename)


'''
Description:
	Summary:  plots a graph using matplotlib and networkX.  The graph layout algorithms don't seem quite as sophisticated as pygraphviz
'''
def plot_matplotlib_networkX(abstract_graph, plot_filename, debugMode=0):
	import networkx as nx 
	from matplotlib import pyplot as plt

	# get node and edge sets
	nodes = abstract_graph[0][0]
	edges = abstract_graph[0][1]

	# create networkX graph
	G = nx.Graph()

	# add nodes from list of nodes (could also add from a set, etc)
	G.add_nodes_from(nodes)

	# add edges from a list
	G.add_edges_from(edges)

	# draw graph
	nx.draw(G)
	#nx.draw_circular(G)
	#nx.draw_spectral(G)
	plt.show()
