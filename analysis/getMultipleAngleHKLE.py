# The work flow of this program:
#  (1) iterate through a list of angles, (1.1) for each angle, read in the work_#/sim-#.nxs file and then call the single slice python script to convert, (1.2) save the resulting files in the initial directory

# NOTE:  you should open this program in a new empty directory which you want to populate with the csv data files from every slice in a particular simulation

import numpy as np
import sys
import getSingleAngleHKLE

def allAnglesHKLE(Ei):
	# get path to directory holding scattering angle folders (e.i, the directory holding all the "work_#" directories)
	path = sys.argv[1]

	minAngle = 0.0
	maxAngle = 1.0

	angle = minAngle
	increment = 3.0
	while angle <= maxAngle:
		stringAngle = str(angle)
		pathExtension = "/work_%s/sim-%s.nxs" % (stringAngle,stringAngle)
		slicePath = path + pathExtension
		getSingleAngleHKLE.rawInfo(slicePath, Ei)
		#print angle
		angle += increment



if __name__ == '__main__':
	# this program executed as script
	# do something:
	
	# define Ei (can later update to take this as a terminal argument maybe)
	Ei = 100.0
	allAnglesHKLE(Ei)
