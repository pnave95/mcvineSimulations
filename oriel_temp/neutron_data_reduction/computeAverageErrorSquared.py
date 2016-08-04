# the purpose of this program is to determine the average squared error of a simulation

import numpy as np
import h5py as hf
import sys


def computeErrors(ifile):

	# debugger:
	#print "ifile = " + ifile
	
	# read in file
	f = hf.File(ifile, "r")

	EarrayHDF = f.get('data/E')
	WarrayHDF = f.get('data/weights')
	Earray = np.array(EarrayHDF)
	Warray = np.array(WarrayHDF)

	N = len(Warray)
	sumErr2 = 0.0
	sumWeights = 0.0

	sumErrNonTrivialDispersion = 0.0
	sumNonTrivialWeights = 0.0

	for i in range(N):
		#if Warray[i] > 0.05:
		err2 = Earray[i][3]*Warray[i]
		sumErr2 += err2
		sumWeights += Warray[i]

		if Earray[i][2] > 0.0:
			sumErrNonTrivialDispersion += sumErr2
			sumNonTrivialWeights += sumWeights

	meanErr2 = sumErr2 / sumWeights
	print "mean error squared = " + str(meanErr2)


	meanNonTrivialErr2 = sumErrNonTrivialDispersion / sumNonTrivialWeights
	print "For inelastic dispersion, mean error squared = " + str(meanNonTrivialErr2)


if __name__ == '__main__':

	ifile = sys.argv[1]
	computeErrors(ifile)




