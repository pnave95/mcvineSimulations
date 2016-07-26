# The purpose of this program is to take the values of H and K for the Mourigal 2D SLAF simulation neutrons, and compute what the expected energy would be for that location.  Then, the MCViNE value of energy and this theoretical value of energy can be compared to compute the error

import numpy as np
import pandas as pd



def E_Q(H, K):
	gamma_k = 0.5*(np.cos(2.0*np.pi*H) + np.cos(2.0*np.pi*K))
	twoTheta = 0.200334842323
	E = 40.0*np.sqrt( (1.0+gamma_k)*(1 - gamma_k*np.cos(twoTheta) ))
	return E



minAngle = -90.0
maxAngle = 90.0

currentAngle = minAngle
increment = 3.0

while(currentAngle <= maxAngle):
	stringAngle = str(currentAngle)
	filename = "sim-" + stringAngle + "-Data.csv"

	# Read csv file into pandas data frame
	data = pd.read_csv(filename)

	# insert new columns
	data.insert(4, 'error_E', 0.0)
	data.insert(5, 'error_E_squared', 0.0)
	data.insert(6, 'scalarQ', 0.0)

	# Now iterate through every row of the data
	for i in range (data['H'].count()):
		Hval = data['H'][i]
		Kval = data['K'][i]
		E = data['E'][i]
		Etrue = E_Q(Hval, Kval)
		error = E - Etrue
		error2 = error**2
		data['error_E'][i] = error
		data['error_E_squared'][i] = error2
		Qsquared = (data['Qx'][i])**2 + (data['Qy'][i])**2 + (data['Qz'][i])**2
		data['scalarQ'][i] = np.sqrt(Qsquared)

	# save new file
	# all data needed should now be available
	newFileName = "results-" + stringAngle + ".csv"
	data.to_csv(newFileName)

	# increment angle
	currentAngle += increment
