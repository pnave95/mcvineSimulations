#  The purpose of this program is to remove beam simulations which didn't work

import numpy as np
import os


# possible fermi chopper packages
fermiChoppers = ['100-1.5-SMI', '700-1.5-SMI', '700-0.5-AST']

# possible fermi chopper frequencies
fermiFrequencies = [300, 360, 420, 480, 540, 600]

# possible T0 chopper frequencies
T0Frequencies = [60, 90, 120]

# range of Ei's to test:
LowEiRange = np.linspace(15.0, 85.0, num=5)
HighEiRange = np.linspace(100.0, 1500.0, num=5)

# define current directory
workdir = os.path.abspath('.')

for chopper in fermiChoppers:
	for fermi in fermiFrequencies:
		for T0 in T0Frequencies:
			for E in LowEiRange:
				# change back to starting directory
				os.chdir(workdir)
				# determine beam sim directory name
				dirname = "beam-E-" + str(E) + "-chopper-" + str(chopper) + "-fermi-" + str(fermi) + "-T0-" + str(T0) + "-INCORRECT"
				# check if directory exists; if so, delete it
				directory = os.path.abspath(dirname)
				if os.path.isdir(directory):
					remove = "rm -r " + str(dirname)
					os.system(remove)				
				







