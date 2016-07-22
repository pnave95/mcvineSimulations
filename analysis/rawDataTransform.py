# The purpose of this program is to convert tof and pixel position data to Qx,Qy,Qz,E  data (instrument coordinates)
# Then, we want to convert the instrument coordinates to reciprocal space coordinates based on the orientation of a crystal sample
#   for a given rotation angle of the sample

# status:  INCOMPLETE

import numpy as np


# compute vf (final speed)
# pixel_d = distance from sample to pixel, sample_d = distance from beam monitor 1 to sample, tof = time of flight of neutron (estimated)
def _vf(pixel_d, sample_d, tof):
	
