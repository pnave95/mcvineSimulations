import numpy as np
ostream = open('submit.sh', 'wt')
for a in np.arange(-90, 90.1, 1.):
    ostream.write('./scripts/sim.py --angle=%s \n' % a)
    continue
ostream.close()