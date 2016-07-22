# The purpose of this program is to experiment with plotting a skewed Gaussian distribution

import numpy as np 
import scipy.special
from matplotlib import pyplot as plt


# define Gaussian cdf
def _Phi(x, mu, sigma):
	cdf = (1/2.0)*(1 + scipy.special.erf((x - mu)/(sigma*np.sqrt(2.0)) ))
	return cdf

# debugger:
#print str(_Phi(3.0))


# define Gaussian pdf
def _gaussian(x, mu, sigma):
	pdf = (1.0/(sigma*np.sqrt(2.0*np.pi) ))*np.exp(-(x - mu)**2/(2.0*sigma**2))
	return pdf

# debugger
#print str(_gaussian(3.0))

# define skew distribution
def _skewNormal(x, mu, sigma, alpha):
	xp = x*alpha
	f = 2.0*_gaussian(x, mu, sigma)*_Phi(xp, mu, sigma)  # Do I need to scale mu and sigma by alpha a well?
	return f



# Define parameters
alpha = 3.0
mu = 0.0
sigma = 0.1


# compute and plotting
x = np.linspace(-3.0, 3.0, num=120)
f = _skewNormal(x, mu, sigma, alpha)

plt.plot(x, f)
plt.ylim(0.0, 10.0)
plt.savefig("skewNormal_distribution_plot.png")