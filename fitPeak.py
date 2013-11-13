import tables
import numpy as np
from matplotlib.pyplot import *
import sys
from scipy.optimize import curve_fit

from functions import *

#################################################################
# Check that correct number of arguments are passed
#################################################################
if (len(sys.argv)!=5):
        print "The script takes 4 arguments: \n"
        print "Arg 1: path to file containing spectrum \n"
        print "Arg 2: start of region where peak is located [channel number] \n"
        print "Arg 3: end of region where peak is located [channel number] \n"
        print "Arg 4: initial guess of fwhm [channels] \n"
        sys.exit();

#################################################################

#################################################################
# Assign input arguments to corresponding variables
#################################################################
# Set path to spectrum file
spath=sys.argv[1];

# Set region where peak is located. 
r1 = int(sys.argv[2]);
r2 = int(sys.argv[3]);  

# Set initial value of fwhm
fwhm_ini=int(sys.argv[4]);

#################################################################
# Begin analysis
#################################################################
# Load spectrum from textfile into an array
hist=np.loadtxt(spath, delimiter=" ", skiprows=2);
spec=hist[:,1];

# Bins are set to channel numbers  
bins=np.linspace(0, spec.shape[0], spec.shape[0]);

# Determine initial value of the centroid of the peak
# It is taken as the bin in the region defined by r1 and r2 with 
# the highest spectrum value
c_ini = r1+spec[r1:r2].argmax();

# Determine start and end point for fit
x1=c_ini-2*fwhm_ini;
x2=c_ini+2*fwhm_ini;

# Determine start area for fit
# The area is estimated by integrating the spectrum from x1 to x2
area_ini = spec[x1:x2].sum();

# Do fit 
# fit_param holds the fit parameters, and pcov the covariance matrix. 
# The diagonal of pcov corresponds to the variance of each fit parameter.
fit_param, pcov = gaussfit(spec, bins, x1, x2, area_ini, c_ini, fwhm_ini, 1);

# Print results
print "Results from fit:";
print "Area = ", fit_param[0], "Stdev = ", np.sqrt(pcov[0,0]);
print "Centroid = ", fit_param[1], "Stdev = ", np.sqrt(pcov[1,1]);
print "FWHM = ", fit_param[2], "Stdev = ", np.sqrt(pcov[2,2]);
print "Background start = ", fit_param [3], "Stdev = ", np.sqrt(pcov[3,3]);
print "Background slope = ", fit_param[4], "Stdev = ", np.sqrt(pcov[4,4]);

show()