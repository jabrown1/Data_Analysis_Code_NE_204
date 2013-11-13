import numpy as np
from scipy.optimize import curve_fit
from matplotlib.pyplot import *


#####################################################################
# Function that does trapezoidal filtering on a number of pulses
def trap_filter(pulses, k, m, M):
# pulses is an array with dimensions [NumberOfPulses SamplesPerPulse]
# k = peaking time
# m = peaking time + gap time
# M = decay correction parameter
# Returns an array containing the shaped pulses, with the same size as pulses.

        #Adders
        add=np.zeros_like(pulses);
        add+=pulses;
        add[:,k:]-=pulses[:,:-k];
        add[:,m:]-=add[:,:-m];

        #Accumulators
        acc=np.zeros_like(pulses);
        acc+=add;
        acc=acc.cumsum(axis=1);
        acc+=M*add;
        acc=acc.cumsum(axis=1);

        # Returns the filtered pulses
        return acc;

###################################################################

# Function that saves a spectrum in a txt file
def save_spectrum(histogram, filename):
# histogram = array to be saved
# filename= name of output file
        spec = open(filename,'w');        
        spec.write("# Col 0: Lower bin edge \n");
        spec.write("# Col 1: Number of counts \n");

        for line in histogram:
                spec.write(" ".join(str(x) for x in line) + "\n");
        spec.close();
        return;

#####################################################################

# Gauss function with linear background
def gauss_lin_bg(x, a, b, c, y1, slope):
        # x = values where function is evaluated (a vector)
        # a = area of the peak
        # b = centroid 
        # c = FWHM 
        # y1 = start y value of linear background
        # slope = slope of linear background
        arg=2.35*(x-b)/c;
        gauss=(0.9375*a/c)*np.exp(-0.5*arg*arg);
        bg=y1+slope*(x-x[0]);
        return gauss+bg;

####################################################################

# Function that fits a gaussian + linear background to data
def gaussfit(spectrum, bins, x1, x2, area, centroid, fwhm, plot_flag):
        # spectrum = histogram data to be fitted 
        # x1 = start of fit region
        # x2 = end of fit region
        # area = start value of area
        # centroid = start value of centroid
        # fwhm = start value of fwhm
        # plot_flag = if set to 1 the result of the fit will be plotted

        # Determine start values at ends of fit region
        y1 = np.mean(spectrum[(x1-fwhm):(x1+fwhm)]);
        y2 = np.mean(spectrum[(x2-fwhm):(x2+fwhm)]);

        # Determine start value for slope
        slope = (y2-y1)/(x2-x1);

        # Do fit
        # fit_param holds the fit parameters, and pcov the covariance matrix. 
        # The diagonal of pcov corresponds to the variance of each fit parameter.
        fit_param, pcov = curve_fit(gauss_lin_bg, bins[x1:x2], spectrum[x1:x2], \
        p0=[area, centroid, fwhm, y1, slope]);

        # Calculate fitted curve
        fitted_curve=gauss_lin_bg(bins[x1:x2], fit_param[0], fit_param[1], fit_param[2], fit_param[3], fit_param[4]);
        # Calculate difference between fitted curve and data
        rest=fitted_curve-spectrum[x1:x2];

        if (plot_flag == 1):
                # Plot results
                figure(1)
                plot(bins,spectrum)
                plot(bins[x1:x2], fitted_curve, 'r');
                plot(bins[x1:x2], rest);

        return fit_param, pcov;

####################################################################

# Function that saves multiple spectra to a text file
def save_multiple_spectra(histogram, kg_array, filename):
# histogram = array to be dumped to file
# filename= name of output file
        spec = open(filename,'w');        
        spec.write("########################### \n");
        spec.write("# File containing multiple spectra created with different shaping times \n");
        spec.write("########################### \n");
        spec.write("# Shaping times: \n");
        spec.write("# Row 1: Peaking time for each spectrum \n");
        spec.write("# Row 1: Gap time for each spectrum \n");
        
        for line in kg_array:
                spec.write(" ".join(str(x) for x in line) + "\n");

        spec.write("########################### \n");
        spec.write("# Bin edges and spectra: \n");
        spec.write("# Col 1: bin lower bin edge \n");
        spec.write("# Col 2-end: number of counts \n");

        for line in histogram:
                spec.write(" ".join(str(x) for x in line) + "\n");

        spec.close();
        return;