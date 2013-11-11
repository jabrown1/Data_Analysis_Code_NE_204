import tables
import sys
import numpy as np;
import time;
from scipy import signal
from matplotlib.pyplot import *
from scipy.optimize import curve_fit

def load_data(filename):
	hf = tables.openFile(filename);
	return hf;

def read_n_m_lines(f,n,m,pl):
	
	g = f.root.RawData.read(n,m);

	for i in range(0,m-n):
		g[i,:]-=np.average(g[i,pl/2-80:pl/2-50]);
		#g[i] = g[i,(pl/2-80):]
        return g;

def create_trap_filter(pl,p,g,cf):
    a=np.zeros((pl,pl));
    for n in range(0,pl):
        if (n<=p):
            v=p;
            for j in range(0,n):
                a[j,n]=cf+v;
                v-=1;
        elif (n<p+g):
            v=p;
            for j in range(0,n):
                a[j,n]=cf+v;
                v-=1;
            for k in range(0,n-p):
                a[k,n]=p;

        elif (n<2*p+g):
            v=p;
            for j in range(n-p,n):
                a[j,n]=cf+v;
                v-=1;
            for k in range(n-p-g,n-p):
                a[k,n]=p;
            v=-(cf)
            for l in range(0,n-p-g):
                a[l,n]=v;
                v+=1;
        else:
            v=p;
            for j in range(n-p,n):
                a[j,n]=cf+v;
                v-=1;
            for k in range(n-p-g,n-p):
                a[k,n]=p;
            v=-(cf);
            for l in range(n-(2*p+g),n-p-g):
                a[l,n]=v;
                v+=1;
    return a; 
def analyze_spectrum(filename,peaking_time,gap_time,cf, resolution,max_bin,mode):
    f = load_data(filename)
    length = f.root.RawData.nrows;
    pulse_length = f.root.RawData.rowsize/4
    print length
    print pulse_length
    
    a = create_trap_filter(pulse_length,peaking_time,gap_time,cf)

    spectra = np.histogram(0,bins=(max_bin/resolution),range = (0,max_bin))
    
    
    
    
    
    
    if (mode == 1):
        g = read_n_m_lines(f,0,10000,pulse_length)
        z = vector_product(g,a)
        z = z/(2*peaking_time)**2
        x = np.histogram(z[:,(pulse_length/2+peaking_time+gap_time-9*gap_time/10)],bins = (max_bin/resolution),range = (0,max_bin))        
        
        
        peaks = np.array(signal.find_peaks_cwt(x[0],np.arange(5/resolution,20/resolution),min_snr=3))
        peaks = resolution * peaks
        print peaks

        #plot(x[1][1:],x[0])
        #show()
        
        #a1 = raw_input("Enter something: ")
        #a2 = raw_input("Enter something: ")
        #a3 = raw_input("Enter something: ")
        #peaks = np.array([float(a1),float(a2),float(a3)])        
        
        a1 = np.array([[peaks[0],1],[peaks[1],1],[peaks[2],1]])
        b1 = np.array([59.5,1173,1332])
        print a
        cor = np.linalg.lstsq(a1,b1)[0]
    

        
        for chunk in range(10000, length, 10000):
            g = read_n_m_lines(f,chunk-10000,chunk,pulse_length)
            z = vector_product(g,a)
            z = z/(2*peaking_time)**2
            z = z*cor[0]+cor[1]
            #z = z*59.5/peaks[0]
	    test = z[:,pulse_length/4]>4
	    z = z[np.invert(test),:]
	    x = np.histogram(z[:,(pulse_length/2+peaking_time+gap_time-9*gap_time/10)],bins = (max_bin/resolution),range = (0,max_bin))
            
            spectra[0][:] = np.add(spectra[0],x[0])

    if (mode == 0):
        
        g = read_n_m_lines(f,0,10000)
        z = vector_product(g,a)
        z = z/(peaking_time+gap_time)**2
        x = np.histogram(np.amax(z,axis=1),bins = (max_bin/resolution),range = (0,max_bin))        


        
        peaks = np.array(signal.find_peaks_cwt(x[0],np.arange(5/resolution,20/resolution),min_snr=3))
        peaks = resolution * peaks

       # plot(x[1][1:],x[0])
       # show()
        
        print peaks
        a1 = raw_input("Enter something: ")
        a2 = raw_input("Enter something: ")
        #a3 = raw_input("Enter something: ")
        peaks = np.array([float(a1),float(a2),float(a3)])

        
        a1 = np.array([[peaks[0],1],[peaks[1],1],[peaks[2],1]])
        b1 = np.array([1173,1332])
        #print a
        cor = np.linalg.lstsq(a1,b1)[0]
        
        
        for chunk in range(10000, length, 10000):
            g = read_n_m_lines(f,chunk-10000,chunk)
            z = vector_product(g,a)
            z = z/(peaking_time+gap_time)**2
            z = z*cor[0]+cor[1]
            x = np.histogram(np.amax(z,axis=1),bins = (max_bin/resolution),range = (0,max_bin))

            spectra[0][:] = np.add(spectra[0],x[0])
                        
    
    #plot(spectra[1][1:],spectra[0])
    #show()
    peaks = np.array(signal.find_peaks_cwt(spectra[0],np.arange(5,20),min_snr=3))
    peaks = resolution * peaks
    print peaks;
    return spectra;
                  
def vector_product(x,y):
    start=time.time();
    z = np.dot(x,y);
    end=time.time();
    print end-start;
    return z;

def trapezoidalFilter(k,m,t,v):
    
    '''
    Takes signal v with rise time k (in samples) and gap time m (in samples) and returns the filtered signal.
    '''
        
    M = t
    
    s = np.zeros_like(v)
    p = np.zeros_like(v)
    r = np.zeros_like(v)
    d = np.zeros_like(v)
    n = np.argmax(v)
    length = np.size(v)
    
    for n in range(2*k+m,length):
        d[n] = v[n] - v[n-k] - v[n-m] + v[n-k-m]
        p[n] = p[n-1] + d[n]
        r[n] = p[n] + M*d[n]
        s[n] = s[n-1] + r[n]
        
    return s

def calib_spectrum(current):
    test = np.array(current)
    a = np.array([[current[0],1],[current[1],1],[current[2],1]])
    b = np.array([59,1173,1332])
    #print a
    cor = np.linalg.lstsq(a,b)[0]
    
    spectrum[1][:] = spectrum[1][:]*cor[0]+cor[1]
    
    
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


# Function that fits a gaussian + linear background to data
def gaussfit(spectrum, bins, x_1, x_2, area, centroid, fwhm, plot_flag):
	# spectrum = histogram data to be fitted 
	# x1 = start of fit region
	# x2 = end of fit region
	# area = start value of area
	# centroid = start value of centroid
	# fwhm = start value of fwhm
	# plot_flag = if set to 1 the result of the fit will be plotted
    print x_1
    x1 = np.where(bins>x_1)[0][0]
    print x1
    x2 = np.where(bins>x_2)[0][0]
	# Determine start values at ends of fit region
    y1 = np.mean(spectrum[(x1-fwhm):(x1+fwhm)]);
    y2 = np.mean(spectrum[(x2-fwhm):(x2+fwhm)]);

	# Determine start value for slope
    slope = (y2-y1)/(x2-x1);

	# Do fit
	# fit_param holds the fit parameters, and pcov the covariance matrix. 
	# The diagonal of pcov corresponds to the variance of each fit parameter.
    fit_param, pcov = curve_fit(gauss_lin_bg, bins[x1:x2], spectrum[x1:x2], 
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


