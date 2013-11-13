
from numba import autojit,codegen
import numpy
import matplotlib.pyplot as plt
import logging
import scipy.optimize as opt
import time
import scipy.stats as stats
codegen.debug.logger.setLevel(logging.INFO)

def readFile(fname,ftype):
    
    '''Reads in file in the current directory of type 'ftype' (string).  Can read HDF5, plain text and ROOT trees.  
    
    fname:  string containing the file name with extension
    
    ftype:  string indicating file type.  Use:
        'HDF5': .h5 file type containing HDF5 formatted information
        'plain': plain text information with ' ' delimeter
        'root': ROOT tree formatted information
        
    '''
    
    if ftype == 'HDF5':
        import h5py
        f = h5py.File(fname,'r',dtype=int)
        data = f['RawData']
        return data
    
    if ftype == 'plain':
        f = numpy.loadtxt(fname)
        return f
    
    if ftype == 'root':
        import ROOT
        f = ROOT.fopen(fname,'r')
        return f

def triangularFilter(k,v):
    
    '''
    Takes signal v with rise time k (in samples) and returns the filtered (triangular) signal.
    '''
    
    import numpy as np

    s = numpy.zeros_like(v)
    for n in range(len(v)):
        try:
            s[n] = s[n-1] + v[n] - 2*v[n-k] + v[n-2*k]
        except:
            s[n] = v[n]

    return s

@autojit
def trapezoidalFilter(k,m,t,v):
    
    '''
    Takes signal v with rise time k (in samples) and gap time m (in samples) and returns the filtered signal.
    '''
    
    M = t
    s = numpy.zeros_like(v)
    p = numpy.zeros_like(v)
    r = numpy.zeros_like(v)
    d = numpy.zeros_like(v)
    length = numpy.size(v)
    
    for n in range(2*k+m,length):
        d[n] = v[n] - v[n-k] - v[n-k-m] + v[n-2*k-m]
        p[n] = p[n-1] + d[n]
        r[n] = p[n] + M*d[n]
        s[n] = s[n-1] + r[n]
        
    return s/k, max(s)/k

def subtractBaseline(v):
    
    '''
    Takes in signal v and subtracts baseline by taking the average of samples before a sharp rise.
    '''
    mu = numpy.mean(v[:100])

    v[:] -= mu*numpy.ones(len(v))
        
    return v

def fitExp(v):
    n = numpy.argmax(v)+25
    xdata = numpy.arange(len(v[n:]))
    ydata = numpy.array(v[n:])
    def exp_func(x,a,b):
        return a*numpy.exp(-x/b)
    fit,fit2 = opt.curve_fit(exp_func,xdata,ydata)
    return int(numpy.round(fit[1]))

def findRise(v, findPeak=False, thresh=100):
    
    Dv = numpy.zeros(len(v))
    Dv2 = numpy.zeros(len(v))
    peak_locs = []
    
    for j in range(len(v)-1):
        Dv[j] = v[j+1]-v[j]
        
    if findPeak:
        try:
            for j in range(len(v)-1):
                Dv2[j] = Dv[j+1]-Dv[j]
            peaks = [sample<-thresh for sample in Dv2]
            peaks = list(numpy.nonzero(peaks)[0])
            peak_locs.append(peaks[0])
            for i in range(1,len(peaks)):
                if numpy.abs(peaks[i]-peaks[i-1])>10:
                    peak_locs.append(peaks[i])
            else:
                return numpy.array(peak_locs)
        except:
            print 'Error: There were no peaks found in this histogram. Note: threshold may be too high; try lowering.'
    else:
        return numpy.argmax(v)

def extractHeight(s,t,k,m):
    
    '''
    Takes a shaped signal, peak rise slope time (sample #) and the shaping parameters and finds the height of the signal for histogramming.
    '''
    
    if (t+numpy.floor((k+m)/2))<len(s):
        peak = s[t+numpy.floor((k+m)/2)]
    else:
        peak = max(s)
        
    return peak

def peakFit(hist,peak_locs):
    
    '''
    Takes the energy histogram and the indices of the peaks, fits the peaks, returns the centroid and the FWHM in a list of tuples.
    '''
    def exp_func(x,a,b,c,m,y0):
        return a*numpy.exp((-5.5451*(x-b)**2.)/(2.*c**2))+m*x+y0
    
    def lin_func(x,a,b):
        return a*x+b
    
    peaks = []
    
    if len(peak_locs) > 0:
        for peak_loc in peak_locs:
            try:
                #baseline_info,garbage = opt.curve_fit(lin_func,numpy.arange(-10,10),numpy.array(hist[peak_loc-10:peak_loc+10]))
                #hist_s = numpy.array(hist[peak_loc-10:peak_loc+10]) - baseline_info[0]*numpy.arange(peak_loc-10,peak_loc+10)\
                #        - baseline_info[1]*numpy.ones(20)
                peak_info,cov = opt.curve_fit(exp_func,numpy.arange(-10,10),numpy.array(hist[peak_loc-10:peak_loc+10]))
                peak_info[1] = peak_info[1] + peak_loc
                peak_info = numpy.insert(peak_info,3,cov[2,2])
                peaks.append(peak_info)
            except:
                print 'Peak fitting failed for peak at ', peak_loc
                continue
    else:
        print 'Error: There were no peaks to fit.'

    return numpy.array(peaks)

def makeSpect(fname,ftype,k = 250,m = 100,pulseHeight=None,energies = [],numbins=2000):
    
    '''
    Takes a file name for data to be processed.  Processes and produces a spectrum.
    '''
    
    if pulseHeight == None:
        f = readFile(fname,ftype)
        M = 0
        pulseHeight = numpy.zeros(len(f))
        numPulses = len(f)
        zeropad = numpy.zeros(2*k+m)
        print 'Reading ', numPulses, ' traces from ', fname
        
        #for i in range(500):
        #    o = fitExp(subtractBaseline(f[i,:]))
        #    M = (M + o)/2
        start = time.time()
        
        i = 0
        numRead = 10000
        iterations = (numPulses-numPulses%numRead)/numRead + 1
        for n in range(iterations):
            n1 = n*10000
            if n >= iterations-1:
                n2 = numPulses
            else:
                n2 = n1 + 10000
            for trace in f[n1:n2,:]:
                trace_s = subtractBaseline(trace)
                s, pulseHeight[i] = trapezoidalFilter(k,m,4467,numpy.append(zeropad,trace_s))
                #t = findRise(trace_s)
                #pulseHeight[i] = extractHeight(s,1000,k,m)
                i += 1
            end = time.time()
        print 'Processing ',n2, ' samples (',n1, ' to ', n2, ') took: ', end-start
        
    #hist,low_range,binsize,extrapoints = stats.histogram(pulseHeight,numbins=max(pulseHeight))
    #plt.plot(hist[0])
    #plt.show()
    #def onclick(event):
    #    return event.xdata
    #cid = fig.mpl_connect('button_press_event',onclick)
    if len(energies) == 0:
        usr_in = ''
        while usr_in != 'done':
            usr_in = raw_input('Please enter energies of the expected peaks.  When finished, enter \'done\': ')
            try:
                energies.append(float(usr_in))
            except:
                if usr_in != 'done':
                    print 'Error: Please enter numbers. When finished, enter \'done\''
                continue
    
    hist,low_range,binsize,extrapoints = stats.histogram(pulseHeight,numbins=int(max(energies)+100))
    histogram = hist*(hist>1)
    histogram = [point for point in histogram if point>0]
    hist,low_range,binsize,extrapoints = stats.histogram(pulseHeight,numbins=int(max(energies)+numbins),defaultlimits=(0.,len(histogram)*binsize))
    tries = 0
    threshold = 10
    peak_indices = []
    while len(energies) != len(peak_indices):
        tries += 1
        thresh=threshold+tries*10
        peak_indices = findRise(hist,findPeak=True,thresh=thresh)
        if tries > 100:
            print 'Error: Could not match number of peaks to energies entered.'
            break
    
    peak_info = peakFit(hist,peak_indices)
    calib,cov = opt.curve_fit(lambda x,m,b: m*x+b, numpy.array([peak[1] for peak in peak_info])*binsize,numpy.array(energies))
    print calib
    bins = numpy.array([binsize*x*calib[0]+calib[1] for x in xrange(len(hist))])
    plt.plot(bins,hist)
    plt.show()
    peak_indices = findRise(hist,findPeak=True,thresh=thresh)
    peak_info = peakFit(hist,peak_indices)
    peak_info[:,1:2] *= binsize*calib[0]
    peak_info[:,1:2] += numpy.ones(peak_info[:,1:2].shape)*calib[1]
    peak_info[:,3] *= calib[0]
    print 'Peak #\tPeak Info'
    print '\tHeight (counts)\t\tCentroid (keV)\t\tFWHM (keV)\t\tsigma_FWHM (keV)'
    for peak_num,peak in enumerate(peak_info):
        print peak_num, list(peak)
    #plt.plot(hist[0])
    #plt.show()
    return hist, bins, peak_info, pulseHeight

def compareTrapParams():
    
    '''
    Script to test rise and gap times on data and determine optimal parameters.
    '''
    
    rise = 410
    gap = 60
    peak_array = numpy.ndarray([15,3,4])
    
    for i,j in enumerate(range(0,150,10)):
        hist,bins,peak_info,pulseHeight = makeSpect('../JB_AP_CO_60_AM_241_2000_samples_092413.h5','HDF5',rise,gap+j,energies=[59.5,1173,1332])
        peak_array[i,:,:] = peak_info
        fname = 'pulseHeight_' + str(rise) + '_' + str(gap+j) + '.txt'
        header = 'Peak #\t\t{Height} [counts]\t\{tCentroid\t\tFWHM\t\tsigma_FWHM} [keV]'+ '\n' + str([str(peak_num)+str(peak) for peak_num,peak in enumerate(peak_info)])
        print peak_num, peak
        numpy.savetxt(fname,pulseHeight,header=header)
            
    return peak_array