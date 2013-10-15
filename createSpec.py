
from numba import autojit,codegen
import numpy
import matplotlib.pyplot as plt
import logging
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
        import numpy as np
        f = np.loadtxt(fname)
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

    s = np.zeros_like(v)
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
    n = numpy.argmax(v)
    length = numpy.size(v)
    
    for n in range(2*k+m,length):
        d[n] = v[n] - v[n-k] - v[n-m] + v[n-k-m]
        p[n] = p[n-1] + d[n]
        r[n] = p[n] + M*d[n]
        s[n] = s[n-1] + r[n]
        
    return s

@autojit
def subtractBaseline(v):
    
    '''
    Takes in signal v and subtracts baseline by taking the average of samples before a sharp rise.
    '''
    for i in range(numpy.size(v,axis=1)):
        v[:,i] -= numpy.mean(v[:,:100],axis=1)
        
    return v

@autojit
def fitExp(v):
    n = numpy.argmax(v)+25
    xdata = numpy.arange(np.size(s[n:]))
    ydata = numpy.array(s[n:])
    fit,fit2 = opt.curve_fit(lambda x,a,b: a*np.exp(-x/b),xdata,ydata)
    return int(np.round(fit[1]))

@autojit
def findRise(v):
    
    for j in range(len(v)-1):
        Dv[j] = v[j+1]-v[j]
        
    return np.argmax(Df)

def extractHeight(s,t,k,m):
    
    '''
    Takes a shaped signal, peak rise slope time (sample #) and the shaping parameters and finds the height of the signal for histogramming.
    '''
    
    if abs(len(s)/2-t)<50:
        peak = s[t+k+np.floor(m/2)]
    else:
        peak = max(s)
        
    return s

def makeSpect(fname,ftype):
    
    '''
    Takes a file name for data to be processed.  Processes and produces a spectrum.
    '''
    
    f = readFile(fname,ftype)
    k = 250
    m = 100
    pulseHeight = np.zeros(len(f[:,1]))
    
    for i in range(100):
        o = fitExp(f[1,:])
        M = (M + o)/2
    
    for i,trace in enumerate(f):
        trace_s = subtractBaseline(trace)
        s = trapezoidalFilter(k,m,M,trace_s)
        t = findRise(trace_s)
        pulseHeight[i] = extractHeight(s,t,k,m)
    
    hist = np.histogram(pulseHeight,bins=1000)
    plt.plot(hist)
    plt.show()