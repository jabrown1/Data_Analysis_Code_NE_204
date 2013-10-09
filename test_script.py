import tables
import sys
import numpy as np;
import time;
from matplotlib.pyplot import *

def create_avg_filter(pl,p,g,cf):
    a=np.zeros((pl,pl));
    for n in range(0,pl):
        if (n<=p):
            for j in range(0,n):
                a[j,n]=p;
        elif (n<p+g):
            for j in range(0,n):
                a[j,n]=0;
            for k in range(0,n-p):
                a[k,n]=p;

        elif (n<2*p+g):
            for j in range(n-p,n):
                a[j,n]=p;
            for k in range(n-p-g,n-p):
                a[k,n]=0;
            for l in range(0,n-p-g):
                a[l,n]=-p;
        else:
            v=p;
            for j in range(n-p,n):
                a[j,n]=cf*p+cf*v;
                v-=1;
            for k in range(n-p-g,n-p):
                a[k,n]=cf*p;
            v=-p;
            for l in range(n-(2*p+g),n-(p+g)):
                a[l,n]=cf*v;
                v+=1;
    return a;

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
            v=-(g+cf)
            for l in range(0,n-p-g):
                a[l,n]=-l;
                v+=1;
        else:
            v=p;
            for j in range(n-p,n):
                a[j,n]=cf+v;
                v-=1;
            for k in range(n-p-g,n-p):
                a[k,n]=p;
            v=-(g+cf);
            for l in range(n-(2*p+g),n-p-g):
                a[l,n]=v;
                v+=1;
    return a;

def load_data():
	hf = tables.openFile("/Users/bnrc/Desktop/NE_204_data/JB_AP_CO_60_2000_samples_092413.h5");
	return hf;

def read_n_m_lines(f,n,m):
	
	g = f.root.RawData.read(n,m);

	for i in range(0,m-n):
		g[i,:]-=np.average(g[i,0:800]);
                  
	return g;
 
                  
def vector_product(x,y):
    start=time.time();
    z = np.dot(x,y);
    end=time.time();
    print end-start;
    return z;


		

