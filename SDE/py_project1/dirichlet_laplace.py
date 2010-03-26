
import numpy
from numpy import *
import random
import math

l2norm = numpy.linalg.norm


def hitting_position( x=(0.,0.), dt=0.1, drift=lambda t:0.0, region=lambda x: l2norm(x)<=1.0 ):
    
    x=array(x)
    sqrt_dt = sqrt(dt)
    zero = zeros(len(x))
    dB_function = lambda : array([ sqrt_dt*random.normalvariate(0.0, 1.0 ) for k in zero ])
    drift_function = drift  
    
    t = 0.0
    X = x
    logM_T = 0.0
    
    while region(X):
        
        dB = dB_function()
        drift = drift_function(t)
        
        dX = dB + drift*dt
        logM_T += -dot(drift,dB) -0.5*dot(drift,drift)*dt
        
        X += dX
        t += dt
    
    M_T = exp(logM_T)
    
    return X, M_T


def hitting_value( x=(0.,0.), f=lambda x:l2norm(x), dt=0.1, drift=lambda t:0.0, region=lambda x: l2norm(x)<=1.0 ):
    
    X, M_T = hitting_position(x,dt,drift,region)
    
    return f(X), M_T


def hitting_value_distribuion( n_samples = 100, measure_sets=[[]], x=(0.,0.), f=lambda x:l2norm(x), dt=0.1, drift=lambda t:0.0, region=lambda x: l2norm(x)<=1.0 ):
    
    distribution = [0. for s in measure_sets ]
    distribution_nocom = [0. for s in measure_sets ]
    
    for k in range(n_samples):
        
        fX, M_T = hitting_value( x, f, dt, drift, region )
        
        for k,E in enumerate(measure_sets):
            if fX in E:
                distribution[k] += M_T/n_samples
                distribution_nocom[k] += 1.0/n_samples
        
    
    return distribution, distribution_nocom


def sig_m(X):
    xm = mean(X)
    return sqrt(sum([ (xi - xm)**2 for xi in X ])/(len(X)*(len(X)-1)))


def solution( x=(0.,0.) , f=lambda x:l2norm(x), n_samples = 100, dt=0.1, drift=lambda t:0.0, region=lambda x: l2norm(x)<=1.0 ):
    
    results = []
    for k in range(n_samples):
        fX, M_T = hitting_value( x,f,dt,drift,region)
        results.append(fX*M_T)
    
    return mean(results), sig_m(results)















