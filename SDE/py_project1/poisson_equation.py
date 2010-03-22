
import numpy
from numpy import *
import random
import math

l2norm = numpy.linalg.norm


def hitting_position( x=(0.,0.), dt=0.1, drift=lambda t:0.0, region=lambda x: l2norm(x)<=1.0 ):
    
    x = array(x)
    dim = len(x)
    zero = zeros(dim)
    sqrt_dt = sqrt(dt)
    
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
    
    return (X, M_T)


def hitting_value( x=(0.,0.), f=lambda x:l2norm(x), dt=0.1, drift=lambda t:0.0, region=lambda x: l2norm(x)<=1.0 ):
    
    (X, M_T) = hitting_position(x,dt,drift,region)
    
    return (f(X),M_T)


def poisson_solution( x=(0.,0.) , f=lambda x:l2norm(x), n_samples = 100, dt=0.1, drift=lambda t:0.0, region=lambda x: l2norm(x)<=1.0 ):
    
    sum = 0.0
    sumM_T = 0.0
    for k in range(n_samples):
        (fX, M_T) = hitting_value( x,f,dt,drift,region)
        sum += fX*M_T
        sumM_T += M_T
    u = sum/sumM_T
    
    return u




if __name__ == "__main__":
    
    f = lambda x: math. atan2(x[0],x[1])
    region=lambda x: l2norm(x)<=1.0
    
    drift=lambda t: array((0.0, 0.0 ))
    
    n_samples = 100
    dt = 0.03
    
    u = lambda x: poisson_solution( x,f,n_samples,dt,drift,region )
    
    for p in range(1,10):
        
        n_samples = 2**p
        
        n_runs = 20
        X = [ u((0., 0. ))  for k in range(n_runs) ]
        
        print n_samples, mean(X), var(X)
    
    
    
    
    













