
from dirichlet_laplace import *
from scipy import stats
import scipy.integrate


_4pi = 12.566370614359172
_2pi = 6.2831853071795865


def cylinder_integrand( z, phi, x ):
    y = array((R*cos(phi), R*sin(phi), z))
    return f(y)*R/(l2norm(x-y)*_4pi)

def cap1_integrand( r, phi, x ):
    y = array((r*cos(phi), r*sin(phi), 0))
    return f(y)*r/(l2norm(x-y)*_4pi)

def cap2_integrand( r, phi, x ):
    y = array((r*cos(phi), r*sin(phi), H))
    return  f(y)*r/(l2norm(x-y)*_4pi)


def potential(x):
    
    x = array(x)
    
    V=0.
    
    v,e = scipy.integrate.dblquad(cylinder_integrand, 0., _2pi, lambda x:0., lambda x:H, args=(x,) )
    V += v
    
    v,e = scipy.integrate.dblquad(cap1_integrand, 0., _2pi, lambda x:0., lambda x:R, args=(x,)  )
    V += v
    
    v,e = scipy.integrate.dblquad(cap2_integrand, 0., _2pi, lambda x:0., lambda x:R, args=(x,) )
    V += v
    
    return V




H = 3.
R = 1.

def f(x):
    z = x[2]
    
    if z <= 0.5:
        return -0.
    
    if z >= H-0.5:
        return 10.
    
    return 0.

def region(x):
    r = l2norm((x[0],x[1]))
    z = x[2]
    return r<=R and 0<=z and z<=H

def drift(t):
    return array((0.,0.,0.))


def check_solution():
    
    n_points = 20
    T = [ float(k)/(n_points-1) for k in range(n_points)]
    #X = [ array(( 0.7*cos(2.*_2pi*t), 0.7*sin(2.*_2pi*t), 0.8*t+0.1 )) for t in T ]
    X = [ array((0.,0., t*H )) for t in T ]
    
    V_quad = [ potential(x) for x in X ]
    
    n_samples = 200
    dt = 0.02
    
    V_sde = [ solution( x, f, n_samples, dt, drift, region) for x in X ]
    
    for k in range(n_points):
        print X[k][2], '\t',  V_quad[k], V_sde[k]


check_solution()























def check_solution_convergence_samplesize():
    
    f = lambda x: math.atan2(x[0],x[1])
    region=lambda x: l2norm(x)<=1.0
    
    drift=lambda t: array((0.7, 1.2 ))
    
    means = []
    sigms = []
    n_sampless = []
    
    for p in range(5,20):
        
        n_samples = int(1.5**p)
        print n_samples
        
        x = (.3, .4 )
        dt = 0.01
        
        mean,sigm = solution( x,f, n_samples ,dt,drift,region)
        
        means.append(mean)
        sigms.append(sigm)
        n_sampless.append(n_samples)
        
        print mean, sigm
    
    gradient, intercept, r_value, p_value, std_err =  stats.linregress([ log(n_samples) for n_samples in n_sampless ], [ log(sigm) for sigm in sigms ])
    
    print "\n"
    print "Result: ", means[-1], "+/-", sigms[-1]
    print "Convergence: error(N) ~ C*N^a, C=", exp(intercept), "a=", gradient
    



def check_solution_convergence_dt():
    
    f = lambda x: math.atan2(x[0],x[1])
    region=lambda x: l2norm(x)<=1.0
    
    drift=lambda t: array((0.0, 0.0 ))
    
    means = []
    sigms = []
    dts = []
    
    for p in range(1,10):
        
        n_samples = 500
        
        x = (.3, .4 )
        dt = 2**-p
        
        mean,sigm = solution( x,f, n_samples ,dt,drift,region)
        
        means.append(mean)
        sigms.append(sigm)
        dts.append(dt)
        
        print mean, sigm
    




    
    




