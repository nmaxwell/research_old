


dir = "/workspace/output/acoustic_propagation/"

import time as time_module
from math import *
from numpy import *
import numpy
from waveprop import *


wall_time_0 = time_module.clock()
wall_time = lambda : time_module.clock()-wall_time_0

class polygon2D:
    
    def __init__(self, vertices ):
        
        self.x_points = [ float(v[0]) for v in vertices ]
        self.y_points = [ float(v[1]) for v in vertices ]
        self.n_points = len(self.x_points)
    
    def interior_test(self, x,y):
        
        i=j=0
        c = bool(False)
        j=self.n_points-1
        while i < self.n_points:
            if ( self.y_points[i] > y ) != ( self.y_points[j] > y ) :
                if x < (self.x_points[j]-self.x_points[i]) * (y-self.y_points[i]) / ( self.y_points[j]-self.y_points[i] ) + self.x_points[i]:
                    c = not c
            j = i
            i += 1
        
        return c
    
    def __contains__(self, x):
        return self.interior_test(x[0],x[1])

class union:
    def __init__(self, set1=None, set2=None, set3=None, set4=None  ):
        self.sets=[]
        for E in [set1,set2,set3,set4]:
            if E != None:
                self.sets.append(E)
    
    def __contains__(self, x):
        for E in self.sets:
            if x in E:
                return True
        return False

class complement:
    def __init__(self, complement_set  ):
        self.complement_set = complement_set
    
    def __contains__(self, x):
        return not ( x in self.complement_set )

class simpleFunction:
    
    def __init__(self, sets_values ):
        self.sets = [ x[0] for x in sets_values ]
        self.values = [ x[1] for x in sets_values ]
        if len(self.sets) != len(self.values):
                print "error in simple_function: len(self.sets) != len(self.values):"
    
    def __call__(self, x ):
        sum = (self.values[0])*0.0
        for k,E in enumerate(self.sets):
            if x in E:
                sum += self.values[k]
        return sum





grid = grid2d( a1 = -500, b1 = 3500, a2=-300, b2=1200, n1=1024, n2=384)

region_1 = polygon2D([ (0,375), (1500,375), (1500,438), (0,438) ])
region_2 = polygon2D([ (1500,438), (3000,373), (3000,438), (1500,495) ])

region = union( region_1, region_2 )

velocity = simpleFunction([ (region,4500.0**2), (complement(region),2000.0**2) ])

def damping(x,y):
    scale = 50.0
    return 10000*( cosh(abs(x-grid.a1)/scale)**-2 + cosh(abs(x-grid.b1)/scale)**-2 + cosh(abs(y-grid.a2)/scale)**-2 + cosh(abs(y-grid.b2)/scale)**-2  )

velocity_grid = grid.evaluate( lambda x,y: velocity((x,y)) )
damping_grid = grid.evaluate(damping )

write_png( velocity_grid, dir+ "vlocity.png", center=mean(velocity_grid), major_scale=std(velocity_grid), red_params=(0.33,0,0,1), green_params=(0.5,0,1,0), blue_params=(0.66,1,0,0,) )
write_png( damping_grid, dir+ "damping.png", major_scale=1.0, red_params=(0.33,0,0,1), green_params=(0.5,0,1,0), blue_params=(0.66,1,0,0,) )


def driving_function_time_general(t, sigma, gamma, tau ):
    return -sqrt(2.0/pi)*sigma*gamma*(sigma-2.0*sigma*gamma*(sigma*t-tau)**2)*exp(-gamma*(sigma*t-tau)**2)

driving_function_time = lambda t: driving_function_time_general(t=t, sigma=1.5*20, gamma=8, tau=1 )

i0 = grid.index1( 1500 )
j0 = grid.index2( 800 )

"""
T = 0.7*array([ float(i)/1000.0 for i in range(1000) ])
F = [ driving_function_time(t) for t in T ]

import pylab as p
p.plot( T, F)
p.show()
"""

print i0,j0

time_step = 0.001
final_time = 0.7
every_nth_step = 10
max_count = inf

count = 0
driving_1 = grid.zeros()
driving_2 = grid.zeros()
u = grid.zeros()
v = grid.zeros()
u[i0,j0] = 0.0

P = propagator1()
P.set( grid, velocity=velocity_grid, damping=damping_grid,  expansion_order=5,   hdaf_m1=8, hdaf_m2=8, hdaf_gamma1=0.8, hdaf_gamma2=0.8 )


if __name__ == "__main__":
    
    time=0

    while time<=final_time and count<max_count:
        
        print "step:", count, "\ttime:", time, "\twall time", wall_time()
        
        driving_1[i0,j0] = driving_function_time(time)
        driving_2[i0,j0] = driving_function_time(time+time_step)
        
        P( time_step, u,v, u,v, driving_1, driving_2 )
        
        if count%every_nth_step==0:
            write_png( u, dir+ "%04d.png"%count, center=mean(u), major_scale=std(u)*3, red_params=(0.33,0,0,1), green_params=(0.5,0,1,0), blue_params=(0.66,1,0,0,) )
            #write_png_resample( u, dir+ "%04d.png"%count, grid_new=grid2d(n1=1500,n2=450, b1=grid.b1, b2=grid.b2), grid_old=grid, center=numpy.mean(u), major_scale=numpy.std(u), red_params=(0.33,0,0,1), green_params=(0.5,0,1,0), blue_params=(0.66,1,0,0,) )
        
        count += 1    
        time += time_step



