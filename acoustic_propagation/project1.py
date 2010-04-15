
from math import *
from numpy import *
import numpy
from waveprop import *

class Polygon2D:
    
    def __inint__(self, x_points, y_points )
        
        self.x_points = numpy.array(x_points)
        self.y_points = numpy.array(y_points)
        self.n_points = len(x_points)
        if len(self.x_points) != len(self.y_points)
            print "error: len(self.x_points) != len(self.y_points)"
    
    def inerior_test(x,y):
        
        i=j=0
        c = bool(False)
        j=self.n_points
        while j < self.n_points-1:
            
            if ((self.y_points[i]>y) != (self.y_points[j]>y)) :
                if (x < (self.x_points[j]-self.x_points[i]) * (y-self.y_points[i]) / (self.y_points[j]-self.y_points[i]) + x_points[i])
                    c = not c;
            i += 1
            j = i
        
        return c:



int Polygon2D::interior_test(const float & x, const float & y)
{
    // from http://www.ecse.rpi.edu/Homepages/wrf/Research/Short_Notes/pnpoly.html
    
    int i, j, c = 0;
    for (i = 0, j = n_points-1; i < n_points; j = i++)
        if ( ((y_points[i]>y) != (y_points[j]>y)) && (x < (x_points[j]-x_points[i]) * (y-y_points[i]) / (y_points[j]-y_points[i]) + x_points[i]) )
			c = !c;
    return c;
}

dir = "/workspace/output/acoustic_propagation/"


grid = grid2d( n1 = 150, n2 = 600, dx1=20.0, dx2=1.5)


velocity = grid.zeros();

mcn = 0
mcn1 = 0

for j_prime in range(grid.n2):
    for i_prime in range(grid.n1):
        
        i = grid.n1-i_prime-1
        j = grid.n2-j_prime-1
        
        j = j_prime
        i = i_prime
        
        if j <= 270:
            velocity[i,j] = 2000.0**2
        
        elif j >= 270 and j <= 307 and i >= 75 and i < 75+2*mcn+1:
            mcn = j-270
            velocity[i,j] = 4500.0**2
        elif j >= 307 and j <= 344:
            if i <= 75:
                velocity[i,j] = 4500**2
            elif i >= 75+2*mcn1+1:
                mcn1=j-307
                velocity[i,j] = 4500.0**2
            else:
                velocity[i,j] = 2000.0**2
        else:
            velocity[i,j] = 2000.0**2

D1 = Polygon2D([ 0.0,0.0,1520.0,1520.0,0.0], [460.5,516.0,516.0,405.0,460.5])
D2 = Polygon2D([ 1520.0,1520.0,3140.0,3020.0,1520.0 ], [ 405.0,460.5,516.0,460.5,405.0 ])



write_png_resample( velocity, dir+ "vlocity.png", grid_new=grid2d(n1=1500,n2=450, b1=grid.b1, b2=grid.b2), grid_old=grid, center = numpy.mean(velocity), major_scale=std(velocity), red_params=(0.33,0,0,1), green_params=(0.5,0,1,0), blue_params=(0.66,1,0,0,) )


def damping(x,y):
    
    i = x/grid.dx1
    j = y/grid.dx2
    
    return 10.0*( cosh(0.18*i)**-2 + cosh(0.18*(grid.n1-i))**-2 + cosh(0.18*(grid.n2-j))**-2 + cosh(0.18*j)**-2 )

damping = grid.evaluate( damping )

write_png_resample( damping, dir+ "damping.png", grid_new=grid2d(n1=1500,n2=450, b1=grid.b1, b2=grid.b2), grid_old=grid, center = numpy.mean(damping), major_scale=std(damping ), red_params=(0.33,0,0,1), green_params=(0.5,0,1,0), blue_params=(0.66,1,0,0,) )




def driving_function_time_general(t, sigma, gamma, tau ):
    return -sqrt(2.0/pi)*sigma*gamma*(sigma-2.0*sigma*gamma*(sigma*t-tau)**2)*exp(-gamma*(sigma*t-tau)**2)

driving_function_time = lambda t: driving_function_time_general(t=t, sigma=1.5*20, gamma=8, tau=1 )

"""
T = 0.7*array([ float(i)/1000.0 for i in range(1000) ])
F = [ driving_function_time(t) for t in T ]

import pylab as p
p.plot( T, F)
p.show()
"""

time_step = 0.001
final_time = 0.7

count = 0
driving_1 = grid.zeros()
driving_2 = grid.zeros()

u = grid.zeros()
v = grid.zeros()

u[75,0] = 1.0

P = propagator1()
P.set( grid, expansion_order=10, velocity=velocity, damping=damping )


every_nth_step = 50
max_count = inf


time=0

while time<final_time and count<max_count   :
    
    print "time=", time
    
    driving_1[75,50] = driving_function_time(time)
    driving_2[75,50] = driving_function_time(time+time_step)
    
    P( time_step, u,v, u,v, driving_1, driving_2 )
    
    if count%every_nth_step==0:
        write_png_resample( u, dir+ "%04d.png"%count, grid_new=grid2d(n1=1500,n2=450, b1=grid.b1, b2=grid.b2), grid_old=grid, center=numpy.mean(u), major_scale=numpy.std(u), red_params=(0.33,0,0,1), green_params=(0.5,0,1,0), blue_params=(0.66,1,0,0,) )
    count += 1
    
    time += time_step

