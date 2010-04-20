
from setup import *

kernel = lambda t,f: exp(-2.0j*pi*t*f)

P = propagator1()
P.set( grid, velocity=velocity_grid, damping=damping_grid,  expansion_order=5,   hdaf_m1=8, hdaf_m2=8, hdaf_gamma1=0.8, hdaf_gamma2=0.8 )


frequency = 47.110022310144111
final_time = 0.7
time_step = .001
every_nth_step = 5

FT_sum = grid.zeros()
u = grid.zeros()
v = grid.zeros()
time = 0
count = 0



while time <= final_time:
    
    print "step:", count, "\ttime:", time, "\twall time", wall_time()
    
    driving_1[i0,j0] = driving_function_time(time)
    driving_2[i0,j0] = driving_function_time(time+time_step)
    
    P( time_step, u,v, u,v, driving_1, driving_2 )
    
    K = grid.ones()*kernel(time,frequency)
    
    FT_sum += u*K*time_step
    
    
    if count%every_nth_step==0:
        write_png( u, dir+ "%04d.png"%count, center=mean(u), major_scale=std(u)*3, red_params=(0.33,0,0,1), green_params=(0.5,0,1,0), blue_params=(0.66,1,0,0,) )
    
    count += 1
    time += time_step




FT_sum_real = real(FT_sum)
FT_sum_imag = imag(FT_sum)

write_png( FT_sum_real, dir+ "FT_real.png", center=mean(FT_sum_real), major_scale=std(FT_sum_real)*3, red_params=(0.33,0,0,1), green_params=(0.5,0,1,0), blue_params=(0.66,1,0,0,) )
write_png( FT_sum_imag, dir+ "FT_imag.png", center=mean(FT_sum_imag), major_scale=std(FT_sum_imag)*3, red_params=(0.33,0,0,1), green_params=(0.5,0,1,0), blue_params=(0.66,1,0,0,) )

