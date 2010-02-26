
from run_sde import *
import matplotlib.pyplot as plt

"""

print 'looping over n_steps, no drift: \n'

for p in range(5,10+1):
    
    n_steps=2**p
    n_runs=500
    
    
    print 'n_steps: ', 2**p, '\nn_runs: ', 500
    
    results = run_sde( n_runs = 500, n_steps=n_steps, X0=0.0, stop_time=1.0, drift=lambda t: 0.0 )
    
    time = results['time']
    w_slice = results['X']
    t_slice = w_slice.transpose()
    
    n_steps = len(t_slice)
    n_runs = len(w_slice)
    
    X_var = var(w_slice,0)
    X_mean = mean(w_slice,0)
    
    print "rms(X_mean): ", rms(X_mean)
    print "rms(X_var-time): ", rms(X_var-time)
    print '\n'



print 'looping over n_runs, no drift: \n'

for p in range(2,10+1):
    
    n_runs=2**p
    n_steps=512
    
    
    print 'n_steps: ', n_steps, '\nn_runs: ', n_runs
    
    results = run_sde( n_runs = 500, n_steps=n_steps, X0=0.0, stop_time=1.0, drift=lambda t: 0.0 )
    
    time = results['time']
    w_slice = results['X']
    t_slice = w_slice.transpose()
    
    n_steps = len(t_slice)
    n_runs = len(w_slice)
    
    X_var = var(w_slice,0)
    X_mean = mean(w_slice,0)
    
    print "rms(X_mean): ", rms(X_mean)
    print "rms(X_var-time): ", rms(X_var-time)
    print '\n'

"""


print 'looping over n_runs, no drift = t**2: \n'

for p in range(2,10+1):
    
    n_runs=2**p
    n_steps=512
    
    drift=lambda t: t**2
    
    print 'n_steps: ', n_steps, '\nn_runs: ', n_runs
    
    results = run_sde( n_runs = 500, n_steps=n_steps, X0=0.0, stop_time=2.0, drift=drift )
    
    time = results['time']
    M = results['M']
    X = results['X']
    
    X2= X*X
    
    PM = mean(M,0)
    PX = mean(X,0)
    QX = mean(X*M,0)
    QX2 = mean(X*X*M,0)
    VX = QX2-QX**2
    
    plt.plot(time, PX, time, [drift(t) for t in time])
    plt.savefig( dir + "out/" + "PX_" + '%02d' % p )
    plt.clf()
    plt.plot(time, QX)
    plt.savefig( dir + "out/" + "QX_" + '%02d' % p )
    plt.clf()
    plt.plot(time, QX2)
    plt.savefig( dir + "out/" + "QX2_" + '%02d' % p )
    plt.clf()
    plt.plot(time, VX)
    plt.savefig( dir + "out/" + "VX_" + '%02d' % p )
    plt.clf()
    plt.plot(time, PM)
    plt.savefig( dir + "out/" + "PM_" + '%02d' % p )
    plt.clf()
    
    print "rms(X_mean): ", rms(QX)
    print "rms(X_var-time): ", rms(VX-time)
    print '\n'






