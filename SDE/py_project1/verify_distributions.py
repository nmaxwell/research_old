
from dirichlet_laplace import *
from scipy import stats


class interval:
    "interval (a,b] of the real numbers."
    
    def __init__(self, a,b ):
        self.a = a
        self.b = b
    
    def __contains__(self, x):
        return self.a < x and x <= self.b



def check_distribution_nball():
    
    dim =3
    n_measure_sets = 2**dim
    measure_sets = [ (k,) for k in range(n_measure_sets) ]
    
    region = lambda x: l2norm(x)<=1.0
    
    a = array([ random.normalvariate(0.5, .2 ) for k in range(dim) ])
    drift = lambda t: array([ random.normalvariate(0.5, .2 ) for k in range(dim) ])
    f = lambda x : sum([(2.0**k)*a for k,a in enumerate([xi>=0 for xi in x])])
    
    dt = 0.03
    x = zeros(dim)
    n_samples = 2000
    
    distribution, distribution_nocom = hitting_value_distribuion( n_samples, measure_sets, x, f, dt, drift, region )
    
    for k in distribution:
        print k
    
    print 'sum:', sum(distribution)
    
    print '\n', mean(distribution), var(distribution), sig_m(distribution)
    
    import pylab as p
    print '\n\n'
    
    left = range(len(distribution))
    
    p.bar(left, distribution, 1 )
    p.plot(left, distribution_nocom, 'ro' )
    p.show()
    
    

def check_distribution():
    
    f = lambda x: math.atan2(x[0],x[1])
    region=lambda x: l2norm(x)<=1.0
    
    drift=lambda t: array((0.7, 0.5 ))
    
    n_samples = 10000
    
    n_measure_sets = 16
    measure_sets = [ interval( 2.*pi*k/n_measure_sets-pi, 2.*pi*(k+1)/n_measure_sets-pi ) for k in range(n_measure_sets) ]
    
    
    x = (0.3, 0. )
    dt = 0.02
    
    distribution, distribution_nocom = hitting_value_distribuion( n_samples, measure_sets, x, f, dt, drift, region )
    
    for k in distribution:
        print k
    
    print 'sum:', sum(distribution)
    
    print '\n', mean(distribution), var(distribution), sig_m(distribution)
    
    
    
    import pylab as p
    print '\n\n'
    
    left = [  (2.*k/n_measure_sets-1)*pi  for k in range(n_measure_sets) ]
    
    p.bar(left, distribution, 2.*pi/n_measure_sets )
    p.plot(left, distribution_nocom, 'ro' )
    p.show()





    
check_distribution()




