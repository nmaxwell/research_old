
from pylab import *
sys.path.append("/workspace/mathlib/tools/python/")
from read_file import *

dir = "/workspace/output/SDE/test/"



time = read_file( dir + "time" )
Xm = read_file( dir + "X_mean" )

Y = [ 0.5*t**2 for t in time ]

Y = [ Xm[k]-Y[k] for k,x in enumerate(Xm) ]

plot( time, Y)
show()

