
from pylab import *
sys.path.append("/workspace/mathlib/tools/python/")
from read_file import *

dir = "/workspace/output/SDE/test/"



X = read_file( dir + "X" )
time = read_file( dir + "time" )
Xm = read_file( dir + "X_mean" )


for k in range(20):
    plot( time, X[k])
 
plot( time, Xm)
show()




