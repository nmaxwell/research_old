
from pylab import *
sys.path.append("/workspace/mathlib/tools/python/")
from read_file import *

dir = "/workspace/output/SDE/test/"



#X = read_file( dir + "X" )
time = read_file( dir + "time" )
Xm = read_file( dir + "X_mean" )

"""
for k in range(5):
    plot( time, X[k])

plot( time, Xm)
show()

"""

plot( time, Xm)
show()



"""
M = read_file( dir + "M" )
for k in range(20):
    plot( time, M[k])

show()




for k in range(20):
    Y  = [ M[k][j]*X[k][j] for j,x in enumerate(X[k]) ]
    plot( time, Y )

show()
"""


