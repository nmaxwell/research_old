
#define ML_NOT_USE_FFTW_MALLOC

#include <mathlib/math/random/ml_random.h>
#include <mathlib/math/SDE/random_walk.h>
#include <mathlib/math/grids/grid1D.h>
#include <mathlib/math/grids/grid2D.h>
//#include <mathlib/math/grids/extra/plots.cpp>

#include <mathlib/link.cpp>
#include <mathlib/non-lib_things.h>

//#include <boost/math/special_functions/erf.hpp>




template< class T >
void output( T * data, int n, const char * fname, const char * delim1="\t", const char * delim2="\n" )
{
    ofstream out;
    out.open(fname, fstream::out);
    assert(out.good() && out.is_open());
    
    for ( int k=0; k<n; k++ )
        out << data[k] << delim2;
    
    out.close();
}

template< class T1, class T2 >
void output( T1 * data1, T2 * data2, int n, const char * fname, const char * delim1="\t", const char * delim2="\n" )
{
    ofstream out;
    out.open(fname, fstream::out);
    assert(out.good() && out.is_open());
    
    for ( int k=0; k<n; k++ )
        out << data1[k] << delim1 << data2[k] << delim2;
    
    out.close();
}



/*

import math
import pylab as p
import sys
sys.path.append("/workspace/mathlib/tools/python/")
from read_file import *

for n in range(1,10):
    p.plot( read_file( "out_" + str(n) ) )

p.plot ( read_file( "U_mean" )  )
p.show()

Y = read_file( "U_mean" )
X = [ float(k) / len(Y) for k in range(len(Y))]
Y = [ math.exp(x) for x in X ]
p.plot(X,Y,X, read_file( "U_mean" )  )
p.show()


*/



int main()
{
    std_setup();
    
    int n_runs = 20;
    int n_steps = 1E6;
    double **W=0, *dW=0;
    
    W = ml_alloc<double > (n_runs, n_steps );
    dW = ml_alloc<double > ( n_steps );
    
    ml_random rng;
    
    double stop_time = 1.0;
    double dt = stop_time/n_steps;
    
    for (int run=0; run<n_runs; run++)
    {
        rng.std_normal_rv( dW, n_steps);
        
        W[run][0] = 0;
        
        for ( int step=0; step<n_steps; step++ )
            W[run][step] = W[run][step-1] + dW[step-1]*sqrt(dt);
        
        
        double ito_sum = 0.0;
        
        for (int step=0; step<n_steps-1; step++)
            ito_sum += W[run][step]*(W[run][step+1]-W[run][step]);   
        
        cout << "ito: " << ito_sum << "\t" << 0.5*( W[run][n_steps-1]*W[run][n_steps-1] - stop_time )  << "\t" << fabs( ito_sum - 0.5*(W[run][n_steps-1]*W[run][n_steps-1] - stop_time) ) <<endl;
        
        
        double strat_sum = 0.0;
        
        for (int step=0; step<n_steps-1; step++)
            strat_sum += (W[run][step] + W[run][step+1] )*(W[run][step+1]-W[run][step])/2;   
        
        cout << "strat: " << strat_sum << "\t" << 0.5*( W[run][n_steps-1]*W[run][n_steps-1] )  << "\t" << fabs( strat_sum - 0.5*(W[run][n_steps-1]*W[run][n_steps-1] ) ) <<endl;
        
    }
    
    ml_free( W, n_runs );
    ml_free( dW );
    
    std_exit();
}



