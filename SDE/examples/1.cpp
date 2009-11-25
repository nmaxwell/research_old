
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







int main()
{
    std_setup();
    
    int n_runs = 1;
    int n_steps = 1E4;
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
        
        for ( int j=1; j<n_steps; j++ )
            W[run][j] = W[run][j-1] + dW[j-1]*sqrt(dt);
            //W[run][j] = dW[j];
    }
    
    for (int run=0; run<n_runs; run++)
    {
        sprintf(fname, "/workspace/output/temp/out_%d", run);
        output( W[run], n_steps, fname );
    }
    
    
    ml_free( W, n_runs );
    ml_free( dW );
    
    std_exit();
}



