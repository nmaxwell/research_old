
#include <mathlib/math/std_math.h>


#include <mathlib/math/SDE/BM.h>





#include <mathlib/link.cpp>
#include <mathlib/non-lib_things.h>




double drift( double t )
{
	return 0.7*t+cos(t*_2pi);
	
	
}

double volatility( double X, double t )
{
	return 1.0;
	
}



int main()
{
	std_setup();
    
    int n_steps = 500;
    int n_runs = 20000;
    
    double X0 = 0.0;
    double stop_time = 2.0;
    
    double dt = stop_time/n_steps;
    
    double **W=0, *time=0, *X_true=0;
    double **X = ml_alloc<double> ( n_runs, n_steps );
    double * X_mean = ml_alloc<double> (n_steps);
    
    double **M = ml_alloc<double> ( n_runs, n_steps );
    
    gen_BM( dt, n_steps, W, n_runs, BM_mode_2 );
    
    time = ml_alloc<double > ( n_steps );
    for (int k=0; k<n_steps; k++ )
        time[k] = dt*k;
     
    
    for ( int run=0; run<n_runs; run++ )
    {
		M[run][0] = 1.0;
		
        for ( int j=1; j<n_steps; j++ )
        {
            M[run][j] = M[run][j-1]*exp(
					-0.5*dt*( drift(dt*j)*drift(dt*j) + drift(dt*j-dt)*drift(dt*j-dt) )
					-0.5* drift(dt*j)*(W[run][j]-W[run][j-1])	);
        }
    }
    
    for ( int run=0; run<n_runs; run++ )
    {
        X[run][0] = X0;
			
        for ( int j=1; j<n_steps; j++ )
        {
         //   X[run][j] = X[run][j-1] + drift(X[run][j-1], dt*j)*dt + volatility(X[run][j-1], dt*j)*(W[run][j]-W[run][j-1]);
			
			X[run][j] = X[run][j-1] + ( drift(dt*j) + drift(dt*j-dt) )*dt/2 + volatility(X[run][j-1], dt*j)*(W[run][j]-W[run][j-1]);
        }
    }
    
    
    for ( int run=0; run<n_runs; run++ )
    for ( int step=0; step<n_steps; step++ )
        X[run][step] *= M[run][n_steps-1];
    
    
    for ( int run=0; run<n_runs; run++ )
    for ( int step=0; step<n_steps; step++ )
        X_mean[step] += (M[run][step])*X[run][step]/n_runs;
    
    
    
    output( X_mean, n_steps, "/workspace/output/SDE/test/X_mean" );    
    output( X, max(n_runs,100), n_steps, "/workspace/output/SDE/test/X" );
    output( time, n_steps, "/workspace/output/SDE/test/time" );
    output( M, max(n_runs,100), n_steps, "/workspace/output/SDE/test/M" );
    
    
    ml_free( X, n_runs );
    ml_free( X_mean );
    ml_free( W, n_runs );
    ml_free( time );
    ml_free( X_mean );
    
    std_exit();
}




