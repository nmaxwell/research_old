#include <mathlib/math/std_math.h>

#include <mathlib/math/SDE/BM.h>

#include <mathlib/link.cpp>
#include <mathlib/non-lib_things.h>

#include <boost/lexical_cast.hpp>

/*
double drift( double t )
{
    //return 0.0;
    
	return (0.7*t+cos(t*_2pi))*1.0;
}*/

int main( int argc, char *argv[] )
{
    //for (int k=1; k<argc; k++ )
    //    cout << argv[k] << endl;
    
    
	std_setup();
    
    // setup:
    
    int n_steps = 100;
    int n_runs = 100;
    double X0 = 0.0;
    double stop_time = 2.0;
    double * drift=0;
    char dir[] = "/workspace/output/SDE/test";
    
    {
        
        // n_runs 
        
        if(argc>=2)
        try {
            n_runs = boost::lexical_cast< int >( argv[1] );
        }
        catch (const boost::bad_lexical_cast &) {
            
        }
        
        // n_steps 
        
        if(argc>=3)
        try {
            n_steps = boost::lexical_cast< int >( argv[2] );
        }
        catch (const boost::bad_lexical_cast &) {
            
        }
        
        // X0 
        
        if(argc>=4)
        try {
            X0 = boost::lexical_cast< double >( argv[3] );
        }
        catch (const boost::bad_lexical_cast &) {
            
        }
        
        // stop_time
        
        if(argc>=4)
        try {
            stop_time = boost::lexical_cast< double >( argv[4] );
        }
        catch (const boost::bad_lexical_cast &) {
            
        }
        
        // drift
        
        if(argc>=6)
        {
            ifstream in;
            sprintf(fname,"%s/%s", dir, argv[5]);
            in.open(fname, ios::in);
            if ( in.good() && in.is_open() );
            {
                drift = ml_alloc<double> (n_steps);
                for (int k=0; k<n_steps; k++)
                    drift[k] = 0.0;
                char str[200];
                
                int k=0;
                while ( in>>str and k<n_steps )
                {
                    try {
                        drift[k] = boost::lexical_cast< double >( str );
                    } catch (const boost::bad_lexical_cast &) {
                        }
                    k++;
                }
                
                in.close();
            }
        }
    }
    
    
    double dt = stop_time/n_steps;
    
    double * time = ml_alloc<double > ( n_steps );
    for (int k=0; k<n_steps; k++ )
        time[k] = dt*k;
    
    double **W=0;
    double **X = ml_alloc<double> ( n_runs, n_steps );
    double **M = ml_alloc<double> ( n_runs, n_steps );
    
    if (drift == 0)
    {
        drift = ml_alloc<double> (n_steps);
        for (int k=0; k<n_steps; k++)
            drift[k] = 0.0;
    }
    
    {
        double sum =0;
        for (int k=0; k<n_steps; k++)
            sum += fabs(drift[k]);
        //cout << "sum drift: " << sum << endl;
    }
    
    // generate the BMs:
    
    gen_BM( dt, n_steps, W, n_runs, BM_mode_2 );
    
    
    // integrate the SDE
    
    for ( int run=0; run<n_runs; run++ )
    {
        X[run][0] = X0;
			
        for ( int j=1; j<n_steps; j++ )
			X[run][j] = X[run][j-1] + ( drift[j] + drift[j-1] )*dt/2 + (W[run][j]-W[run][j-1]);
    }
    
    // compute the Girsanov change of measure at every t
    
    for ( int run=0; run<n_runs; run++ )
    for ( int k=0; k<n_steps-1; k++ )
    {
        double sum1 = 0;
        for (int j=0; j<=k; j++)
            sum1 += drift[j]*drift[j]*dt;
        
        double sum2 = 0;
        for (int j=0; j<=k; j++)
            sum2 += drift[j]*(W[run][j+1]-W[run][j]);
        
        M[run][k] = exp(-sum2 -0.5*sum1);
    }
    
    sprintf(fname, "%s/X", dir );
    output( X, max(n_runs,min(100,n_runs)), n_steps, "/workspace/output/SDE/test/X" );
    sprintf(fname, "%s/M", dir );
    output( M, max(n_runs,min(100,n_runs)), n_steps, "/workspace/output/SDE/test/M" );
    
    ml_free( X, n_runs );
    ml_free( W, n_runs );
    ml_free( M, n_runs );
    
    std_exit();
}




