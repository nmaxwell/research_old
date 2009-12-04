



#include <mathlib/math/std_math.h>
#include <mathlib/math/kernels/FDkernels.h>
#include <mathlib/math/hdaf/hdaf.h>
#include <mathlib/math/transforms/FFT.h>
#include <mathlib/math/random/ml_random.h>
#include <mathlib/link.cpp>
#include <mathlib/non-lib_things.h>



int main()
{
    
    int n = 16;
    double a = 0.0;
    double b = 1.0;
    double h = (b-a)/n;
    
    double * ker = ml_alloc<double > ( n );
    complex<double> * ker_ft = ml_alloc< complex<double> > ( n/2 );
    
    for ( int k=0; k<n/2; k++ )
        ker_ft[k] = 0;
    
    for ( int k=0; k<n/2; k++ )
        ker_ft[k] = 0.5* ( - (ml_4pi2)*( k*k ) );
        
        //ker_ft[k] = -(ml_4pi2*k*k)/((b-a)*(b-a));
    
    IFFT( ker_ft, ker, n );
    
    output( ker, n, "/workspace/output/temp/A" );
    
    for ( int k=0; k<n; k++ )
        ker[k] = cos( ((b-a)*((double)k)/(n)) * 4.0*_2pi/(b-a) ) * 16.0*(ml_4pi2)*(-1) ;
    
    output( ker, n, "/workspace/output/temp/B" );
}




