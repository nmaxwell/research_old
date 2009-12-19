
#include <mathlib/math/std_math.h>
#include <mathlib/math/kernels/FDkernels.h>
//#include <mathlib/math/hdaf/hdaf.h>
#include <mathlib/math/transforms/FFT.h>
//#include <mathlib/math/random/ml_random.h>
#include <mathlib/link.cpp>
#include <mathlib/non-lib_things.h>

#define mp_0 ((mp_real)(0.0))
#define mp_1 ((mp_real)(1.0))
#define mp_2 ((mp_real)(2.0))
#define mp_pi (atan(mp_1)*4.0)
#define mp_2pi (mp_pi*2.0)
#define mp_4pi2 (mp_2pi*mp_2pi)
#define mp_sqrt2 (sqrt(mp_2))
#define mp_iu ( complex<mp_real > (mp_0, mp_1 ) )


inline int periodic_mod(int const & k, int const & n)
{
    if (k>=0)
        return k%n;
    else
        return k+n; // do better;
        //return -k-div_up(-k,n);
}

void periodic_convolve ( int n, double *& in, double *& out, int M, double *& kernel, double scale )
{
    // convplution: periodic boundaries, (2*M+1 wide kernel), n points.
    
    if (!out) out = ml_alloc<double > (n);
    
    for (int k=0; k<n; k++)
    {
        double sum = 0.0;
        
        for (int j=-M; j<=M; j++)
            sum += in[periodic_mod(k-j,n)]*kernel[j];
        
        out[k] = sum*scale;
    }
}



double l2_error( double *& analytic, double *& numeric, int n, double a, double b)
{
    double sum1 = 0.0;	    
    for (int k=0; k<n; k++)
	sum1 += (numeric[k]-analytic[k])*(numeric[k]-analytic[k]);
    
    double sum2 = 0.0;
    for (int k=0; k<n; k++)
	sum2 += analytic[k]*analytic[k];
    
    if ( sum2 > 0 )
        return sqrt(sum1/sum2);
    else
        return sqrt(sum1);
}


double l2_norm( double *& data, int n, double a, double b)
{
    double sum = 0.0;
    for (int k=0; k<n; k++)
        sum += data[k]*data[k];
    
    return sqrt(sum/n);
}






