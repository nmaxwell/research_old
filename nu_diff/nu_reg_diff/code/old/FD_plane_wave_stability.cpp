
#define INCLUDE_FD_D1_LR
//#define INCLUDE_FD_D2_LR

#include <mathlib/math/std_math.h>
#include <mathlib/math/kernels/FDkernels.h>
#include <mathlib/math/hdaf/hdaf.h>
#include <mathlib/math/transforms/FFT.h>
#include <mathlib/math/random/ml_random.h>
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

void hdaf_diff( double a, double b, int n, int m, double gamma, double *& in, double *& out )
{
    double h = (b-a)/n;
    
    double sigma = h*sqrt(2*m+1)/(pi*gamma);
    
    ml_poly<mp_real > P;
    make_hdaf_ml_poly(P,m );    
    differentiate_hdaf_poly( P, 1 );
    
    int w = (int)ceil(hdaf_truncate_point (1E-15, 5E-15, 8,sigma, 1 )/h);
    
    double *ker = ml_alloc<double > (2*w+1);
    
    mp_real H = h;
    mp_real s = H*H/((mp_2*sigma)*sigma);
    
    for (int j=-w; j<=w; j++ )
	ker[j+w] = dble( exp(-s*(j*j)) *P(dble((H*j)/(mp_sqrt2*sigma))) /((mp_2*sigma)*sigma) );
    
    double * p = &(ker[w]);
    periodic_convolve( n, in, out, w, p, h );
    
    ml_free (ker);
}

void FD (double d, int w, double a, double b, int n, double *&in, double *&out )
{
    // differentiation: periodic boundaries, dth derivative, mth order (2*m+1 wide kernel), period a to b, n points.
    
    double scale = pow((b-a)/n,-d);
    
    double *ker = ml_alloc<double > (2*w+1);
    
    for (int j=-w; j<=w; j++ )
	ker[j+w] = FD_Dn_LR(d,w,w,j);
    
    double * p = &(ker[w]);
    periodic_convolve( n, in, out, w, p, scale );
    
    ml_free (ker);
}


void test_function(double **& D, int m, int n, double a, double b, double k )
{
    /*
     * m derivatives, n points, k frequency
     */
    
    if (!D) D = ml_alloc<double > ( m, n );
    
    mp::mp_init(30);
    
    mp_real L = b; L-= a;
    mp_real h = L/n;
    mp_real K = k;
    
    for ( int j=0; j<n; j++ )
    {
        mp_real x = a + h*j;
        
        for ( int d=0; d<m; d++ )
        {
            if (d%2) // sin
                D[d][j] = dble( -pow(K,d)*sin(k*x)* pow(-1.0,floor(((double)(d+1))/2)) );
            else // cos
                D[d][j] = dble( -pow(K,d)*cos(k*x)* pow(-1.0,floor(((double)(d+1))/2)) );
            
            //D[d][j] = real(pow(ml_iu*k,d)*exp(ml_iu*k*dble(x)));
        }
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
    
    return sqrt(sum1/sum2);
}

double rms( double *& x, int n, double a, double b)
{
    double sum = 0.0;	    
    for (int k=0; k<n; k++)
        sum += x[k]*x[k];
    
    return sqrt(sum/n);
}


/*




from math import *
import pylab as p
import sys
sys.path.append("/workspace/mathlib/tools/python/")
from read_file import *

data = read_text( text )

for n in range (30):
    p.plot( [ log10(line[n]) for line in data ] )

p.show()





*/





int main()
{
    std_setup();

    mp::mp_init(100);
    
    
    int n_D = 30;
    int n = 1025;
    double a = 0.0;
    double b = 1.0;
    int m = 24;
    
    double energy[n/2][n_D];
    
    for ( int k_ = 1; k_<n/2; k_++ )
    {
        double k = (_2pi*k_)/(b-a);
        
        //--------
        
        double h = (b-a)/n;
        double **D = ml_alloc<double> ( 1, n );
        double **G = ml_alloc<double> ( n_D, n );
        
        test_function( D, 1, n, a, b, k );
        
        for ( int j=0; j<n; j++ )
            G[0][j] = D[0][j];
        
        for ( int j=1; j<n_D; j++ )
            FD (1, m, a, b, n, G[j-1], G[j] );
        
        for ( int j=1; j<n_D; j++ )
            energy[k_][j] = rms( G[j], n, a, b )*pow(k,-j)*sqrt2;
        
        for ( int j=1; j<n_D; j++ )
            cout << log10(energy[k_][j]) << ", ";
        cout << ";\n";
        
        //sprintf( fname, "/workspace/output/temp/%d", k_ );
        //output( D[0], n, fname );
        
        ml_free(D, 1);
        ml_free(G, n_D);
    }
    
    
    
    
    
    std_exit();
}







    /*
    for ( int k_ = 1; k_<n/2; k_++ )
    {
        for ( int j=1; j<n_D; j++ )
            cout << error[k_][j] << "\t";
        
    }
      */  




/*

        double **E = ml_alloc<double> ( n_D, n );
        
        for ( int j=0; j<n_D; j++ )
            error( D[j], G[j], E[j], n, a, b );
        
        for ( int j=1; j<n_D; j++ )
        {
            sprintf( fname, "/workspace/output/temp/%d", j );
            output( E[j], n, fname );
        }
        
        ml_free(E, n_D);



*/

/*

void error(double *& analytic, double *& numeric, double *& out, int n, double a, double b)
{
    if (!out) out = ml_alloc<double > (n);
    
    double sum = 0.0;

    for (int k=0; k<n; k++)
    {
        out[k] = (fabs( analytic[k]-numeric[k] ));
  //      if ( fabs(analytic[k]-numeric[k])>1E-18 )
  //          out[k] = log10(fabs( analytic[k]-numeric[k] ));
 //       else
 //           out[k] = -2;
    }
}

*/
