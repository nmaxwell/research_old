

#include <mathlib/math/hdaf/apply_hdaf.h>
#include <mathlib/math/differentiation/differentiation.h>



template<class St1,class Xt1,class St2,class Xt2 >
void Del2_fft(
	grid2D<double,St1,Xt1 > & in,
	grid2D<double,St2,Xt2 > & out )
 )
{
    out.copy(in);
    
    grid2D<complex<double> > Q;
    
    double sN=0.0, sum2=0.0;
    
    
    
}




void D_fft (double d, double a, double b, int n, double *&in, double *&out, double eps = 1E-15 )
{
    complex<double > *ft = ml_alloc<complex<double> > ( n/2+1 );
    
    FFT( in, ft, n );
    
    double sN=0.0, sum2=0.0;
    
    for (int k=0; k<n/2+1; k++ )
        sN += abs(ft[k]);
    
    double mark = log(eps)+log(sN);
    
    int k_max = n/2+1;
    double sk = abs(ft[k_max]);
    
    while ( log(sk) < mark and k_max >= 16 )
    {
        k_max--;
        sk += abs(ft[k_max]);
    }
    
    //cout << k_max << endl;
    
    for (int k=0; k<n/2+1; k++ )
        if (k <= k_max)
            ft[k] = ft[k] * complex<double> ( 0.0, ml_2pi*k )/((b-a)*n);
        else ft[k] = 0.0;
    
    IFFT( ft, out, n );
    
    ml_free (ft);
}

















