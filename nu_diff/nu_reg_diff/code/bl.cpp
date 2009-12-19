


#define N_FFT_THREADS 1
#define INCLUDE_FD_D1_LR

#include <mathlib/math/random/ml_random.h>

#include "common.cpp"






void test_function(double **& D, int n_d, int n, double a, double b, double k_max )
{    
    ml_random rng;
    
    complex<double > ** Q = ml_alloc<complex<double> > ( n_d, n/2+2);
    
    for (int k=0; k<n/2+1; k++)
	if ( k*_2pi/(b-a) < k_max && k!=0 )
	    Q[0][k] = exp(iu*_2pi*rng.gen_double())*rng.gen_double();
	else
	    Q[0][k] = 0;
    
    for (int j=1; j<n_d; j++)
    for (int k=0; k<n/2+1; k++)
        Q[j][k] = Q[0][k]*pow( iu*((double)k)*_2pi/(b-a), j ) ;
	
    for (int j=0; j<n_d; j++)
        IFFT( Q[j], D[j], n );
        
    for (int j=0; j<n_d; j++)
    for (int k=0; k<n; k++)
        D[j][k] /= n;
    
    ml_free(Q, n_d);
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
    
    cout << k_max << endl;
    
    for (int k=0; k<n/2+1; k++ )
        if (k <= k_max)
            ft[k] = ft[k] * complex<double> ( 0.0, ml_2pi*k )/((b-a)*n);
        else ft[k] = 0.0;
    
    IFFT( ft, out, n );
    
    ml_free (ft);
}


int main()
{
    mp::mp_init(30);
    
    int n = 512;
    double a = 0.0;
    double b = 1.0;
    double h = (b-a)/n;
    
    int n_d = 100;
    
  //  double **error = ml_alloc<double> ( n/2, n_d );
    
    {
        double k_max = (pi*n/(b-a))*(0.25);
        
        //--------
        
        double **ana = ml_alloc<double> ( n_d, n );
        double **num = ml_alloc<double> ( n_d, n );
        double *mag = ml_alloc<double> ( n_d );
        double **err = ml_alloc<double> ( n_d, n );
        double **temp = ml_alloc<double> ( n_d, n );
        complex<double> **ft_temp = ml_alloc<complex<double> > ( n_d, n/2+1 );
        
        test_function( ana, n_d, n, a, b, k_max );
        
        for ( int j=0; j<n; j++ )
            num[0][j] = ana[0][j];
        
        for ( int j=1; j<n_d; j++ )
            D_fft (1, a, b, n, num[j-1], num[j] );
        
        for ( int j=0; j<n_d; j++ )
            mag[j] = l2_norm( ana[j], n, a, b ) + 1E-18;
        
        for ( int j=1; j<n_d; j++ )
        for ( int i=0; i<n; i++ )
            err[j][i] = log10( fabs( num[j][i] - ana[j][i] ) + 1E-18) -log10(mag[j]);
        
        output( ana, n_d, n, "/workspace/output/temp/ana" );
        output( num, n_d, n, "/workspace/output/temp/num" );
        output( err, n_d, n, "/workspace/output/temp/err" );
        
        
        for ( int j=0; j<n_d; j++ )
        for ( int i=0; i<n; i++ )
            err[j][i] = fabs( num[j][i] - ana[j][i] )/mag[j];
        
        for ( int j=0; j<n_d; j++ )
            FFT( err[j], ft_temp[j], n );
        
        for ( int j=0; j<n_d; j++ )
        for ( int i=0; i<n/2+1; i++ )
            temp[j][i] = log10(abs(ft_temp[j][i]) +1E-18 );
        
        output( temp, n_d, n/2, "/workspace/output/temp/ft_err" );
        
        
        for ( int j=0; j<n_d; j++ )
            FFT( ana[j], ft_temp[j], n );
        
        for ( int j=0; j<n_d; j++ )
        for ( int i=0; i<n/2+1; i++ )
            temp[j][i] = log10(abs(ft_temp[j][i]) +1E-18 );
        
        output( temp, n_d, n/2, "/workspace/output/temp/ft_ana" );
        
        
        for ( int j=0; j<n_d; j++ )
            FFT( num[j], ft_temp[j], n );
        
        for ( int j=0; j<n_d; j++ )
        for ( int i=0; i<n/2+1; i++ )
            temp[j][i] = log10(abs(ft_temp[j][i]) +1E-18 );
        
        output( temp, n_d, n/2, "/workspace/output/temp/ft_num" );
        
        
        
        ml_free( mag );
        ml_free(ana, n_d);
        ml_free(num, n_d);
        ml_free(err, n_d);
        ml_free(temp, n_d);
        ml_free(ft_temp, n_d);
    }
    
    
    
    //sprintf( fname, "/workspace/output/temp/%d", g );
  //  output ( error, k_max, n_d-1,  "/workspace/output/temp/out" );
    
  //  ml_free( error, k_max );    
    
}


