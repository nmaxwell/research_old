


#define N_FFT_THREADS 1
#define INCLUDE_FD_D1_LR

#include "common.cpp"



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



void D_fft (double d, double a, double b, int n, double *&in, double *&out )
{
    // differentiation: periodic boundaries, dth derivative, mth order (2*m+1 wide kernel), period a to b, n points.
    
    double scale = pow((b-a)/n,-d);
    
    complex<double > *ft = ml_alloc<complex<double> > ( n/2+1 );
    
    FFT( in, ft, n );
    
    for (int k=0; k<n/2+1; k++ )
        ft[k] = ft[k] * complex<double> ( 0.0, ml_2pi*k )/((b-a)*n);
    
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
    
    int k_max = n/2;
    
    int n_d = 20;
    
    double **error = ml_alloc<double> ( n/2, n_d );
    
    for ( int k_ = 0; k_<k_max; k_++ )
    {
        cout << k_ << endl;
        double k = (_2pi*k_)/(b-a);
        
        //--------
        
        double **ana = ml_alloc<double> ( n_d, n );
        double **num = ml_alloc<double> ( n_d, n );
        double *mag = ml_alloc<double> ( n_d );
        double **err = ml_alloc<double> ( n_d, n );
        complex<double> **ft_err = ml_alloc<complex<double> > ( n_d, n/2+1 );
        
        test_function( ana, n_d, n, a, b, k );
        
        for ( int j=0; j<n; j++ )
            num[0][j] = ana[0][j];
        
        for ( int j=1; j<n_d; j++ )
            D_fft (1, a, b, n, num[j-1], num[j] );
        
        
        
        for ( int j=1; j<n_d; j++ )
            error[k_][j-1] = l2_error( ana[j], num[j], n, a, b );
        
        for ( int j=0; j<n_d; j++ )
            mag[j] = l2_norm( ana[j], n, a, b ) + 1E-18;
        
        
        for ( int j=1; j<n_d; j++ )
        for ( int i=0; i<n; i++ )
            err[j][i] = fabs( num[j][i] - ana[j][i] )/mag[j] ;
        
        for ( int j=1; j<n_d; j++ )
            FFT( err[j], ft_err[j], n );
        
        for ( int j=1; j<n_d; j++ )
        for ( int i=0; i<n/2+1; i++ )
            err[j][i] = log(abs(ft_err[j][i]) +1E-18 );
        
        sprintf( fname, "/workspace/output/temp/%d", k_ );
        output( err, n_d, n/2+1, fname );
        
        
        
        
        
        ml_free( mag );
        ml_free(ana, n_d);
        ml_free(num, n_d);
        ml_free(err, n_d);
    }
    
    
    
    //sprintf( fname, "/workspace/output/temp/%d", g );
    output ( error, k_max, n_d-1,  "/workspace/output/temp/out" );
    
    ml_free( error, k_max );    
    
}


