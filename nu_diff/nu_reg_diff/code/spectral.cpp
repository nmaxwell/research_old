
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



void D_fft (double d, int w, double a, double b, int n, double *&in, double *&out )
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
    
    int k_max = 64;
    
    int n_d = 20;
    
    double **error = ml_alloc<double> ( n/2, n_d );
    
    for ( int k_ = 0; k_<k_max; k_++ )
    {
        double k = (_2pi*k_)/(b-a);
        
        //--------
        
        double **D = ml_alloc<double> ( n_d, n );
        double **G = ml_alloc<double> ( n_d, n );
        
        test_function( D, n_d, n, a, b, k );
        
        for ( int j=0; j<n; j++ )
            G[0][j] = D[0][j];
        
        for ( int j=1; j<n_d; j++ )
            D_fft (1, 16, a, b, n, G[j-1], G[j] );
            
        for ( int j=1; j<n_d; j++ )
            error[k_][j-1] = l2_error( D[j], G[j], n, a, b );
        
        sprintf( fname, "/workspace/output/temp/%d", k_ );
        output( D, n_d, n, fname );
        
        ml_free(D, n_d);
        ml_free(G, n_d);
    }
    
    //sprintf( fname, "/workspace/output/temp/%d", g );
    output ( error, k_max, n_d-1,  "/workspace/output/temp/out" );
    
    
    
    ml_free( error, k_max );    
    
    
    
    
    
    
    
    
}


