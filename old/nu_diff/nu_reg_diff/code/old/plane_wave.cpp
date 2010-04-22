
    #define INCLUDE_FD_D1_LR
    #define INCLUDE_FD_D2_LR
    #define INCLUDE_FD_D3_LR
//  #define INCLUDE_FD_D4_LR
//  #define INCLUDE_FD_D5_LR
//  #define INCLUDE_FD_D6_LR

#include <mathlib/math/std_math.h>
#include <mathlib/math/kernels/FDkernels.h>
#include <mathlib/math/hdaf/hdaf.h>
#include <mathlib/math/transforms/FFT.h>
#include <mathlib/math/random/ml_random.h>
#include <mathlib/link.cpp>
#include <mathlib/non-lib_things.h>

#define mp_1 ((mp_real)(1.0))
#define mp_2 ((mp_real)(2.0))
#define mp_pi (atan(mp_1)*4.0)
#define mp_2pi (mp_pi*2.0)
#define mp_4pi2 (mp_2pi*mp_2pi)
#define mp_sqrt2 (sqrt(mp_2))


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


void test_function(double *& d0, double *& d1, double *& d2, double k, int n, double a, double b )
{
    if (!d0) d0 = ml_alloc<double > (n);
    if (!d1) d1 = ml_alloc<double > (n);
    if (!d2) d2 = ml_alloc<double > (n);
    
    double h = (b-a)/n;
    
    for ( int j=0; j<n; j++ )
    {
        double x = a + h*j;
        
        d0[j] = cos(x*k);
        d1[j] = -k*sin(x*k);
        d2[j] = -k*k*cos(x*k);
    }
}

void output(double *& out_array, int n, double a, double b)
{
    static int count = 0;
    
    ofstream out;
    sprintf(fname,"%s/out_%d.py", pwd, count);
    out.open(fname, fstream::out);
    assert(out.good() && out.is_open());
    
    out << "\nimport pylab as p\n";
    out << "X = [";
    for (int k=0; k<n; k++)
        out << a+ ((b-a)*k)/n << ",";
    out << "];\n";
    out << "Y = [";
    for (int k=0; k<n; k++)
        out << out_array[k] << ",";
    out << "];\n";
    out << "p.plot(X,Y) \np.show()\n\n";
    
    count++;
    
    out.close();
}

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

int main()
{
    std_setup();

    mp::mp_init(100);

    if(1)
    {
        int n = 1024;
        double a = 0.0;
        double b = 4.0;
        int m = 6;
        double h = (b-a)/n;
        double k = _2pi*2.0;
        
        double *G0=0,*G1=0,*G2=0,*G3=0,*G4=0,*G5=0,*G6=0,*G7=0;
        
        test_function( G0, G1, G2, k, n, a, b );
        
        FD (1, m, a, b, n, G0, G3 );
        FD (2, m, a, b, n, G0, G4 );
        
        error( G1, G3, G5, n, a, b );
        error( G2, G4, G6, n, a, b );
        
        output( G0, n, a, b );
        output( G5, n, a, b );
        output( G6, n, a, b );
        
        G7 = ml_alloc<double > (n);        
    
        for ( int j=0; j<n; j++ )
        {
            double d_x = a + h*j;
            mp_real mp_x = ((mp_real)a) + (( ((mp_real)b-(mp_real)a)/n )*j);
            
            cout << d_x << "\t" << dble((mp_real)(cos(d_x*k))- cos( mp_x * k )) << endl;
            G7[j] = dble((mp_real)(cos(d_x*k))- cos( mp_x * k ));
        }
        
        output( G7, n, a, b );
        
        ml_free(G1);
        ml_free(G2);
        ml_free(G3);
        ml_free(G4);
        ml_free(G5);
        ml_free(G6);
        ml_free(G7);
    }
    
    std_exit();
}



