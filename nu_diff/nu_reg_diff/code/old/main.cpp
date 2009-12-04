
    #define INCLUDE_FD_D1_LR
    #define INCLUDE_FD_D2_LR
//  #define INCLUDE_FD_D3_LR
//  #define INCLUDE_FD_D4_LR
//  #define INCLUDE_FD_D5_LR
//  #define INCLUDE_FD_D6_LR

#include <mathlib/math/std_math.h>
#include <mathlib/math/kernels/FDkernels.h>
#include <mathlib/math/hdaf/hdaf.h>
#include <mathlib/link.cpp>

#define mp_1 ((mp_real)(1.0))
#define mp_pi (atan(mp_1)*4.0)
#define mp_2pi (mp_pi*2.0)
#define mp_4pi2 (mp_2pi*mp_2pi)


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
        
        for ( int j=-M; j<=M; j++ )
            sum += in[periodic_mod(k-j,n)]*kernel[j];
        
        out[k] = sum*scale;
    }
}

void down_sample( int n_in, double *& in, int rate, int n_out, double *& out,
		    int l, int r, double * kernel, double scale )
{
    if (!out) out = ml_alloc<double > ( n_out );
    
    for (int k=0; k<n_out; k++)
    {
        double sum = 0.0;
        
        for (int j=-r; j<=l; j++)
            sum += in[periodic_mod(rate*k-j,n_in)]*kernel[j];
        
        out[k] = sum*scale;
    }    
}

void hdaf_down_sample(
    int m, double a, double b, int n_in, double *& in, int rate, int n_out, double *& out )
{
    if (!out) out = ml_alloc<double > ( n_out );
    
    double h = (b-a)/n_in;
    double gamma = 1.0/rate;    
    double sigma = h*sqrt(2*m+1)/(pi*gamma);
    hdaf_delta<mp_real > delta(m,sigma);
    
    int K = (int)ceil(hdaf_truncate_point (1E-14, 5E-14, 8,sigma, 0 )/h);
    double *ker = ml_alloc<double > (2*K+1);
    
    for (int j=-K; j<=K; j++ )
	ker[j+K] = dble( delta(h*j) );
    
    double * p = &(ker[K]);
    down_sample( n_in, in, rate, n_out, out,  K,K , p, h );
    
    ml_free (ker);
}


void hdaf_lp( double a, double b, int n, int m, double gamma, double *& in, double *& out )
{
    double h = (b-a)/n;
    
    double sigma = h*sqrt(2*m+1)/(pi*gamma);
    hdaf_delta<mp_real > delta(m,sigma);
    
    int K = (int)ceil(hdaf_truncate_point (1E-8, 5E-8, 8,sigma, 0 )/h);
    
    double *ker = ml_alloc<double > (2*K+1);
    
    for (int j=-K; j<=K; j++ )
	ker[j+K] = dble( delta(h*j) );
    
    double * p = &(ker[K]);
    periodic_convolve( n, in, out, K, p, h );
	
    ml_free (ker);
}

void differentiate (double d, double m, double a, double b, int n, double *&in, double *&out )
{
    // differentiation: periodic boundaries, dth derivative, mth order (2*m+1 wide kernel), period a to b, n points.
    
    if (!out) out = ml_alloc<double > (n);
    
    double scale = pow((b-a)/n,-d);
    
    
    for (int k=0; k<n; k++)
    {
        mp_real sum = 0.0;
        
        for (int j=-m; j<=m; j++)
            sum += (mp_real)(in[periodic_mod(k-j,n)])*(mp_real)(FD_Dn_LR(d,m,m,j));
        
        out[k] = dble(sum*scale);
    }
}




void filtered_differentiate (
	double d, double w, double gamma, int m,
	double a, double b, int n, double *&in, double *&out )
{
    if (!out) out = ml_alloc<double > (n);
    double *temp = ml_alloc<double> (n);
    
    differentiate (d, w, a, b, n, in, temp );
    
    if(1)
	hdaf_lp( a, b, n, m, gamma, temp, out );
    
    ml_free(temp);
}


void test_function(double *& out, int n, double a, double b)
{
    if (!out) out = ml_alloc<double > (n);
    
    for (int k=0; k<n; k++)
    {
        double x = a+ ((b-a)*k)/n;
        
        out[k] = exp(cos(x*_2pi)-1);
    }
}

void test_function_d1(double *& out, int n, double a, double b)
{
    if (!out) out = ml_alloc<double > (n);
    
    for (int k=0; k<n; k++)
    {
        double x = a+ ((b-a)*k)/n;
        
        out[k] = dble(-mp_2pi*sin(mp_2pi*x)*exp(cos(x*mp_2pi)-mp_1));
    }
}

void test_function_d2(double *& out, int n, double a, double b)
{
    if (!out) out = ml_alloc<double > (n);
    
    for (int k=0; k<n; k++)
    {
        double x = a+ ((b-a)*k)/n;
        
        out[k] = dble(exp(-mp_1)*mp_4pi2*(sin(x*mp_2pi)*sin(x*mp_2pi)-cos(x*mp_2pi))*exp(cos(x*mp_2pi)));
    }
}

void output(double *& out, int n, double a, double b)
{
    cout << "\nimport pylab as p\n";
    cout << "X = [";
    for (int k=0; k<n; k++)
        cout << a+ ((b-a)*k)/n << ",";
    cout << "];\n";
    cout << "Y = [";
    for (int k=0; k<n; k++)
        cout << out[k] << ",";
    cout << "];\n";
    cout << "p.plot(X,Y) \np.show()\n\n";
}

void error(double *& analytic, double *& numeric, double *& out, int n, double a, double b)
{
    if (!out) out = ml_alloc<double > (n);
    
    double sum = 0.0;
    
    for (int k=0; k<n; k++)
		sum += fabs(analytic[k]);//*analytic[k];
	
	//sum = sqrt(sum);
	    
    for (int k=0; k<n; k++)
    {
        double x = a+((b-a)*k)/n;
        
        out[k] = (analytic[k]-numeric[k])/sum	;
    }
    
    
 /*   double sum = 0.0;
	    
    for (int k=0; k<n; k++)
    {
		//cout << k << "\t" << analytic[k]-numeric[k] << "\t" << analytic[k] << endl;
        if (fabs(analytic[k]-numeric[k])>1E-18 && fabs(analytic[k])>1E-18 )
			out[k] = (analytic[k]-numeric[k])/fabs(analytic[k]);
		else
			//out[k] = (analytic[k]-numeric[k]);
			out[k] = 1E-17;
    }*/
}


double error2(double *& analytic, double *& numeric,  int n, double a, double b)
{    
    double sum1 = 0.0;
    
    for (int k=0; k<n; k++)
	sum1 += (analytic[k]-numeric[k])*(analytic[k]-numeric[k]);
    
    double sum2 = 0.0;    
    
    for (int k=0; k<n; k++)
	sum2 += analytic[k]*analytic[k];
    
    return sqrt(sum1/sum2);
}

int main()
{
    std_setup();
    
    mp::mp_init(100);
    
    
    int n = 128;
    double a = -1.0;
    double b = 1.0;
    int os_rate = 4;
    for ( int k=2; k<=14; k++ )
    {
	n = pow(2,k);
	
	double *G1=0,*G2=0,*G3=0,*G4=0,*G5=0,*G6=0;
	
	test_function( G1, n, a, b );
	test_function( G4, n*4, a, b );
	
	differentiate (1, 8, a, b, n, G1, G2 );
	differentiate (1, 8, a, b, n*os_rate, G4, G5 );
	
	hdaf_down_sample( 20, a, b, n*os_rate, G5, os_rate, n, G6 );
	
        /*filtered_differentiate (
	    1, 20, 0.25, 6,
	    a, b, n, G1, G2 );*/
	
	test_function_d1( G3, n, a, b );
	
	cout << k << "\t" << error2( G3, G2,  n, a, b) << endl;
	cout << k << "\t" << error2( G3, G6,  n, a, b) << endl;
	cout << endl;
	
	//error( G3, G2, G4, n, a, b );
	
	//output( G4, n, a, b  );  
	
	ml_free(G1);
	ml_free(G2);
	ml_free(G3);
	ml_free(G4);
	ml_free(G5);
	ml_free(G6);
    }
    std_exit();
}















/*   if (!out) out = ml_alloc<double > (n);
    
    int n1 = n;
    int n2 = n*oversample_rate;
    
    double h1 = (b-a)/n1;
    double h2 = (b-a)/n2;
    
    double *Y=0,*X = ml_alloc<double> (n2);
    
    for (int k=0; k<n2; k++ )
	{
		double x = ((double)k)/oversample_rate;
		
		double sum = 0.0;
		
		for (int j=0; j<n1; j++ )
			sum += in[j]*norm_sinc(x-j);
		
		X[k] = sum;
	}
	
	differentiate (d, m, a, b, n2, X, Y );
	
	double gamma = 1.0/oversample_rate;
	int m_hdaf = 10;
	
	double sigma = h2*sqrt(2*m_hdaf+1)/(pi*gamma);
	hdaf_delta<mp_real > delta(m,sigma);
	
	int K = (int)ceil(hdaf_truncate_point (1E-8, 5E-8, 8,sigma, 0 )/h2);
//	cout << K << endl;
	
    double *ker = ml_alloc<double > (2*K+1);
    
    for (int j=-K; j<=K; j++ )
		ker[j+K] = dble( delta(h2*j) );
    
    double * p = &(ker[K]);
    periodic_convolve( n2, Y, X, K, p, h2 );
    
	ml_free (ker);
	
	for (int k=0; k<n1; k++ )
		out[k] = X[k*oversample_rate];
	
	ml_free (X);
	ml_free (Y);
	*/
