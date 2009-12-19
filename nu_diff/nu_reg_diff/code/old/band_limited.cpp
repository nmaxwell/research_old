
    #define INCLUDE_FD_D1_LR
    #define INCLUDE_FD_D2_LR
//  #define INCLUDE_FD_D3_LR
//  #define INCLUDE_FD_D4_LR
//  #define INCLUDE_FD_D5_LR
//  #define INCLUDE_FD_D6_LR

#include <mathlib/math/std_math.h>
#include <mathlib/math/kernels/FDkernels.h>
#include <mathlib/math/hdaf/hdaf.h>
#include <mathlib/math/transforms/FFT.h>
#include <mathlib/math/random/ml_random.h>
#include <mathlib/link.cpp>

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


void test_function(double *& d0, double *& d1, double *& d2, int n, double a, double b, double k_max )
{    
    ml_random rng;
    
    complex<double > * Q1 = ml_alloc<complex<double> > (n/2+2);
    complex<double > * Q2 = ml_alloc<complex<double> > (n/2+2);
    
    for (int k=0; k<=n/2+1; k++)
	if ( k*_2pi/(b-a) < k_max )
	    Q1[k] = exp(iu*_2pi*rng.gen_double())*rng.gen_double();
	else
	    Q1[k] = 0;
    
    for (int k=0; k<=n/2+1; k++)
	Q2[k] = Q1[k];
    
    IFFT( Q2, d0, n );
    
    for (int k=0; k<=n/2+1; k++)
	Q1[k] *= iu*((double)k)*_2pi/(b-a);
    
    for (int k=0; k<=n/2+1; k++)
	Q2[k] = Q1[k];
	
    IFFT( Q2, d1, n );
    
    for (int k=0; k<=n/2+1; k++)
	Q1[k] *= iu*((double)k)*_2pi/(b-a);
    
    for (int k=0; k<=n/2+1; k++)
	Q2[k] = Q1[k];
	
    IFFT( Q2, d2, n );
    
    ml_free(Q1);
    ml_free(Q2);
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
	//cout << k << "\t" <<  ((analytic[k]-numeric[k])/fabs(analytic[k])) << "\t" <<  log10(fabs((analytic[k]-numeric[k])/analytic[k])) <<  endl;
        if (fabs(analytic[k]-numeric[k])>1E-18 && fabs(analytic[k])>1E-18 )
	    out[k] = log10( fabs((analytic[k]-numeric[k])/analytic[k]) );
	else
	    //out[k] = (analytic[k]-numeric[k]);
	    out[k] = -17;
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



    if (1)
    {
	int n = 2048;
	double a = -1.0;
	double b = 1.0;
	double gamma = 0.25;
	int w = 20;
	int n_trials = 10;
	
	
	for (gamma = 0.05; gamma <=1.01; gamma += 0.05)
	{
	    double *G0=0,*G1=0,*G2=0,*G3=0,*G4=0,*G5=0,*G6=0;
	    
	    test_function( G0, G1, G2, n, a, b, gamma*(pi*n)/(b-a) );
        
	    hdaf_diff( a, b, n, 30, 0.35, G0, G3 );
	    
	    double error = 0.0;
	    for (int t=0; t<n_trials; t++)
		error += log10(l2_error(G1,G3, n,a,b ));
	    error /= n_trials;
	    
	    cout << log2(n) << "\t" << gamma << "\t" << error << endl;
	    
	    ml_free(G1);
	    ml_free(G2);
	    ml_free(G3);
	    ml_free(G4);
	    ml_free(G5);
	    ml_free(G6);
	}
    }



    if(0)
    {
	int n = 256*8;
	double a = -1.0;
	double b = 1.0;
	int m = 12;
	
	double *G0=0,*G1=0,*G2=0,*G3=0,*G4=0,*G5=0,*G6=0;
	
	test_function( G0, G1, G2, n, a, b, 0.01*(pi*n)/(b-a) );
	
	hdaf_diff( a, b, n, m, 0.25, G0, G3 );
	//FD (1, w, a, b, n, G0, G3 );
	
	error( G1, G3, G6, n, a, b );
	
	output( G0, n, a, b );
	output( G1, n, a, b );
	output( G3, n, a, b );
	output( G6, n, a, b );
	
	ml_free(G1);
	ml_free(G2);
	ml_free(G3);
	ml_free(G4);
	ml_free(G5);
	ml_free(G6);
    }


    if(0)
    {
	int n = 512;
	double a = -1.0;
	double b = 1.0;
	double gamma = 0.25;
	int w = 20;
	int n_trials = 4;
	
	for (int w=1; w<=20; w++)
	for (int k=7; k<=12; k++)
	{
	    double error = 0.0;
	    int trials = 0;
	    
	    for (gamma = 0.02; gamma <=0.25; gamma += 0.02)
	    {
		n = pow(2,k);
		double *G0=0,*G1=0,*G2=0,*G3=0,*G4=0,*G5=0,*G6=0;
		
		test_function( G0, G1, G2, n, a, b, gamma*(pi*n)/(b-a) );
		
		FD (1, w, a, b, n, G0, G3 );	    
		
		double error2 = 0.0;
		for (int t=0; t<n_trials; t++)
		    error2 += log10(l2_error(G1,G3, n,a,b ));
		error2 /= n_trials;
		
		error += error2;
		trials++;		
		//cout << log2(n) << "\t" << gamma << "\t" << error << endl;
		
		ml_free(G1);
		ml_free(G2);
		ml_free(G3);
		ml_free(G4);
		ml_free(G5);
		ml_free(G6);
	    }
	    
	    error /= trials;
	    
	    cout << w << "\t" << k << "\t" << error << endl;
	}
    }


    if (1)
    {
	int n = pow(2,11);
	double a = -1.0;
	double b = 1.0;
	double gamma = 0.25;
	int w = 20;
	int n_trials = 10;
	
	
	for (gamma = 0.05; gamma <=0.75; gamma += 0.05)
	{
	    double *G0=0,*G1=0,*G2=0,*G3=0,*G4=0,*G5=0,*G6=0;
	    
	    test_function( G0, G1, G2, n, a, b, gamma*(pi*n)/(b-a) );
	    
	    FD (1, w, a, b, n, G0, G3 );	    
	    
	    double error = 0.0;
	    for (int t=0; t<n_trials; t++)
		error += log10(l2_error(G1,G3, n,a,b ));
	    error /= n_trials;
	    
	    cout << log2(n) << "\t" << gamma << "\t" << error << endl;
	    
	    ml_free(G1);
	    ml_free(G2);
	    ml_free(G3);
	    ml_free(G4);
	    ml_free(G5);
	    ml_free(G6);
	}
    }
    
	
    std_exit();
}



