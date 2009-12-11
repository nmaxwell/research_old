


//#define FFTW_PLAN_MODE FFTW_PATIENT

#define N_FFT_THREADS  3




#include <mathlib/math/std_math.h>
#include <mathlib/math/grids/grid2D.h>
#include <mathlib/math/grids/grid_file.cpp>
#include <mathlib/math/transforms/fft.h>
#include <mathlib/math/laplacian/laplacian.h>
#include <mathlib/math/grids/extra/plots.cpp>

#include <mathlib/non-lib_things.h>

#include <mathlib/link.cpp>


void output(double t, grid2D< > & u, grid2D< > & v)
{
	static int n = 0;
    sprintf(fname,"/workspace/output/acoustic_propagate_2d/out_dat/u_%05d.dat", n);
    
	assert(!writeFile(u,fname));
	n++;
}

void propagate(double t, grid2D<> & C2, int exp_order, grid2D<> & u0, grid2D<> & v0, grid2D<> & u1, grid2D<> & v1)
{
    u1.copy(u0);
    v1.copy(v0);
    
    int n1 = u0.n1;
    int n2 = u0.n2;
    double L1 = u0.b1-u0.a1;
    double L2 = u0.b2-u0.a2;
    
    static int exp_order_ = 0;
    static grid2D<> * U = 0;
    static grid2D<> * V = 0;
    
    if (exp_order_ != exp_order)
    {
        exp_order_ = exp_order;
        
        if (U) delete [] U;
        if (V) delete [] V;
        
        U = new grid2D<> [exp_order+2];
        V = new grid2D<> [exp_order+1];
    }
    
    U[0] = u0;
    V[0] = v0;
    
    for (int n = 1; n<=exp_order+1; n++)
    {
        laplacian_2d_fft ( U[n-1].array, U[n].array, n1, n2, L1, L2 );
        
        U[n] *= C2;
        U[n] *= t*t/((4.0*n-2.0)*n);
    }
    
    for (int n = 1; n<=exp_order; n++)
    {
        laplacian_2d_fft ( V[n-1].array, V[n].array, n1, n2, L1, L2 );
        
        V[n] *= C2;
        V[n] *= t*t/((4.0*n-2.0)*n);
    }
    
    u1 = 0.0;
    v1 = 0.0;
    
    for (int n = 0; n<=exp_order; n++)
    GRID2D_LINLOOP( u0 )
    {
        u1[i] += U[n][i];
        u1[i] += V[n][i]* t/(2.0*n+1.0);
        
        v1[i] += U[n+1][i]*((4.0*n+6.0)*n+2.0)/(t*(2.0*n+1.0));
        v1[i] += V[n][i];
    }
}


class u_0_func : public functor2<double,double >
{	
public:
	double operator() (double const & x1, double const & x2) const 
	{
		//return 1.0*exp(-(x*x+y*y));
        double k1 = 1.0;
        double k2 = 0.0;
        
        return cos( k1*x2+k2*x2 );
        
	/*	double s = sqrt(x*x+y*y);
		if (s <= 0.5)
			return cos(_2pi*s)+1.0;
		else return 0.0;*/
	}
};




class u_t_func : public functor2<double,double >
{	
public:
    double t;
    u_t_func( double t):t(t) {}
public:
	double operator() (double const & x1, double const & x2) const 
	{
		//return 1.0*exp(-(x*x+y*y));
        double k1 = 1.0*_2pi;
        double k2 = 0.0*_2pi;
        double c = 1.0;
        
        return cos( k1*x2+k2*x2 )*cos(t*c*sqrt(k1*k1+k2*k2));
        
	/*	double s = sqrt(x*x+y*y);
		if (s <= 0.5)
			return cos(_2pi*s)+1.0;
		else return 0.0;*/
	}
};


int main()
{
    std_setup();
    
    int n = 256;
    double x0 = 1.0;
    
    double tf = 1.0;
    double dt = 0.01;
    
    grid2D<double,double,double > grid( n,x0,-x0, n,x0,-x0 ),u,v,C2,damp;
    grid = 0.0;
    
    u = grid;
    v = grid;
    C2 = grid;
    
    C2 = 1.0;
    u = u_0_func();
    
    sprintf(fname,"/workspace/output/acoustic_propagate_2d/out_dat/C2.dat", n);
    writeFile(C2,fname);
    
    double t = 0.0;
    
    while (t <= tf)
    {
        u = u_t_func(t);
        v = 0;
        output(t,u,v);
        
        double t1 = get_real_time();
        
        propagate (dt, C2, 1, u,v,u,v );
        
        double t2 = get_real_time();        
        cout << t << "\t" << t2-t1 << endl;
        
        t += dt;
    }
    u = u_t_func(t);
        v = 0;
    
    output(t,u,v);
}
















/*







int main()
{
    std_setup();
    
    int n = 128;
    double x0 = 10.0;
    
    double xi = 25.0;
    
    int n1 = n;
    int n2 = n;
    
    double a1 = -x0;
    double b1 = +x0;
    
    double a2 = -x0;
    double b2 = +x0;
    
    double ia1 = -xi;
    double ib1 = +xi;
    
    double ia2 = -xi;
    double ib2 = +xi;
    
    double tf = 10.0;
    double dt = 0.1;
    
    grid2D<double,double,double > grid( n1,a1,b1, n2,a2,b2 ),u,v,C2,damp;
    grid = 0.0;
    
    u = grid;
    v = grid;
    C2 = grid;
    
    C2 = C2_func();
    u = u_0_func();
    
    damp = grid;
    damp = damping_function(1.0, a1,b1,a2,b2, ia1,ib1,ia2,ib2, 4, 2  );
    
    sprintf(fname,"/workspace/output/acoustic_propagate_2d/out_dat/C2.dat", n);
    writeFile(C2,fname);
    
    sprintf(fname,"/workspace/output/acoustic_propagate_2d/out_dat/damp.dat", n);
    writeFile(damp,fname);
    
    double t = 0.0;
    
    while (t <= tf)
    {
        output(t,u,v);
        
        double t1 = get_real_time();
        
        propagate (dt, C2, 3, u,v,u,v );
        
        //u *= damp;
        //v *= damp;
        
        double t2 = get_real_time();        
        cout << t << "\t" << t2-t1 << endl;
        
        t += dt;
    }
    
    output(t,u,v);
}








double step(double x)
{
	//return (pow(x,6)/6.0+pow(x,4)/2.0+pow(x,2)+1.0)*exp(-x*x);
	if (fabs(x) < 1.0)
		return 1.0;
	else
		return 0.0; 
}

class C2_func : public virtual functor2<double,double >
{	
public:
	double operator() (double const & x,double const & y) const
	{
    	double c = 1.0;
		//return (1.0- (step((x-7.0)/2.0)*step(y/5.0))*0.7 )*c*c;
        return c;
	}
};

double damping_profile(double t, double x, int n, double a, double b, int m )
{
    // returns the profile for a damping function of the form
    // f(x) = exp(-t*x^n)
    // this is transformed such that
    // f(a) = 1, f(b) = eps, f(x<=a) = 1, 0 < f(x>a) < 1 , at t=1
    // where eps = 10^(-m)
    
    if ( a < b && x <= a ) return 1.0;
    if ( b < a && x >= a ) return 1.0;
    return exp( -t*pow((x-a)/(b-a),n)*log(10.0)*m );
}

class damping_function : public functor2<double,double >
{
public:
	double t;
	int n,m;
	double xa,xb,ya,yb; // dimensions of grid
    double ixa,ixb,iya,iyb; // dims of interiod where no damping
    
	damping_function(
	    double t,
        double xa, double xb, double ya, double yb,
        double ixa, double ixb, double iya, double iyb,
        int n, int m )
        :t(t),xa(xa),xb(xb),ya(ya),yb(yb),ixa(ixa),ixb(ixb),iya(iya),iyb(iyb),n(n),m(m) {};
	
public:
	double operator() (double const & x,double const & y) const
	{
    	double R = 1.0;
    	R *= damping_profile(t, x, n, ixa, xa, m );
    	R *= damping_profile(t, y, n, iya, ya, m );
    	R *= damping_profile(t, x, n, ixb, xb, m );
    	R *= damping_profile(t, y, n, iyb, yb, m );
		return R;
	}
};




*/


