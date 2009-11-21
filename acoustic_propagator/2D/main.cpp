
#define N_FFT_THREADS 1


#include <mathlib/math/grids/grid2D.h>
#include <mathlib/math/grids/grid_file.cpp>
#include <mathlib/math/grids/extra/plots.cpp>

#include "Del2_hdaf.cpp"

#include <mathlib/math/grids/grid2D_Del2_FD.h>


#include <mathlib/link.cpp>

color3 cmap(double x)
{
	x *= 10;
	float s = atan(x)/pi+0.5;
	return color3(s,s,s);
}

void output(double t, grid2D< > & u, grid2D< > & v)
{
	static int n = 0;
	//sprintf(fname,"%s/out_dat/u_%05d.dat",pwd,n);
    sprintf(fname,"%s/out_png/u_%05d.png",pwd,n);
	plotGrid2D_1(u,fname,cmap);
	//writeFile(u,fname);
	n++;
}


void P(double t, grid2D<> & C2, int M, double k_cutoff, grid2D<> & u0, grid2D<> & v0, grid2D<> & u1, grid2D<> & v1)
{
    int diff_order = 16;
    //cout << "\t" << L2norm(u0) << "\t" << L2norm(v0) << endl;
    u1.copy(u0);
    v1.copy(v0);
    
    static int M_ = 0;
    static grid2D<> * U = 0;
    static grid2D<> * V = 0;
    
    if (M_ != M)
    {
        M_ = M;
        
        if (U) delete [] U;
        if (V) delete [] V;
        
        U = new grid2D<> [M+2];
        V = new grid2D<> [M+1];
    }
    
    int m_x = 16;
    int m_y = 16;
    
    
    U[0] = u0;
    V[0] = v0;
    
    double t1 = get_real_time();
    
    for (int n = 1; n<=M+1; n++)
    {
        //Del2_FD(U[n-1], U[n], diff_order, GRID_NONPER_BD_1, 24 );        
        Del2_hdaf( U[n-1], U[n], k_cutoff, m_x, m_y );
        
        U[n] *= C2;
        U[n] *= t*t/(4.0*n*n-2.0*n);
    }
    
    for (int n = 1; n<=M; n++)
    {
        //Del2_FD(V[n-1], V[n], diff_order, GRID_NONPER_BD_1, 24 );        
        Del2_hdaf( V[n-1], V[n], k_cutoff, m_x, m_y );
        
        V[n] *= C2;
        V[n] *= t*t/(4.0*n*n-2.0*n);
    }
    
    double t2 = get_real_time();
    
    //cout << "\t" << t2-t1 << endl;	
    
    u1 = 0.0;
    v1 = 0.0;
    for (int n = 0; n<=M; n++)
    GRID2D_LINLOOP( u0 )
    {
        u1[i] += U[n][i];
        u1[i] += V[n][i]* t/(2.0*n+1.0);
        
        v1[i] += U[n+1][i]* (4.0*(n*n)+6.0*n+2.0)/(t*(2.0*n+1.0));
        v1[i] += V[n][i];
    }
}






class u_0_func : public functor2<double,double >
{	
public:
	double operator() (double const & x, double const & y) const 
	{
	//	return 1.0*exp(-(x*x+y*y));
        
		double s = sqrt(x*x+y*y);
		if (s <= 0.5)
			return cos(_2pi*s)+1.0;
		else return 0.0;
	}
};


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
    	double c = 3.0;
		return (1.0- (step((x-7.0)/2.0)*step(y/5.0))*0.7 )*c*c;
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



int main()
{
    std_setup();
    
    int N = 512;
    
    double xa = -40.0;
    double xb = +40.0;
    
    double ya = -40.0;
    double yb = +40.0;
    
    double ixa = -25.0;
    double ixb = +25.0;
    
    double iya = -25.0;
    double iyb = +25.0;
    
    double tf = 10.0;
    double dt = 0.1;
    
    grid2D<double,double,double > grid(N,xa,xb, N,xa,xb),u,v,C2,damp;
    grid = 0.0;
    
    u = grid;
    v = grid;
    C2 = grid;
    
    C2 = C2_func();
    u = u_0_func();
    
    damp = grid;
    damp = damping_function(1.0, xa,xb,ya,yb, ixa,ixb,iya,iyb, 4, 2  );
        
    sprintf(fname,"%s/out_dat/C2.dat",pwd);	
    writeFile(C2,fname);
    
    sprintf(fname,"%s/damp.png",pwd);
    plotGrid2D_1(damp,fname,cmap);
    
    double k_cutoff = 20.0;
    
    double t = 0.0;
    
    while (t <= tf)
    {
        output(t,u,v);
        
            double t1 = get_real_time();
        
        P(dt, C2, 12, k_cutoff, u,v,u,v );
        
        u *= damp;
        v *= damp;
        
            double t2 = get_real_time();        
            cout << t << "\t" << t2-t1 << endl;
        
        t += dt;
    }
    
    output(t,u,v);
}

