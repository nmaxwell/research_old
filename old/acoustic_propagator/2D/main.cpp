
/*

mencoder 'mf://*.png' -mf  fps=25:type=png -ovc copy -oac copy -o output.avi

*/

// #define FFTW_PLAN_MODE FFTW_PATIENT

#define N_FFT_THREADS  3

/*
#include <mathlib/math/std_math.h>

#include <mathlib/math/grids/grid2D.h>
#include <mathlib/math/grids/grid_file.cpp>
#include <mathlib/math/grids/extra/plots.cpp>

#include <mathlib/math/PDE/laplacian_hdaf.h>

#include <mathlib/link.cpp>./a.
#include <mathlib/non-lib_things.h>
*/


#include <mathlib/math/grids/grid2D.h>
#include <mathlib/math/grids/grid_file.cpp>
#include <mathlib/math/grids/extra/plots.cpp>

//#include <mathlib/math/PDE/laplacian_hdaf.h>
#include "/workspace/waveprop_lib/waveprop.h"
#include <mathlib/non-lib_things.h>


#include <mathlib/link.cpp>

void output(double t, grid2D< > & u, grid2D< > & v)
{
	static int n = 0;
    sprintf(fname,"/workspace/output/acoustic_propagate_2d/out_dat/u_%05d.dat", n);
    
	assert(!writeFile(u,fname));
	n++;
}

ml_color cmap_damp(double x)
{
    x *= 10;
	float s = atan(x)/pi+0.5;
	return ml_color(s,s,s);
}

/*
class acoustic_propagator
{
public:
    grid2D<double, double > * U;
    grid2D<double, double > * V;
    grid2D<double, double > C2;
    grid2D<double, double > damping;
    int exp_order;
    laplacian_2d_hdaf Del2;
    
    acoustic_propagator( ):exp_order(0),U(0),V(0),Del2(),C2(),damping() {}
    
    ~acoustic_propagator() {
        if (U) delete [] U;
        if (V) delete [] V;
        damping=0;
        U=0; V=0; exp_order=0; }
    
public:
    
    void init( grid2D<> & grid, int new_order, grid2D<> & C2_, grid2D<> & damping_ )
    {
        exp_order = new_order;
            
        if (U) delete [] U;
        if (V) delete [] V;
        U = new grid2D<> [exp_order+2];
        V = new grid2D<> [exp_order+1];
        
        for (int k=0; k<=exp_order+1; k++ )
            U[k] = grid;
        
        for (int k=0; k<=exp_order; k++ )
            V[k] = grid;
        
        int n1 = grid.n1;
        int n2 = grid.n2;
        double L1 = grid.b1-grid.a1;
        double L2 = grid.b2-grid.a2;
        
        Del2.init( n1, n2, L1, L2, 8, 8,  0.8, 0.8 );
        
        C2 = C2_;
        damping = damping_;
        
        
        sprintf(fname,"/workspace/research/acoustic_propagator/2D/damp.png" );
        plotGrid2D_1(damping,fname,cmap_damp);
    }
    
    void operator() ( double t, grid2D<> & u0, grid2D<> & v0, grid2D<> & u1, grid2D<> & v1, bool damp = true )
    {
        propagate( t, u0, v0, u1, v1, damp );
    }
    
    void propagate( double t, grid2D<> & u0, grid2D<> & v0, grid2D<> & u1, grid2D<> & v1, bool damp = true )
    {
     //   u1.copy(u0);
     //   v1.copy(v0);
        
        U[0] = u0;
        V[0] = v0;
        
        for (int n = 1; n<=exp_order+1; n++)
        {
            Del2.execute( U[n-1].array, U[n].array );
            
            U[n] *= C2;
            U[n] *= t*t/(4*n*n-2*n);
        }
        
        for (int n = 1; n<=exp_order; n++)
        {
            Del2.execute( V[n-1].array, V[n].array );
            
            V[n] *= C2;
            V[n] *= t*t/(4*n*n-2*n);
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
        
        if ( damp )
        {
            for (int i=0; i<u1.n1; i++ )
            for (int j=0; j<u1.n2; j++ )
                u1(i,j) *= exp(-t*damping(i,j));
        }
        
    }
};
*/


class u_0_func : public functor2<double,double >
{	
public:
	double operator() (double const & x, double const & y) const 
	{
    //    double k1 = 1.0;
    //    double k2 = 0.0;
    //    return cos( k1*x2+k2*x2 );
        
    //      return exp(-x1*x1)*exp(-x2*x2);
        
        double x1 = x - 7.0;
        double x2 = y + 0;
        
        
		double s = sqrt(x1*x1+x2*x2);
        double a = 0.5;
        
		if (s <= 0.5*a)
			return cos(ml_2pi*s/a)+1.0;
		else return 0.0;
	}
};

class C2_func : public functor2<double,double >
{	
public:
	double operator() (double const & x, double const & y) const 
	{
        return 10.0;
        
        double xx = x - 5;
        double yy = y;
    
        return 1.0/( 2.0 + xx*xx )+ 3.0;
        
        
	}
};

class damping_func : public functor2<double,double >
{	
public:
	double operator() (double const & x, double const & y) const 
	{
        double ax = 0.4;
        double ay = 0.4;        
        
        return  1500.0*( 1.0/(cosh((x-10)/ax)) +  1.0/(cosh((x+10)/ax)) + 1.0/(cosh((y-10)/ay)) +  1.0/(cosh((y+10)/ax))  );
        //return  5.0/(cosh((x-10)/ax));
        // +cosh((y-10)/ay)+cosh((y+10)/ay)
	}
};






int main()
{
    std_setup();
    
    int n = 256;
    double x0 = 10;
    
    double tf = 10.0;
    double dt = 0.05;
    int expansion_order = 7;
    int hdaf_order = 8;
    
    grid2D<double,double,double > grid( n,-10,10.0, n,-10,10.0 ),u,v,C2,damping;
    grid = 0.0;
    
    u = grid;
    v = grid;
    C2 = grid;
    damping = grid;
    
    u = u_0_func();
    C2 = C2_func();
    damping = damping_func();
    
    sprintf(fname,"/workspace/output/acoustic_propagate_2d/out_dat/C2.dat" );
    writeFile(C2,fname);
    
    
    void * pdata = 0;
    method1_init( &pdata, grid.n1, grid.n2, grid.dx1(), grid.dx2(), C2.array, damping.array, 7, hdaf_order, hdaf_order, 0.8, 0.8 );
    
    //acoustic_propagator P;
    //P.init( grid, 7, C2, damping );
    
    double t = 0.0;
    
    cout << "step: " << t << "\t" << L2norm(u) << endl;
    
    n=0;
    double * b_times = new double [ int(tf/dt+5) ];
    
    while (t <= tf)
    {
    
        output(t,u,v);
        
        double t1 = get_real_time();
        
        method1_execute( pdata, dt, u.array,v.array,u.array,v.array );
        //P( dt, u,v,u,v );
        
        double t2 = get_real_time();
        
        b_times[n] = t2-t1;
        double mean = 0;
        for (int j=0; j<n; j++)
            mean += b_times[j]/n;
        cout << mean << endl;
        n++;
        
        
        double mag = L2norm(u);
        cout << "step: " << t << "\t" << t2-t1 << "\t" << mag << endl;
        
        t += dt;
        
        if ( mag > 1E10 ) break;
    }
       
    output(t,u,v);
}




/*

void output(double t, grid2D< > & u, grid2D< > & v)
{
	static int n = 0;
    sprintf(fname,"/workspace/output/acoustic_propagate_2d/out_dat/u_%05d.dat", n);
    
	assert(!writeFile(u,fname));
	n++;
}

ml_color cmap_damp(double x)
{
    x *= 10;
	float s = atan(x)/pi+0.5;
	return ml_color(s,s,s);
}

class acoustic_propagator
{
public:
    grid2D<double, double > * U;
    grid2D<double, double > * V;
    grid2D<double, double > C2;
    grid2D<double, double > damping;
    int exp_order;
    laplacian_2d_hdaf Del2;
    
    acoustic_propagator( ):exp_order(0),U(0),V(0),Del2(),C2(),damping() {}
    
    ~acoustic_propagator() {
        if (U) delete [] U;
        if (V) delete [] V;
        damping=0;
        U=0; V=0; exp_order=0; }
    
public:
    
    void init( grid2D<> & grid, int new_order, grid2D<> & C2_, grid2D<> & damping_ )
    {
        exp_order = new_order;
            
        if (U) delete [] U;
        if (V) delete [] V;
        U = new grid2D<> [exp_order+2];
        V = new grid2D<> [exp_order+1];
        
        for (int k=0; k<=exp_order+1; k++ )
            U[k] = grid;
        
        for (int k=0; k<=exp_order; k++ )
            V[k] = grid;
        
        int n1 = grid.n1;
        int n2 = grid.n2;
        double L1 = grid.b1-grid.a1;
        double L2 = grid.b2-grid.a2;
        
        Del2.init( n1, n2, L1, L2, 8, 8,  0.8, 0.8 );
        
        C2 = C2_;
        damping = damping_;
        
        
        sprintf(fname,"/workspace/research/acoustic_propagator/2D/damp.png" );
        plotGrid2D_1(damping,fname,cmap_damp);
    }
    
    void operator() ( double t, grid2D<> & u0, grid2D<> & v0, grid2D<> & u1, grid2D<> & v1, bool damp = true )
    {
        propagate( t, u0, v0, u1, v1, damp );
    }
    
    void propagate( double t, grid2D<> & u0, grid2D<> & v0, grid2D<> & u1, grid2D<> & v1, bool damp = true )
    {
        u1.copy(u0);
        v1.copy(v0);
        
        U[0] = u0;
        V[0] = v0;
        
        for (int n = 1; n<=exp_order+1; n++)
        {
            Del2.execute( U[n-1].array, U[n].array );
            
            U[n] *= C2;
            U[n] *= t*t/(4*n*n-2*n);
        }
        
        for (int n = 1; n<=exp_order; n++)
        {
            Del2.execute( V[n-1].array, V[n].array );
            
            V[n] *= C2;
            V[n] *= t*t/(4*n*n-2*n);
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
        
        if ( damp )
        {
            for (int i=0; i<u1.n1; i++ )
            for (int j=0; j<u1.n2; j++ )
                u1(i,j) *= exp(-t*damping(i,j));
        }
        
    }
};


class u_0_func : public functor2<double,double >
{	
public:
	double operator() (double const & x, double const & y) const 
	{
    //    double k1 = 1.0;
    //    double k2 = 0.0;
    //    return cos( k1*x2+k2*x2 );
        
    //      return exp(-x1*x1)*exp(-x2*x2);
        
        double x1 = x - 7.0;
        double x2 = y + 0;
        
        
		double s = sqrt(x1*x1+x2*x2);
        double a = 0.5;
        
		if (s <= 0.5*a)
			return cos(ml_2pi*s/a)+1.0;
		else return 0.0;
	}
};

class C2_func : public functor2<double,double >
{	
public:
	double operator() (double const & x, double const & y) const 
	{
        return 10.0;
        
        double xx = x - 5;
        double yy = y;
    
        return 1.0/( 2.0 + xx*xx )+ 3.0;
        
        
	}
};

class damping_func : public functor2<double,double >
{	
public:
	double operator() (double const & x, double const & y) const 
	{
        double ax = 0.4;
        double ay = 0.4;        
        
        return  1500.0*( 1.0/(cosh((x-10)/ax)) +  1.0/(cosh((x+10)/ax)) + 1.0/(cosh((y-10)/ay)) +  1.0/(cosh((y+10)/ax))  );
        //return  5.0/(cosh((x-10)/ax));
        // +cosh((y-10)/ay)+cosh((y+10)/ay)
	}
};

*/








        
/*        if ( -9 <= x and x <= -6   and   -9 <= y and y <= -6   )
            return 2.0;
        
        if ( 0 <= x and x <= 8   and   -2 <= y and y <= 3   )
            return 3.0;
        
        if ( 0 <= x and x <= 8   and   4 <= y and y <= 5   )
            return 6.0;*/
        
        
 /*       if (  18.0 <= y and y <= 19.5 )
            return 7.0;
        
        if (  9 <= x and y <= 18.0 and y >= 12.75+0.25*x )
            return 9.0;*/
        
       // return 4.0;

















/*




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
        
        
        /*
	}
};









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


