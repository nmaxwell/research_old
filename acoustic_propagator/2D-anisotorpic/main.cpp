
/*

mencoder 'mf://*.png' -mf  fps=25:type=png -ovc copy -oac copy -o output.avi

*/

// #define FFTW_PLAN_MODE FFTW_PATIENT

#define N_FFT_THREADS  1

#include <mathlib/math/std_math.h>

#include <mathlib/math/grids/grid2D.h>
#include <mathlib/math/grids/grid_file.cpp>
#include <mathlib/math/grids/extra/plots.cpp>

#include <mathlib/math/PDE/grad_hdaf.h>

#include <mathlib/link.cpp>
#include <mathlib/non-lib_things.h>


typedef grid2D<double, double > grid;

struct grid_vec
{
    grid _1;
    grid _2;
    
    void copy( grid_vec & rhs )
    {
        _1.copy( rhs._1 );
        _2.copy( rhs._2 );
    }
    
    void copy_data( grid_vec & rhs )
    {
        _1.copyData( rhs._1 );
        _2.copyData( rhs._2 );
    }
    
    void operator= ( grid_vec & rhs )
    {
        this->copy_data( rhs );
    }
    
    void operator *= ( double rhs )
    {
        _1 *= rhs;
        _2 *= rhs;
    }
    
    void operator = ( double rhs )
    {
        _1 = rhs;
        _2 = rhs;
    }
    
    void operator *= ( grid &rhs )
    {
        _1 *= rhs;
        _2 *= rhs;
    }
};

struct grid_mat
{
    grid _11;
    grid _12;
    grid _21;
    grid _22;
    
};

struct grid_symmat
{
    grid _11;
    grid _12;
    grid _22;    
    
};

class strain_tensor
{
public:
    grad_2d_hdaf grad;
    
public:
    
    void init( grid & g, int m1, int m2, double gam1, double gam2 )
    {
        grad.init( g.n1, g.n2, g.b1-g.a1, g.b2-g.a2, m1, m2, gam1, gam2 );
    }
    
    void execute( grid_vec & u, grid_symmat & del_u )
    {
        del_u._11.copy(u._1);
        del_u._12.copy(u._1);
        del_u._22.copy(u._1);
        
        grad.linear_strain_tensor( u._1.array, u._2.array, del_u._11.array, del_u._12.array, del_u._22.array );
    }
};


class stress_tensor
{
public:
    strain_tensor del;
    grid_symmat eps;
    grid stiffness_tensor[6];
    
    //    stiffness_tensor[0] = C_1111
    //    stiffness_tensor[1] = C_1122
    //    stiffness_tensor[2] = C_1112
    //    stiffness_tensor[3] = C_1222
    //    stiffness_tensor[4] = C_2212
    //    stiffness_tensor[5] = C_1212
    
public:
    
    void init( grid & g, grid * C, int m1, int m2, double gam1, double gam2 )
    {
        del.init( g, m1, m2, gam1, gam2 );
        
        for (int k=0; k<6; k++ )
            stiffness_tensor[k] = C[k];
    }
    
    void execute( grid_vec & u, grid_symmat & sig )
    {
        // sig = C : eps
        
        int N = u._1.n1*u._1.n2;
        
        sig._11.copy(u._1);
        sig._12.copy(u._1);
        sig._22.copy(u._1);
        
        del.execute( u, eps );
        
        for ( int n=0; n<N; n++ )
            sig._11[n] =
                stiffness_tensor[0][n]*eps._11[n] +
                stiffness_tensor[1][n]*eps._22[n] +
                stiffness_tensor[2][n]*eps._12[n]*2;
        
        for ( int n=0; n<N; n++ )
            sig._22[n] =
                stiffness_tensor[1][n]*eps._11[n] +
                stiffness_tensor[3][n]*eps._22[n] +
                stiffness_tensor[4][n]*eps._12[n]*2;
        
        for ( int n=0; n<N; n++ )
            sig._12[n] =
                stiffness_tensor[2][n]*eps._11[n] +
                stiffness_tensor[4][n]*eps._22[n] +
                stiffness_tensor[5][n]*eps._12[n]*2;
    }
};


class anis_diff_op
{
public:
    stress_tensor del_sig;
    grad_2d_hdaf del;
    grid_symmat sig;
    
public:
    
    void init( grid & g, grid * C, int m1, int m2, double gam1, double gam2 )
    {
        del_sig.init( g, C, m1, m2, gam1, gam2  );
        del.init( g.n1, g.n2, g.b1-g.a1, g.b2-g.a2, m1, m2, gam1, gam2 );
        
        sig._11.copy( g);
        sig._12.copy( g);
        sig._22.copy( g);
    }
    
    void execute( grid_vec & u, grid_vec & del_u )
    {
        del_sig.execute( u, sig );
        del.grad_dot_sym( sig._11.array, sig._12.array, sig._22.array, del_u._1.array, del_u._2.array );
    }
};


class anisotropic_propagator
{
public:
    grid_vec * U;
    grid_vec * V;
    int exp_order;
    anis_diff_op Del;
    
    anisotropic_propagator( ):exp_order(0),U(0),V(0),Del() {}
    
    ~anisotropic_propagator() {
        if (U) delete [] U;
        if (V) delete [] V;
        U=0; V=0; exp_order=0; }
    
public:
    
    void init( grid2D<> & g, int new_order, grid * C, int m1, int m2, double gam1, double gam2 )
    {
        exp_order = new_order;
        
        if (U) delete [] U;
        if (V) delete [] V;
        U = new grid_vec [exp_order+2];
        V = new grid_vec [exp_order+1];
        
        for (int k=0; k<exp_order+2; k++ )
        {
            U[k]._1 = g;
            U[k]._2 = g;
        }
        
        for (int k=0; k<exp_order+1; k++ )
        {
            V[k]._1 = g;
            V[k]._2 = g;
        }
        
        int n1 = g.n1;
        int n2 = g.n2;
        double L1 = g.b1-g.a1;
        double L2 = g.b2-g.a2;
        
        Del.init( g, C, m1, m2, gam1, gam2  );
    }
    
    void operator() ( double t, grid_vec & u0, grid_vec & v0, grid_vec & u1, grid_vec & v1 )
    {
        propagate( t, u0, v0, u1, v1 );
    }
    
    void propagate ( double t, grid_vec & u0, grid_vec & v0, grid_vec & u1, grid_vec & v1 )
    {
        u1.copy(u0);
        v1.copy(v0);
        
        U[0] = u0;
        V[0] = v0;
        
        for (int n = 1; n<=exp_order+1; n++)
        {
            Del.execute( U[n-1], U[n] );
            
            U[n] *= t*t/(4*n*n-2*n);
        }
        
        for (int n = 1; n<=exp_order; n++)
        {
            Del.execute( V[n-1], V[n] );
            
            V[n] *= t*t/(4*n*n-2*n);
        }
        
        
        u1 = 0.0;
        v1 = 0.0;
        
        for (int n = 0; n<=exp_order; n++)
        GRID2D_LINLOOP( u0._1 )
        {
            u1._1[i] += U[n]._1[i];
            u1._1[i] += V[n]._1[i]* t/(2.0*n+1.0);
            
            v1._1[i] += U[n+1]._1[i]*((4.0*n+6.0)*n+2.0)/(t*(2.0*n+1.0));
            v1._1[i] += V[n]._1[i];
            
            
            u1._2[i] += U[n]._2[i];
            u1._2[i] += V[n]._2[i]* t/(2.0*n+1.0);
            
            v1._2[i] += U[n+1]._2[i]*((4.0*n+6.0)*n+2.0)/(t*(2.0*n+1.0));
            v1._2[i] += V[n]._2[i];
        }
    }
};







void output(double t, grid_vec & u, grid_vec & v)
{
	static int n = 0;
    
    sprintf(fname,"/workspace/output/anisotropic_propagate_2d/out_dat/u1_%05d.dat", n);
	assert(!writeFile(u._1,fname));
    
    sprintf(fname,"/workspace/output/anisotropic_propagate_2d/out_dat/u2_%05d.dat", n);
	assert(!writeFile(u._2,fname));
    
	n++;
}


class u1_func : public functor2<double,double >
{	
public:
	double operator() (double const & x, double const & y) const 
	{
        double x1 = x - 3;
        double x2 = y + 5;
        
		double s = sqrt(x1*x1+x2*x2);
        double a = 0.5;
        
		if (s <= 0.5*a)
			return cos(ml_2pi*s/a)+1.0;
		else return 0.0;
	}
};

class u2_func : public functor2<double,double >
{	
public:
	double operator() (double const & x, double const & y) const 
	{
        double x1 = x + 3;
        double x2 = y - 5;
        
		double s = sqrt(x1*x1+x2*x2);
        double a = 0.5;
        
		if (s <= 0.5*a)
			return cos(ml_2pi*s/a)+1.0;
		else return 0.0;
	}
};






class C11_func : public functor2<double,double >
{
public:
	double operator() (double const & x, double const & y) const 
	{
        // C_1111
        
        return 2.0;
	}
};

class C12_func : public functor2<double,double >
{
public:
	double operator() (double const & x, double const & y) const 
	{
        // C_1122
        
        return 0.0;
	}
};

class C13_func : public functor2<double,double >
{
public:
	double operator() (double const & x, double const & y) const 
	{
        // C_1112
        
        return 0.0;
	}
};

class C22_func : public functor2<double,double >
{
public:
	double operator() (double const & x, double const & y) const 
	{
        // C_1222
        
        return 0.0;
	}
};

class C23_func : public functor2<double,double >
{
public:
	double operator() (double const & x, double const & y) const 
	{
        // C_2212
        
        return 0.0;
	}
};

class C33_func : public functor2<double,double >
{
public:
	double operator() (double const & x, double const & y) const 
	{
        // C_1212
        
        return 2.0;
	}
};





int main()
{
    std_setup();
    
    int n = 512;
    double x0 = 10; 
    
    double tf = 10.0;
    double dt = 0.04;
    
    grid G( n,-x0,x0, n,-x0,x0 );
    grid_vec u,v;
    
    G = 0;
    
    u._1 = G;
    u._2 = G;
    v._1 = G;
    v._2 = G;
    
    grid C[6];
    for (int k=0; k<6; k++ )
        C[k] = G;
    C[0] = C11_func();
    C[1] = C12_func();
    C[2] = C13_func();
    C[3] = C22_func();
    C[4] = C23_func();
    C[5] = C33_func();
    
    u._1 = u1_func();
    u._2 = u2_func();
    v = 0.0;
    
    //sprintf(fname,"/workspace/output/anisotropic_propagate_2d/out_dat/C.dat" );
    //writeFile(C2,fname);
    
    anisotropic_propagator P;
    P.init( G, 10, C, 8, 8, 0.6, 0.6 );
    
    double t = 0.0;
    
    while (t <= tf)
    {
        output(t,u,v);
        
        double t1 = get_real_time();
        
        P( dt, u,v, u,v );
        
        double t2 = get_real_time();
        
        double mag = L2norm(u._1) + L2norm(u._2);     
        cout << "step: " << t << "\t" << t2-t1 << "\t" << mag << endl;
        if ( mag > 1E10 ) break;
        
        t += dt;
        
        
    }
    
    output(t,u,v);
}


