

#define INCLUDE_FD_D2_LR
#define N_FFT_THREADS 1


#include <mathlib/math/grids/grid2D.h>
#include <mathlib/math/grids/grid_file.cpp>
#include <mathlib/math/grids/extra/plots.cpp>

#include "../Del2_hdaf.cpp"

#include <mathlib/math/grids/grid2D_Del2_FD.h>



#include <mathlib/link.cpp>

color3 cmap(double x)
{
	x *= 10;
	float s = atan(x)/pi+0.5;
	return color3(s,s,s);
}




double kx = 5.0*pi;
double ky = 2.0*pi;


class test_func : public functor2<double,double >
{	
public:
	double operator() (double const & x, double const & y) const 
	{
       // return cos(x*kx);
        return real( exp(iu*(x*kx+y*ky)) );
	}
};

class Del2_test_func : public functor2<double,double >
{	
public:
	double operator() (double const & x, double const & y) const 
	{
       // return -kx*kx*cos(x*kx);
        return real( -exp(iu*(x*kx+y*ky))*(kx*kx+ky*ky) );
	}
};



int main()
{
    std_setup();
    
    int N = 200;
    
    double xa = -1.0;
    double xb = +1.0;
    
    double ya = -1.0;
    double yb = +1.0;
    
    grid2D<double > grid(N,xa,xb, N,xa,xb);
    grid = 0.0;
    
    grid2D<double > G1,G2,G3,err;
    err = grid;
    G1 = grid;
    G2 = grid;
    G3 = grid;
    
    for (double c = 10.0; c<=20; c+=1.0)
    {
        G1 = test_func();
        G3 = Del2_test_func();
        
        
        Del2_hdaf(G1,G2, pi*c, 16,16 );
        //Del2_FD(G1, G2, 16, GRID_NONPER_BD_1, 24 ); 
        
        err = G3;
        err -= G2;
        
        cout << c << "    " << L2norm(err) << endl;
    }
    
    if(0)
    {
        sprintf(fname,"%s/out/G1.png",pwd);
        plotGrid2D_1(G1,fname,cmap);
        
        sprintf(fname,"%s/out/G2.png",pwd);
        plotGrid2D_1(G2,fname,cmap);
        
        sprintf(fname,"%s/out/err.png",pwd);
        plotGrid2D_1(err,fname,cmap);
    }
}
    
	
