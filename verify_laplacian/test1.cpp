


#include <mathlib/math/std_math.h>
#include <mathlib/math/random/ml_random.h>


#include "common.h"



void rand_tpoly_2d( cgrid2D & G, int k1_max, int k2_max )
{
 /*   ml_random rng;
    
    for ( int k1=0; k1<G.n1; k1++ )
    for ( int k2=0; k2<G.n2; k2++ )
        if ( (k1 <= k1_max ) and k2<k2_max )
            G(k1,k2) = exp(iu*_2pi*rng.gen_double())*rng.gen_double();
        else
            Q[0][k] = 0;*/
}



void test1_2d( void (*del2)( double *, double *, int, int ), int n1, int n2, int n_powers, int gamma )
{
    cgrid2D cgrid( n1, 0,1, n2, 0,1 );
    
    cgrid2D * Q = ml_alloc<cgrid2D > ( n_powers+1 );
    
    for ( int k=0; k<=n_powers; k++ )
        Q[k].copy( cgrid );
    
    //rand_tpoly_2d( Q[0], gamma*n1, gamma*n2 );
    
    
    
    
    /*
    {
        
        
        for ( int k1=0; k1<n1; k1++ )
        for ( int k2=0; k2<n2/2+1; k2++ )
            if ( k1 < k1_max and k2<k2_max )
                ary2( Q[0], k1, k2, n1, n2/2+1 );
                Q[0][] = exp(iu*_2pi*rng.gen_double())*rng.gen_double();
            else
                Q[0][k] = 0;
    
    for (int j=1; j<n_d; j++)
    for (int k=0; k<n/2+1; k++)
        Q[j][k] = Q[0][k]*pow( iu*((double)k)*_2pi/(b-a), j ) ;
	
    for (int j=0; j<n_d; j++)
        IFFT( Q[j], D[j], n );
        
    for (int j=0; j<n_d; j++)
    for (int k=0; k<n; k++)
        D[j][k] /= n;
    
    ml_free(Q, n_d);*/
}







/*






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
    
	*/
