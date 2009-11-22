
#include <mathlib/math/random/ml_random.h>
#include <mathlib/math/SDE/random_walk.h>
#include <mathlib/math/grids/grid1D.h>
#include <mathlib/math/grids/grid2D.h>
//#include <mathlib/math/grids/extra/plots.cpp>

#include <mathlib/link.cpp>
#include <mathlib/non-lib_things.h>

//#include <boost/math/special_functions/erf.hpp>

class gaussian : functor<double, double >
{
public:
    double mu,sig;
    
    gaussian(double mu, double sig ):mu(mu),sig(sig) {};
    
    double operator() (double const & x) const
    {
        return exp(-(x-mu)*(x-mu)/(2.0*sig*sig))/(sqrt2pi*sig);
    }
};


template< class T >
void output( T * data, int n, const char * fname, const char * delim1="\t", const char * delim2="\n" )
{
    ofstream out;
    out.open(fname, fstream::out);
    assert(out.good() && out.is_open());
    
    for ( int k=0; k<n; k++ )
        out << data[k] << delim2;
    
    out.close();
}

template< class T1, class T2 >
void output( T1 * data1, T2 * data2, int n, const char * fname, const char * delim1="\t", const char * delim2="\n" )
{
    ofstream out;
    out.open(fname, fstream::out);
    assert(out.good() && out.is_open());
    
    for ( int k=0; k<n; k++ )
        out << data1[k] << delim1 << data2[k] << delim2;
    
    out.close();
}

void picture_BM_2D( double stop_time, double time_step, grid2D<double,double > & f, int n_runs )
{
    int n_steps = stop_time/time_step;
    
    double **W=0;
    gen_BM( time_step, n_steps, W, n_runs );
    
    f = grid2D<double,double > ( n_steps, 0.0, stop_time,  n_steps, -stop_time/2, stop_time/2 );
    f = 0.0;
    
    double s = time_step*time_step*2.0;
    cout << "plotting...\n";
    for ( int n=0; n<n_runs; n++)
    for ( int k=0; k<n_steps; k++ )
    {
        for ( int j=0; j<f.n2; j++ )
            f(k,j) += exp(-(f.x2(j)-W[n][k])*(f.x2(j)-W[n][k])/s);
    }
    
}



ml_color cmap( double I)
{
    float s = I;
    return ml_color(s,s,s);
}




int main()
{
    std_setup();
    
    int n_runs = 500;
    double stop_time = 10.0;
    double h = 0.001/32;
    int n_steps = stop_time/(h*32);
    
    double **W=0;
    gen_BM( h, n_steps, W, n_runs );
    
    
        double Dt = 1.0;
        double s = h*h*1024*2.0;
        
        for (int K=0; K<stop_time/Dt; K++ )
        {
            int k = floor(((K*Dt)/(h*32)));
            double t = h*32*k;
            
            grid1D<double, double > f( 100, -3.0*sqrt(t), +3.0*sqrt(t) );
            f = 0.0;
            
            for ( int r=0; r<n_runs; r++)
            for ( int j=0; j<f.n1; j++)
                f(j) += exp(-(f.x1(j)-W[r][k])*(f.x1(j)-W[r][k])/s);
            
            double * X = ml_alloc<double > (f.n1);
            for (int k=0; k<f.n1; k++)
                X[k] = f.x1(k);
            
            sprintf(fname, "/workspace/output/X_%d", K);
            output( X, f.n1, fname );
            
            sprintf(fname, "/workspace/output/out_%d", K);
            output( f.array, f.n1, fname );
            
            ml_free(X);
        }
    
    //picture_BM_2D( stop_time, time_step, f, n_runs );
    
    //plotGrid2D_1( f,  "/workspace/output/temp/BM.png", cmap );
    
    //output( W[9], n_steps, "/workspace/output/temp/out" );
    
    std_exit();
}




