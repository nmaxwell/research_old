

//#define FFTW_PLAN_MODE FFTW_PATIENT

#define N_FFT_THREADS  1




#include <mathlib/math/std_math.h>
#include <mathlib/math/grids/grid2D.h>
#include <mathlib/math/transforms/fft.h>
#include <mathlib/math/laplacian/laplacian.h>
#include <mathlib/math/laplacian/laplacian_FD.h>
#include <mathlib/math/laplacian/laplacian_hdaf.h>

#include <mathlib/math/grids/extra/plots.cpp>

#include "common.h"

#include <mathlib/non-lib_things.h>

#include <mathlib/link.cpp>



ml_color cmap(double x)
{
	x *= 10;
	float s = atan(x)/pi+0.5;
	return ml_color(s,s,s);
}

double gauss( double x, double y )
{
    return exp(-x*x)*exp(-y*y);
}

double del2_gauss( double x, double y )
{
    double r2 = x*x+y*y;
    return 4.0*(r2-1.0)*exp(-x*x)*exp(-y*y);
}

double del4_gauss( double x, double y )
{
    double r2 = x*x+y*y;
    return 16.0*(r2*(r2-4.0)+2.0)*exp(-x*x)*exp(-y*y);
}

double del6_gauss( double x, double y )
{
    double r2 = x*x+y*y;
    return 64.0*( r2*(r2*(r2-9.0)+18.0) -6.0 )*exp(-x*x)*exp(-y*y);
}


int main()
{
    std_setup();
    
    double X0 = 40;
    
    int n1 = 512;
    int n2 = n1;
    double a1 = -X0;
    double b1 = X0;
    double a2 = -X0;
    double b2 = X0;
    
    rgrid2D grid( n1, a1, b1,  n2, a2, b2 );
    
    rgrid2D A0(grid), A1(grid), A2(grid), A3(grid);
    rgrid2D N1(grid), N2(grid), N3(grid);
    
    A0 = gauss;
    A1 = del2_gauss;
    A2 = del4_gauss;
    A3 = del6_gauss;
    
    laplacian_2d_fd Del2;
    Del2.init( n1, n2, b1-a1, b2-a2, 24, 24 );
    
  //  laplacian_2d_hdaf Del2;
  //  Del2.init( n1, n2, b1-a1, b2-a2, 24, 24,  0.9, 0.9 );
    
    Del2.execute( A0.array, N1.array );
    Del2.execute( N1.array, N2.array );
    Del2.execute( N2.array, N3.array );
    
    //laplacian_2d_fft( A0.array, N1.array, n1, n2, b1-a1, b2-a2 );
    //laplacian_2d_fft( N1.array, N2.array, n1, n2, b1-a1, b2-a2 );
    //laplacian_2d_fft( N2.array, N3.array, n1, n2, b1-a1, b2-a2 );
    
    cout << endl;
    cout << l2_error( A1.array, N1.array, n1, n2 ) << endl;
    cout << l2_error( A2.array, N2.array, n1, n2 ) << endl;
    cout << l2_error( A3.array, N3.array, n1, n2 ) << endl;
    
    sprintf(fname, "/workspace/output/temp/out1.png" );
    plotGrid2D_1(A0,fname,cmap);
    sprintf(fname, "/workspace/output/temp/out2.png" );
    plotGrid2D_1(N1,fname,cmap);
    sprintf(fname, "/workspace/output/temp/out3.png" );
    plotGrid2D_1(N2,fname,cmap);
    sprintf(fname, "/workspace/output/temp/out4.png" );
    plotGrid2D_1(N3,fname,cmap);
    
    
    
    
    
    for (int j=2; j<100; j++)
    {
        Del2.execute( N1.array, N2.array );
        N1 = N2;
        
        cout << j << "\t" << log10(fabs(max(N2.array,n1*n2))) << endl;
        //-log10(dfactorial[j*2])
    }
    
    
    
    
    
    
    
    
    
    
    
    double t1,t2;
    int n=10;
    
    t1 = get_real_time();
    
    for (int j=0; j<n; j++)
    {
        //laplacian_2d_fft( A0.array, N1.array, n1, n2, b1-a1, b2-a2 );
        Del2.execute( A0.array, N1.array );
    }
    
    t2 = get_real_time();
    
    cout << (t2-t1)/n << endl;
    
    std_exit();
}








