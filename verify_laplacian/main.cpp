
#include <mathlib/math/std_math.h>
#include <mathlib/math/grids/grid2D.h>
#include <mathlib/math/transforms/fft.h>
#include <mathlib/math/laplacian/laplacian.h>
#include <mathlib/math/grids/extra/plots.cpp>

#include "test1.cpp"

#include <mathlib/non-lib_things.h>

#include <mathlib/link.cpp>



ml_color cmap(double x)
{
	x *= 10;
	float s = atan(x)/pi+0.5;
	return ml_color(s,s,s);
}


double gauss(  double x, double y )
{
    return exp(-x*x)*exp(-y*y);
}

double del2_gauss(  double x, double y )
{
    return (4.0*x*x+4.0*y*y-4.0)*gauss(x,y);
}




int main()
{
    std_setup();
    
    int n1 = 512;
    int n2 = 512;
    int a1 = -10;
    int b1 = 10;
    int a2 = -10;
    int b2 = 10;
    
    rgrid2D grid( n1, a1, b1,  n2, a2, b2 );
    
    rgrid2D G0(grid), G1(grid), G2(grid);
    
    G0 = gauss;
    G1 = del2_gauss;
    
    laplacian( G0.array, G2.array, n1, n2, b1-a1, b2-a2 );
    
    
    sprintf(fname, "/workspace/output/temp/out1.png" );
    plotGrid2D_1(G0,fname,cmap);
    sprintf(fname, "/workspace/output/temp/out2.png" );
    plotGrid2D_1(G1,fname,cmap);
    sprintf(fname, "/workspace/output/temp/out3.png" );
    plotGrid2D_1(G2,fname,cmap);
    
    
    
    
    std_exit();
}








