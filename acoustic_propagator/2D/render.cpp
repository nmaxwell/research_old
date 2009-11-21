		
#include <mathlib/math/grids/grid2D.h>
#include <mathlib/math/grids/extra/plots.cpp>
#include <mathlib/math/grids/grid_file.cpp>

#include <mathlib/link.cpp>

color3 cmap(double x)
{
	x *= 10;
	float s = atan(x)/pi+0.5;
	return color3(s,s,s);
	
}

int main()
{
	std_setup();
	
	grid2D<> G;
	
	sprintf(fname,"%s/out_dat/C2.dat",pwd);
	readFile(G,fname);
	sprintf(fname,"%s/out_png/C2.png",pwd);
	plotGrid2D_1(G,fname,cmap);
	
	for( int n = 0; n<9000; n+=1)
	{
		cout << n << endl;
		sprintf(fname,"%s/out_dat/u_%05d.dat",pwd,n);
		if ( readFile(G,fname) ) break;
		
		sprintf(fname,"%s/out_png/%05d.png",pwd,n);	    
		plotGrid2D_1(G,fname,cmap);
	}	
}

