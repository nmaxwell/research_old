		
#include <mathlib/math/grids/grid2D.h>
#include <mathlib/math/grids/extra/plots.cpp>
#include <mathlib/math/grids/grid_file.cpp>

#include <mathlib/link.cpp>

#include <mathlib/non-lib_things.h>


ml_color cmap(double x)
{
	x *= 10;
	float s = atan(x)/pi+0.5;
	return ml_color(s,s,s);
	
}


ml_color cmap_C2(double x)
{
    //x -= 2800*2800;
    x *= 10;
	float s = atan(x)/pi+0.5;
	return ml_color(s,s,s);
	
}


int main()
{
	std_setup();
	
	grid2D<> G;
	
	sprintf(fname,"/workspace/output/acoustic_propagate_2d/out_dat/C2.dat" );
	readFile(G,fname);
	sprintf(fname,"/workspace/output/acoustic_propagate_2d/out_png/C2.png" );
	plotGrid2D_1(G,fname,cmap_C2);
	
    
	for( int n = 0; n<9000; n+=1)
	{
		sprintf(fname,"/workspace/output/acoustic_propagate_2d/out_dat/u_%05d.dat" ,n);
        cout << fname << endl;
		if ( readFile(G,fname) ) break;
		
		sprintf(fname,"/workspace/output/acoustic_propagate_2d/out_png/%05d.png" ,n);	    
		plotGrid2D_1(G,fname,cmap);
	}	
}

