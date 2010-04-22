
#include <mathlib/math/std_math.h>

#include <mathlib/math/grids/grid_file.h>
#include <mathlib/math/grids/grid2D.h>
#include <mathlib/math/grids/extra/plots.cpp>


#include <mathlib/link.cpp>
#include <mathlib/non-lib_things.h>







int main()
{
    std_setup();
    
    
	sprintf(fname,"%s/C_%d%d.bmp", pwd,i,j );
	cout << fname << endl;
	
	bmp_data bmp;
	assert( !read_bmp( fname, bmp ) );
		
	grid2D<double, double > grid( bmp.nx, 0,1, bmp.ny, 0,1 );
    
	for (int k = 0; k<bmp.nx*bmp.ny; k++ )
		grid[k] = ((double)bmp.data[k])/(256.0);
		
	sprintf(fname,"%s/C_%d%d.dat",pwd,i,j );
		
	assert(!writeFile(grid,fname));
    
    
    
    
    
    std_exit();
}




