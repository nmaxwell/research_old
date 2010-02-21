
#include <mathlib/math/std_math.h>
#include <mathlib/tools/graphing/plot2D.h>

#include <mathlib/math/random/ml_random.h>
#include <mathlib/math/grids/grid2D.h>
#include <mathlib/math/grids/extra/plots.cpp>
#include <mathlib/math/PDE/laplacian_hdaf.h>

#include <mathlib/link.cpp>
#include <mathlib/non-lib_things.h>

#include "Polygon2D.cpp"





bool straddle_test( float const & x1, float const & x2, float const & x3, float const & x4, float const & y1, float const & y2, float const & y3, float const & y4 )
{
    if ( max(x1,x2) < min(x3,x4) )
        return false;
    if ( max(x3,x4) < min(x1,x2) )
        return false;
    if ( max(y1,y2) < min(y3,y4) )
        return false;
    if ( max(y3,y4) < min(y1,y2) )
        return false;
    
    float z1 = (x3-x1)*(y2-y1)-(y3-y1)*(x2-x1);
    float z2 = (x4-x1)*(y2-y1)-(y4-y1)*(x2-x1);
    
    //if ( z1==0.f && z2==0.f ) return true;
    
    if ( (z1<0) ^ (z2<0) ) return true;
    
    return false;
}

int hitting_index( Polygon2D &boundary, float X1_0, float X2_0, double dt )
{
    static ml_random rng;
    
    float W1 = X1_0;
    float W2 = X2_0;
    float dW1,dW2;
    
    float sqrt_dt = sqrt(dt);
    int max_steps = 10000;
    int step=0;
    
    while(true )
    {
        step++;
        if(step>max_steps)
            return -1;
        
        rng.std_normal_rv(dW1,dW2);
        
        if ( !boundary.interior_test(W1+dW1, W2+dW2) )
            break;
        
        W1 += dW1*sqrt_dt;
        W2 += dW2*sqrt_dt;
    }
    
    for (int k=0; k<boundary.n_points-1; k++ )
    {
        if ( straddle_test(
            boundary.x_point(k), boundary.x_point(k+1), W1, W1+dW1,
            boundary.y_point(k), boundary.y_point(k+1), W2, W2+dW2 ) )
        {
            return k;
        }
    }
    
    return -1;
}

int fails = 0;

double solve_laplace( double x, double y, Polygon2D &boundary, double *f, double dt, int n_trials )
{
    int n_bd = boundary.n_points-1;
    int trial = 0;
    double sum = 0;
    
    while ( trial < n_trials )
    {
        int index = hitting_index( boundary, x, y, 0.1 );
        if ( index >=0 and index<n_bd )
        {
            trial++;
            sum += f[index];
        }
        else
            fails += 1;
    }
    
    return sum/n_trials;
}



float boundary_def(float & x, float & y, double & f, double t)
{
	// t \in [0,1]
	
	x = cos(t*_2pi)*(cos(t*_2pi*5)*2.5+6);
	y = sin(t*_2pi)*(cos(t*_2pi*5)*2.5+6);
	
	f = cos(_2pi*t*15)*2.0 + sin(_2pi*t*13)*3.0+0*sin(_2pi*t*60)*7.0 ;
    f *= 2.0;
}

ml_color color_map_green( double z )
{
	float x = atan((double)z)/pi+0.5;
	return ml_color(0,x,0);
}

ml_color color_map_blue( double z )
{
	float x = atan((double)z)/pi+0.5;
	return ml_color(0,0,x);
}

ml_color color_map_red( double z )
{
	float x = atan((double)z)/pi+0.5;
	return ml_color(x,0,0);
}

int main()
{
    std_setup();
    
    
    
    
    int n_bd=150;
    Polygon2D boundary;
    double * f = 0;
    
    if (1)
    {
		// setup boundary
		
		//boundary.set_inscribed_circle(3, 7);
		
		f = ml_alloc<double> (n_bd);
		boundary.init(n_bd+1);
		
		for (int k=0; k<n_bd; k++)
		{
			double t1 = double(k)/n_bd;
			double t2 = double(k+1)/n_bd;
			
			boundary_def( boundary.x_point(k),boundary.y_point(k), f[k], t1 );
			boundary_def( boundary.x_point(k+1),boundary.y_point(k+1), f[k+1], t2 );
		}
	}
    else
    {
        n_bd = 4;
        
        f = ml_alloc<double> (n_bd);
        boundary.init(n_bd+1);
        
        f[0] = 0.0;
        f[1] = 7.0;
        f[2] = 0.0;
        f[3] = 0.0;
        
        boundary.x_point(0) = -8.0;
        boundary.y_point(0) = +8.0;
        
        boundary.x_point(1) = +8.0;
        boundary.y_point(1) = +8.0;

        boundary.x_point(2) = +8.0;
        boundary.y_point(2) = -8.0;
        
        boundary.x_point(3) = -8.0;
        boundary.y_point(3) = -8.0;
        
        boundary.x_point(4) = -8.0;
        boundary.y_point(4) = +8.0;
    }
    
    
    
    
    
    plot2D plt( -10.0, 10.0, -10.0, 10.0, 512, 512 );
    plt = ml_white;
    for (int k=0; k<n_bd; k++)
        plt.ptLineThick( boundary.x_point(k),boundary.y_point(k), boundary.x_point(k+1), boundary.y_point(k+1), 5, color_map_red(f[k]) );
    
    for (int i=0; i<plt.NX; i++)
    for (int j=0; j<plt.NY; j++)
    {
        cout << i << "\t " << j << endl;
        double x = plt.to_x(i);
        double y = plt.to_y(j);
        
        if (boundary.interior_test(x,y))
        {
            double u = solve_laplace( x, y, boundary, f, 0.5, 200 );
            
            float s = atan( u )/pi+0.5;
            plt.set_px(i,j, ml_color(0.f,0.f,s) );
        }
    }
    
    plt.png("/workspace/output/scratch/result.png");
        
    cout << "fails: " << fails << endl;
    
    std_exit();
}


