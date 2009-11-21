
#include <mathlib/math/grids/grid2D.h>
#include <mathlib/math/grids/grid_file.h>
#include <mathlib/math/grids/grid2D_Del2_FD.h>
#include <mathlib/math/random/ml_random.h>
#include <mathlib/math/grids/extra/plots.cpp>

#include <mathlib/link.cpp>


#define DEL2( G,i,j ) ((G(i+1,j)+G(i-1,j)-2.0*G(i,j))/(G.dx1()*G.dx1())+(G(i,j+1)+G(i,j-1)-2.0*G(i,j))/(G.dx2()*G.dx2()))


#define DEl2_center_coefficient (-2.0)

color3 cmap(int x)
{
    if (x==0) return c3_white;
    
    if (x <= 9)
    {
        float s = 0.1*x;
        if (s > 1.0) s = 1.0;
            return color3(0,0,s);
    }
    
    x -= 10;
    if (x <= 9)
    {
        float s = 0.1*x;
        if (s > 1.0) s = 1.0;
            return color3(0,s,0);
    }
    
    x -= 10;
    if (x <= 9)
    {
        float s = 0.1*x;
        if (s > 1.0) s = 1.0;
            return color3(s,0,0);
    }
    return c3_red;
}

color3 cmap(double x)
{
	float s = atan(x)/pi+0.5;
	return color3(s,s,s);
}


        
template<class T >
void solve_poisson( grid2D<T,T > & phi, grid2D<T,T > & B )
{
    int n_runs = 40;
    int walk_max = (phi.n1*phi.n2);
    int n_walks = phi.n1*phi.n2;
    
    grid2D<T,T,T > sum;
    grid2D<grid2D<int,T,T> ,T,T > count;
    sum.copy(phi);
    count.copy(phi);
    
    for (int k=0; k<phi.n1; k++)
    {
        count(k,0).copy(phi);
        count(k,count.n2-1).copy(phi);
    }
    for (int k=0; k<phi.n2; k++)
    {
        count(0,k).copy(phi);
        count(count.n1-1,k).copy(phi);
    }
    
    int i[walk_max+1];
    int j[walk_max+1];
    
    int n_series = 2*walk_max+2;
    int n_series32 = div_up(n_series,32);
    bool series[n_series];
    int series32[n_series32];
    ml_random rng(n_series32);
    
    double time0 = get_real_time();
    
    for (int i0=1; i0<phi.n1-1; i0++)
    {
        cout << i0 << "\t" << get_real_time() - time0 <<      endl;
        
        for (int j0=1; j0<phi.n2-1; j0++)
        {
            int run = 0;
            
            for ( ; run < n_runs; )
            {
                
                for (int k=0; k<n_series32; k++)
                   series32[k] = rng.next_uint();
                int32_to_bool_decomp( n_series32, series32, series );
                
                i[0] = i0;
                j[0] = j0;
                
                int k,c = 0;
                
                for ( k=0; k<walk_max; )
                {
                    if ( series[c++] )
                    {
                        if ( series[c++] )
                            i[k+1] = i[k]+1;
                        else
                            i[k+1] = i[k]-1;
                            
                            j[k+1] = j[k];
                    }
                    else
                    {
                        if ( series[c++] )
                            j[k+1] = j[k]+1;
                        else
                            j[k+1] = j[k]-1;
                            
                            i[k+1] = i[k];
                    }
                    
                    k++;
                    
                    if ( i[k] == 0 or i[k] == phi.n1-1 or j[k] == 0 or j[k] == phi.n2-1 )
                        break;
                }
               //     cout << k << " " ;
                if ( i[k] == 0 or i[k] == phi.n1-1 or j[k] == 0 or j[k] == phi.n2-1 )
                {
                    for ( int n=0; n<k; n++ )
                        count(i[k],j[k])(i[n],j[n])++;
                    run ++;
                }
            }
        }
    }
    
    
    
    for (int i=1; i<phi.n1-1; i++)
    for (int j=1; j<phi.n2-1; j++)
    {
        uint sum = 0;
        
        for (int k=0; k<phi.n1; k++)
        {
            sum += count(k,0)(i,j);
            sum += count(k,phi.n2-1)(i,j);
        }
        
        for (int k=0; k<phi.n2; k++)
        {
            sum += count(0,k)(i,j);
            sum += count(phi.n1-1,k)(i,j);
        }
        
        phi(i,j) = 0;
        
        for (int k=0; k<phi.n1; k++)
        {
            phi(i,j) += B(k,0)*(((double)count(k,0)(i,j))/sum);
            phi(i,j) += B(k,phi.n2-1)*(((double)count(k,phi.n2-1)(i,j))/sum);
        }
        
        for (int k=0; k<phi.n2; k++)
        {
            phi(i,j) += B(0,k)*(((double)count(0,k)(i,j))/sum);
            phi(i,j) += B(phi.n1-1,k)*(((double)count(phi.n1-1,k)(i,j))/sum);
        }
    }
    
    sprintf(fname,"%s/out/1/",pwd);
    mkdir(fname );
        
    for (int i=0; i<phi.n1; i++)
    {
        sprintf(fname,"%s/out/1/%06d_%06d.dat",pwd,i,0);
        writeFile(count(i,0),fname);
    }
    
    sprintf(fname,"%s/out/2/",pwd);
    mkdir(fname );
        
    for (int i=0; i<phi.n1; i++)
    {        
        sprintf(fname,"%s/out/2/%06d_%06d.dat",pwd,i,phi.n2-1);
        writeFile(count(i,phi.n2-1),fname);
    }
    
    sprintf(fname,"%s/out/3/",pwd);
    mkdir(fname );
    
    for (int i=0; i<phi.n2; i++)    
    {
        sprintf(fname,"%s/out/3/%06d_%06d.dat",pwd,0,i);
        writeFile(count(0,i),fname);
    }
    
    sprintf(fname,"%s/out/4/",pwd);
    mkdir(fname );
    
    for (int i=0; i<phi.n2; i++)    
    {
        sprintf(fname,"%s/out/4/%06d_%06d.dat",pwd,phi.n1-1,i);
        writeFile(count(phi.n1-1,i),fname);
    }
    
    for (int k=0; k<phi.n1*phi.n2; k++)
    {
        delete [] count[k].array;
        count[k].array = 0;
    }
    
    delete [] count.array;
    count.array = 0;
    
}



int main()
{
    std_setup();
    
    int n1=500 ,n2=n1;
    grid2D<double,double > grid( n1,-5,5, n2,-5,5 );
    grid = 0.0;
    
    grid2D<double,double > B(grid);
    B = 0.0;
    
    grid2D<double,double > phi(grid);
    phi = 0.0;
    
    for (int i=0; i<n1; i++)
    {
        B(i,0) = 2.0*cos(_2pi*B.x1(i)*0.5);
        B(i,n2-1) = 2.0*sin(_2pi*B.x1(i)*0.5)+15.0*cos(_2pi*B.x1(i)*0.1);
    }
    
    for (int i=0; i<n2; i++)
    {
        B(0,i) = 1.0*sin(_2pi*B.x2(i)*1.0);
        B(n1-1,i) = 1.0*cos(_2pi*B.x2(i)*1.0);        
    }
    
    solve_poisson( phi, B );
    
    grid2D<double,double > D2phi(phi);
    D2phi = 0.0;
   
    for (int i=1; i<phi.n1-1; i++)
    for (int j=1; j<phi.n2-1; j++)
        D2phi(i,j) = DEL2(phi,i,j);
        
    cout << L2norm(D2phi,2) << endl;

	sprintf(fname,"%s/out/Del2_phi1.png",pwd);
    plotGrid2D_1(D2phi,fname,cmap);
        
    Del2_FD(
	    phi,
	    D2phi,
	    10,
	    GRID_NONPER_BD_1,
	    10 );
	    
	cout << L2norm(D2phi,2) << endl;	
    
	sprintf(fname,"%s/out/phi.png",pwd);
    plotGrid2D_1(phi,fname,cmap);
    
	sprintf(fname,"%s/out/Del2_phi2.png",pwd);
    plotGrid2D_1(D2phi,fname,cmap);
    
    sprintf(fname,"%s/out/phi.dat",pwd);
    writeFile(phi,fname);
    
    double sum = 0;
    for (int i=0; i<n1; i++)
    for (int j=0; j<n2; j++)
        sum += D2phi(i,j);
    cout << sum << endl;
}









