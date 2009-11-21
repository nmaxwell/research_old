
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
    
    // blue,cyan,green,yellow,magenta,red
    int levels[] = {20,40,60,80,200,4000};
    
    if (x <= levels[0])
    {
        float s = ((float)x/(2*levels[0]));
        return color3(0,0,s+0.5f);
    }    
    x -= levels[0];
    
    if (x <= levels[1])
    {
        float s = ((float)x/(2*levels[1]));
        return color3(0,s+0.5f,s+0.5f);
    }    
    x -= levels[1];
    
    if (x <= levels[2])
    {
        float s = ((float)x/(2*levels[2]));
        return color3(0,s+0.5f,0);
    }    
    x -= levels[2];
    
    if (x <= levels[3])
    {
        float s = ((float)x/(2*levels[1]));
        return color3(s+0.5f,s+0.5f,0);
    }    
    x -= levels[3];
    
    if (x <= levels[4])
    {
        float s = ((float)x/(2*levels[4]));
        return color3(s+0.5f,0,s+0.5f);
    }    
    x -= levels[4];
    
    if (x <= levels[5])
    {
        float s = ((float)x/(2*levels[5]));
        return color3(s+0.5f,0,0);
    }
    x -= levels[5];

    return c3_red;
}

color3 cmap(double x)
{
	float s = atan(x)/pi+0.5;
	return color3(s,s,s);
}



int main()
{
    std_setup();
    double time0 = get_real_time();    
    
    grid2D< > phi,D2phi;
    
    sprintf(fname,"%s/out/phi.dat",pwd);
    readFile(phi,fname);
    
    Del2_FD(
	    phi,
	    D2phi,
	    24,
	    GRID_NONPER_BD_1,
	    24 );
    
    double sum = 0;
    for (int i=0; i<phi.n1; i++)
    for (int j=0; j<phi.n2; j++)
        //cout << D2phi(i,j) << "\n";
        //cout << endl;
        sum += fabs(D2phi(i,j));
    cout << sum/(phi.n1*phi.n2) << endl;
    
     sum = 0;
    for (int i=0; i<phi.n1; i++)
    for (int j=0; j<phi.n2; j++)
        //cout << D2phi(i,j) << "\n";
        //cout << endl;
        sum += fabs(phi(i,j));
    cout << sum/(phi.n1*phi.n2) << endl;
    
    cout << endl;
    
    grid2D<int, int > temp;
    
    
    sprintf(fname,"%s/out/analysis/1/",pwd);
    mkdir(fname );
    sprintf(fname,"%s/out/analysis/2/",pwd);
    mkdir(fname );
    sprintf(fname,"%s/out/analysis/3/",pwd);
    mkdir(fname );
    sprintf(fname,"%s/out/analysis/4/",pwd);
    mkdir(fname );
    
    
    
    
    for (int i=0; i<phi.n1; i++)
    {
        sprintf(fname,"%s/out/1/%06d_%06d.dat",pwd,i,0);
        readFile(temp,fname);
        sprintf(fname,"%s/out/analysis/1/%06d_%06d.png",pwd,i,0);
        plotGrid2D_1(temp,fname,cmap);
        cout << i << "\t" << get_real_time() - time0 <<      endl;
    }
    
    for (int i=0; i<phi.n1; i++)
    {
        sprintf(fname,"%s/out/2/%06d_%06d.dat",pwd,i,phi.n2-1);
        readFile(temp,fname);
        sprintf(fname,"%s/out/analysis/2/%06d_%06d.png",pwd,i,phi.n2-1);
        plotGrid2D_1(temp,fname,cmap);
        cout << i << "\t" << get_real_time() - time0 <<      endl;
    }
   
    for (int i=0; i<phi.n2; i++)    
    {
        sprintf(fname,"%s/out/3/%06d_%06d.dat",pwd,0,i);
        readFile(temp,fname);
        sprintf(fname,"%s/out/analysis/3/%06d_%06d.png",pwd,0,i);
        plotGrid2D_1(temp,fname,cmap);
        cout << i << "\t" << get_real_time() - time0 <<      endl;
    }
    
    for (int i=0; i<phi.n2; i++)    
    {        
        sprintf(fname,"%s/out/4/%06d_%06d.dat",pwd,phi.n1-1,i);
        readFile(temp,fname);
        sprintf(fname,"%s/out/analysis/4/%06d_%06d.png",pwd,phi.n1-1,i);
        plotGrid2D_1(temp,fname,cmap);
        cout << i << "\t" << get_real_time() - time0 <<      endl;
    }
}









