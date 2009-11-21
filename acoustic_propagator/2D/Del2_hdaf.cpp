

#include <mathlib/math/hdaf/apply_hdaf.h>
#include <mathlib/math/differentiation/differentiation.h>



template<class St1,class Xt1,class St2,class Xt2 >
void Del2_hdaf(
	grid2D<double,St1,Xt1 > & in,
	grid2D<double,St2,Xt2 > & out,
	double k_cutoff,
	int m_x,
    int m_y
 )
{
    out.copy(in);
    
    double sigma_x = sqrt(2*m_x+1)/k_cutoff;
    double sigma_y = sqrt(2*m_y+1)/k_cutoff;
    
    double * temp = ml_alloc<double> ( max(in.n1, in.n2 ) );
    
    double eps_min = 1E-9;
    double eps_max = 1E-10;
    
    out = 0.0;
    
    for (int i=0; i<in.n1; i++)
    {
        apply_hdaf_reg(
            m_y, sigma_y, 2, in.dx2(), eps_min, eps_max,
            &( in(i,0) ),
            temp, in.n2,
            (int)(&(in(i,1))-&(in(i,0)) )  , 1 );            
         
       /*   differentiate(
	&( in(i,0) ),
    temp,
    in.n2,
    in.dx2(),
    2,
    (int)(&(in(i,1))-&(in(i,0)) ) , 1 )  ;*/
        
        for (int j=0; j<in.n2; j++)
            out(i,j) += temp[j];
    }
    
    for (int j=0; j<in.n2; j++)
    {        
        apply_hdaf_reg(
            m_x, sigma_x, 2, in.dx1(), eps_min, eps_max,
            &( in(0,j) ),
            temp, in.n1,
            (int)(&( in(1,j) ) - &( in(0,j) )), 1);
           
  /*    differentiate(
            &( in(0,j) ),
            temp,
            in.n1,
            in.dx1(),
            2,
            (int)(&( in(1,j) ) - &( in(0,j) )), 1 )  ;*/
        
       for (int i=0; i<in.n1; i++)
           out(i,j) += temp[i];
    }
    
    ml_free( temp );
}





















