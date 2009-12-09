
#define rgrid2D grid2D<double, double ,double >
#define cgrid2D grid2D<complex<double>, double ,double >




double l2_error( double *& exact, double *& approx, int n, double a, double b)
{
    double sum1 = 0.0;	    
    for (int k=0; k<n; k++)
	sum1 += (approx[k]-exact[k])*(approx[k]-exact[k]);
    
    double sum2 = 0.0;
    for (int k=0; k<n; k++)
	sum2 += exact[k]*exact[k];
    
    if ( sum2 > 0 )
        return sqrt(sum1/sum2);
    else
        return sqrt(sum1);
}


double l2_norm( double *& data, int n, double a, double b)
{
    double sum = 0.0;
    for (int k=0; k<n; k++)
        sum += data[k]*data[k];
    
    return sqrt(sum/n);
}


