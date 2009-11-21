
#include <mathlib/math/random/ml_random.h>
#include <mathlib/tools/graphing/plot2D.h>

#include <mathlib/link.cpp>

#include <boost/math/special_functions/erf.hpp>

class gaussian : functor<double, double >
{
public:
    double mu,sig;
    
    gaussian(double mu, double sig ):mu(mu),sig(sig) {};
    
    double operator() (double const & x) const
    {
        return exp(-(x-mu)*(x-mu)/(2.0*sig*sig))/(sqrt2pi*sig);
    }
};

int main()
{
    std_setup();
    
    if (0)
    {
        plot2D plot(-4.0,4.0,-4.0,4.0, 600,600 );
        plot = c3_creme;
        plot.axes_1(0,0, 1, c3_creme,c3_black,plot2D_stdFormat,c3_black, 10,10,1,10  );
        
        
       // plot.plot( plot_schauder ,c3_blue);
        
        sprintf(fname,"%s/out/out.png",pwd);
        plot.png(fname);    
    }
    
    if(0)
    {
        int n_steps = 1000;
        float stop_time = 1.0;
        float time_step = stop_time/n_steps;
        float x_step = 0.1;
        
        for (int trial =0; trial<100; trial++)
        {
            double W[n_steps];
            double t[n_steps];
            
            for (int step = 0; step<n_steps; step++)
            {
                t[step] = time_step*step;
                W[step] = 0;
            }
            
            ml_random rng;
            
            int n_series = (n_steps/32+1);
            bool series[n_series*32];
            
            int series32[n_series];
            for (int k=0; k<n_series; k++)
                series32[k] = rng.gen_int();
            
            int32_to_bool_decomp( n_series, series32, series );
            
            W[0] = 0.0;
            
            for (int step = 1; step<n_steps; step++)
                if ( series[step] )
                    W[step] = W[step-1] + x_step;
                else
                    W[step] = W[step-1] - x_step;
            
                   
            plot2D plot(-0.1,stop_time*1.1,-10,10, 900,600 );
            plot = c3_creme;
            plot.axes_1(0,-7, 1, c3_creme,c3_black,plot2D_stdFormat,c3_black, 0.1,10,1,10  );
            
            plot.plot<double, double >(&(t[0]),&(W[0]),n_steps ,c3_black);
            
            sprintf(fname,"%s/out/%03d.png",pwd,trial);
            plot.png(fname);
        }
    }
    
    
    
    if(1)
    {
        uint n_trials = 5e7;
        int n_steps = 200;
        int sample_mod = 10;
        
        int n_samples = (n_steps-n_steps%sample_mod)/sample_mod;
        
        int range = 200;
        int origin = range/2;
        
        uint32 bins[n_samples][range];
        
        for (int s=0; s<n_samples; s++)
        for (int x=0; x<range; x++)
            bins[s][x] = 0;
        
        int n_series = (n_steps/32+1);
        bool series[n_series*32];
        int series32[n_series];
        ml_random rng(n_series*n_trials*1.2);
        
        for (uint trial=0; trial<n_trials; trial++)
        {
            //cout << trial << "\t" << ((float)trial)/n_trials <<  endl;
            
            for (int k=0; k<n_series; k++)
                series32[k] = rng.next_uint();
            
            int32_to_bool_decomp( n_series, series32, series );
            
            int x = origin;
            
            for (int step = 0; step<n_steps; step++)
            {
                if ( series[step] )
                    x++;
                else
                    x--;
                
                if (  x >= 0 && x < range ) // !(step%sample_mod) &&
                    bins[step/sample_mod][x] ++;
            }
        }
        
        cout << endl << endl;
        
        double density[n_samples][range];

        for (int s=0; s<n_samples; s++)
        {
            uint32 sum = 0;
            for (int x = 0; x<range; x++)
                sum += bins[s][x];
            
            for (int x = 0; x<range; x++)
                density[s][x] = ((double)bins[s][x])/sum;
        }
        
        double expectation[n_samples];
        double variance[n_samples];
        
        {
            for (int s=0; s<n_samples; s++)
            {
                expectation[s] = 0;
                
                for (int x = 0; x<range; x++)
                    expectation[s] += density[s][x]*x;
            }
            
            for (int s=0; s<n_samples; s++)
            {
                variance[s] = 0;
                
                for (int x = 0; x<range; x++)
                    variance[s] += density[s][x]*x*x;
            }
            
            for (int s = 0; s<n_samples; s++)
                variance[s] -= expectation[s]*expectation[s];
            
          /*  double mean = 0;
            for (int s=1; s<n_samples-1; s++)
            {
                mean += (variance[s+1]-variance[s])/(n_samples-2);
                cout << expectation[s] << "\t" << variance[s+1]-variance[s] << endl;
            }
            
            cout << mean << endl;*/
        }
        
        
        for (int s=0; s<n_samples; s++)
        {
            plot2D plot( 0, range , -0.1,.21, 900,600 );
            plot = c3_blue;
            //plot.axes_1(0,0, 1, c3_blue,c3_creme,plot2D_stdFormat,c3_creme, 10.0,1,0.2,10  );
            
            cout << expectation[s] << "\t" << variance[s] << endl;
            
            gaussian G(expectation[s], sqrt(variance[s]) );
            
            plot.plot(G, c3_white );
            
            for (int x = 0; x<range; x++)
                plot.ptDot( x , density[s][x], 2, c3_white );
            
            sprintf(fname,"%s/out/%03d.png",pwd,s*sample_mod);
          //  cout << fname << endl;
            plot.png(fname);
        }
        
        
    }
    
    
}




