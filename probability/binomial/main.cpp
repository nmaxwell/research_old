
#include <mathlib/math/std_math.h>
#include <mathlib/math/random/ml_random.h>

#include <mathlib/link.cpp>
#include <mathlib/non-lib_things.h>

template <class T >
void output(T *& X, T *& Y1, T *& Y2, int n, const char * name )
{
    ofstream out;
    char fname[200];
    //sprintf(fname,"/workspace/output/out_%d.py", pwd, count);
    sprintf(fname,"%s/%s.py", pwd, name);
    out.open(fname, fstream::out);
    assert(out.good() && out.is_open());
    
    out << "\nimport pylab as p\n";
    out << "X = [";
    for (int k=0; k<n; k++)
        out << X[k] << ",";
    out << "];\n";
    out << "Y1 = [";
    for (int k=0; k<n; k++)
        out << Y1[k] << ",";
    out << "];\n";
    out << "Y2 = [";
    for (int k=0; k<n; k++)
        out << Y2[k] << ",";
    out << "];\n";
    out << "p.plot(X,Y1,X,Y2) \np.show()\n\n";
    
    out.close();
}

template <class T >
void output( T *& Y1, T *& Y2, int n, const char * name )
{
    ofstream out;
    char fname[200];
    //sprintf(fname,"/workspace/output/out_%d.py", pwd, count);
    out.open(name, fstream::out);
    assert(out.good() && out.is_open());
    
    out << "\nimport pylab as p\n";
    out << "Y1 = [";
    for (int k=0; k<n; k++)
        out << Y1[k] << ",";
    out << "];\n";
    out << "Y2 = [";
    for (int k=0; k<n; k++)
        out << Y2[k] << ",";
    out << "];\n";
    out << "p.plot(Y1) \np.plot(Y2) \np.show()\n\n";
    
    out.close();
}

inline int opp(int const & x)
{
    if (x == 0) return 1;
    else return 0;
}

int main()
{
    std_setup();
    
    int n = 1E6;
    int n2 = 1E2;
    
    ml_random R(4*n/32);
    
    bool * x = ml_alloc<bool > (n);
    
    int n_bins = 50;    
    uint32 * bins_0 = ml_alloc <uint32 > (n_bins);
    uint32 * bins_1 = ml_alloc <uint32 > (n_bins);
    
    for ( int k=0; k<n_bins; k++)
    {
        bins_0[k] = 0;
        bins_1[k] = 0;
    }
    
    int total_score = 0;
    uint32 guesses = 0;
    
    for (int t=0; t<n2; t++ )
    {
       // cout << t << endl;
        R.refresh_pool();
        int32_to_bool_decomp( n/32, (int * )(R.pool), x );
        
        int l = 1;
        bool s = x[0];
        bool prediction = 0;
        
        int score = 0;
        
        for ( int k=1; k<n; k++)
        {
            if (l >= 10)
            {
                prediction = !s;                
                
                if (x[k] == prediction) score ++;
                    else score --;
                    
                guesses++;
            }
            
            
            if ( x[k] == s ) l++;
            else
            {
                if (s==0) {
                    if ( l<n_bins-1 ) bins_0[l] ++;
                    else bins_0[n_bins-1] ++; }
                
                if (s==1) {
                    if ( l<n_bins-1 ) bins_1[l] ++;
                    else bins_1[n_bins-1] ++; }     
                
                l=1;
            }
            
            s = x[k];
        }
        
        cout << score << endl;
        
        total_score += score;
    }
    
    int64 count = n;
    count *= n2;
    
    cout << "\ntotal_score: " << total_score << "\t" << guesses << endl;
    cout << (double(total_score))/guesses << endl;
    
    

    

 
    if(0)
     {
             for ( int l=0; l<n_bins; l++)
        cout << l << "\t" << (bins_0[l]) << ", " << (bins_1[l]) << endl;
        
        int streak_count_0 = 0;
        int streak_count_1 = 0;
        
        for ( int l=0; l<n_bins; l++)
            streak_count_0 += (bins_0[l]);
        
        for ( int l=0; l<n_bins; l++)
            streak_count_1 += (bins_1[l]);
        
        cout  << streak_count_0 << "\t" << streak_count_1 << endl; 
     }
    
    output( bins_0, bins_1, n_bins, "/workspace/output/out.py" );
    
    std_exit();
}

