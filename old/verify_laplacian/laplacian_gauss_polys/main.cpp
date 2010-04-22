
#include <typeinfo>

#include <mathlib/math/std_math.h>

#include <mathlib/link.cpp>




double P0 [] = { 1.0 };

double P1 [] = { 4.0, -4.0 };

double P2 [] = { 16.0, -64.0, 32.0 };

double P3 [] = { 64.0, -576.0, 1152.0, -384.0 };

double P4 [] = { 256.0, -4096.0, 18432.0, -24576.0, 6144.0 };

double P5 [] = { 1024.0, -25600.0, 204800.0, -614400.0, 614400.0, -122880.0 };

double P6 [] = { 4096.0, -147456.0, 1843200.0, - 9830400.0,  + 22118400.0,  - 17694720.0,  + 2949120.0 };

double P7 [] = { 16384.0,  - 802816.0,  + 14450688.0,  - 120422400.0,  + 481689600.0,  - 867041280.0,  + 578027520.0,  - 82575360.0 };

double P8 [] = { 65536.0,  - 4194304.0,  + 102760448.0,  - 1233125376.0,  + 7707033600.0,  - 24662507520.0,  + 36993761280.0,  - 21139292160.0,  + 2642411520.0 };

double P9 [] = { 262144.0,  - 21233664.0,  + 679477248.0,  - 11098128384.0,  + 99883155456.0,  - 499415777280.0,  + 1331775406080.0,  - 1712282664960.0,  + 856141332480.0,  - 95126814720.0 };

double P10 [] = { 1048576.0, 0 - 104857600.0,  + 4246732800.0,  - 90596966400.0,  + 1109812838400.0,  - 7990652436480.0,  + 33294385152000.0,  - 76101451776000.0,  + 85614133248000.0,  - 38050725888000.0,  + 3805072588800.0 };

double P11 [] = { 4194304.0, 2 - 507510784.0, 0 + 25375539200.0,  - 685139558400.0,  + 10962232934400.0,  - 107429882757120.0,  + 644579296542720.0,  - 2302068916224000.0,  + 4604137832448000.0,  - 4604137832448000.0,  + 1841655132979200.0,  - 167423193907200.0 };

double P12 [] = { 16777216.0, 4 - 2415919104.0, 2 + 146163105792.0, 0 - 4872103526400.0,  + 98660096409600.0,  - 1262849234042880.0,  + 10313268744683520.0,  - 53039667829800960.0,  + 165748961968128000.0,  - 294664821276672000.0,  + 265198339149004800.0,  - 96435759690547200.0,  + 8036313307545600.0 };

double P13 [] = { 67108864.0, 6 - 11341398016.0, 4 + 816580657152.0, 2 - 32935419838464.0, 0 + 823385495961600.0,  - 13338845034577920.0,  + 142281013702164480.0,  - 995967095915151360.0,  + 4481851931618181120.0,  - 12449588698939392000.0,  + 19919341918303027200.0,  - 16297643387702476800.0,  + 5432547795900825600.0,  - 417888291992371200.0 };

double P14 [] = { 268435456.0, 8 - 52613349376.0, 6 + 4445828022272.0, 4 - 213399745069056.0, 2 + 6455342288338944.0, 0 - 129106845766778880.0,  + 1742942417851514880.0,  - 15935473534642421760.0,  + 97604775399684833280.0,  - 390419101598739333120.0,  + 976047753996848332800.0,  - 1419705823995415756800.0,  + 1064779367996561817600.0,  - 327624420922019020800.0,  + 23401744351572787200.0 };

double P15 [] = { 1073741824.0,  - 241591910400.0, 8 + 23676007219200.0, 6 - 1333748406681600.0, 4 + 48014942640537600.0, 2 - 1161961611901009920.0, 0 + 19366026865016832000.0,  - 224092596580909056000.0,  + 1792740772647272448000.0,  - 9760477539968483328000.0,  + 35137719143886539980800.0,  - 79858452599742136320000.0,  + 106477936799656181760000.0,  - 73715494707454279680000.0,  + 21061569916415508480000.0,  - 1404104661094367232000.0 };

double * P [] = { &(P0[0]), &(P1[0]), &(P2[0]), &(P3[0]), &(P4[0]), &(P5[0]), &(P6[0]), &(P7[0]), &(P8[0]), &(P9[0]), &(P10[0]), &(P11[0]), &(P12[0]), &(P13[0]), &(P14[0]), &(P15[0]) };



double gauss_2d_lap_poly_0 [] = { 1.000000000000000000e+00 };
double gauss_2d_lap_poly_1 [] = { -4.000000000000000000e+00, 4.000000000000000000e+00 };
double gauss_2d_lap_poly_2 [] = { 3.200000000000000000e+01, -6.400000000000000000e+01, 1.600000000000000000e+01 };
double gauss_2d_lap_poly_3 [] = { -3.840000000000000000e+02, 1.152000000000000000e+03, -5.760000000000000000e+02, 6.400000000000000000e+01 };
double gauss_2d_lap_poly_4 [] = { 6.144000000000000000e+03, -2.457600000000000000e+04, 1.843200000000000000e+04, -4.096000000000000000e+03, 2.560000000000000000e+02 };
double gauss_2d_lap_poly_5 [] = { -1.228800000000000000e+05, 6.144000000000000000e+05, -6.144000000000000000e+05, 2.048000000000000000e+05, -2.560000000000000000e+04, 1.024000000000000000e+03 };
double gauss_2d_lap_poly_6 [] = { 2.949120000000000000e+06, -1.769472000000000000e+07, 2.211840000000000000e+07, -9.830400000000000000e+06, 1.843200000000000000e+06, -1.474560000000000000e+05, 4.096000000000000000e+03 };
double gauss_2d_lap_poly_7 [] = { -8.257536000000000000e+07, 5.780275200000000000e+08, -8.670412800000000000e+08, 4.816896000000000000e+08, -1.204224000000000000e+08, 1.445068800000000000e+07, -8.028160000000000000e+05, 1.638400000000000000e+04 };
double gauss_2d_lap_poly_8 [] = { 2.642411520000000000e+09, -2.113929216000000000e+10, 3.699376128000000000e+10, -2.466250752000000000e+10, 7.707033600000000000e+09, -1.233125376000000000e+09, 1.027604480000000000e+08, -4.194304000000000000e+06, 6.553600000000000000e+04 };
double gauss_2d_lap_poly_9 [] = { -9.512681472000000000e+10, 8.561413324800000000e+11, -1.712282664960000000e+12, 1.331775406080000000e+12, -4.994157772800000000e+11, 9.988315545600000000e+10, -1.109812838400000000e+10, 6.794772480000000000e+08, -2.123366400000000000e+07, 2.621440000000000000e+05 };
double gauss_2d_lap_poly_10 [] = { 3.805072588800000000e+12, -3.805072588800000000e+13, 8.561413324800000000e+13, -7.610145177600000000e+13, 3.329438515200000000e+13, -7.990652436480000000e+12, 1.109812838400000000e+12, -9.059696640000000000e+10, 4.246732800000000000e+09, -1.048576000000000000e+08, 1.048576000000000000e+06 };
double gauss_2d_lap_poly_11 [] = { -1.674231939072000000e+14, 1.841655132979200000e+15, -4.604137832448000000e+15, 4.604137832448000000e+15, -2.302068916224000000e+15, 6.445792965427200000e+14, -1.074298827571200000e+14, 1.096223293440000000e+13, -6.851395584000000000e+11, 2.537553920000000000e+10, -5.075107820000000000e+08, 4.194304000000000000e+06 };
double gauss_2d_lap_poly_12 [] = { 8.036313307545600000e+15, -9.643575969054720000e+16, 2.651983391490048000e+17, -2.946648212766720000e+17, 1.657489619681280000e+17, -5.303966782980096000e+16, 1.031326874468352000e+16, -1.262849234042880000e+15, 9.866009640960000000e+13, -4.872103526400000000e+12, 1.461631057940000000e+11, -2.415919100000000000e+09, 1.677721600000000000e+07 };
double gauss_2d_lap_poly_13 [] = { -4.178882919923712000e+17, 5.432547795900825600e+18, -1.629764338770247680e+19, 1.991934191830302720e+19, -1.244958869893939200e+19, 4.481851931618181120e+18, -9.959670959151513600e+17, 1.422810137021644800e+17, -1.333884503457792000e+16, 8.233854959616000000e+14, -3.293541983846200000e+13, 8.165806571560000000e+11, -1.134139801000000000e+10, 6.710886400000000000e+07 };
double gauss_2d_lap_poly_14 [] = { 2.340174435157278720e+19, -3.276244209220190208e+20, 1.064779367996561818e+21, -1.419705823995415757e+21, 9.760477539968483328e+20, -3.904191015987393331e+20, 9.760477539968483328e+19, -1.593547353464242176e+19, 1.742942417851514880e+18, -1.291068457667788800e+17, 6.455342288338946000e+15, -2.133997450690520000e+14, 4.445828022278000000e+12, -5.261334936800000000e+10, 2.684354560000000000e+08 };
double gauss_2d_lap_poly_15 [] = { -1.404104661094367232e+21, 2.106156991641550848e+22, -7.371549470745427968e+22, 1.064779367996561818e+23, -7.985845259974213632e+22, 3.513771914388653998e+22, -9.760477539968483328e+21, 1.792740772647272448e+21, -2.240925965809090560e+20, 1.936602686501683200e+19, -1.161961611901009920e+18, 4.801494264053760000e+16, -1.333748406681594000e+15, 2.367600721920800000e+13, -2.415919104000000000e+11, 1.073741824000000000e+09 };

double * gauss_2d_lap_poly [] = { &(gauss_2d_lap_poly_0[0]), &(gauss_2d_lap_poly_1[0]), &(gauss_2d_lap_poly_2[0]), &(gauss_2d_lap_poly_3[0]), &(gauss_2d_lap_poly_4[0]), &(gauss_2d_lap_poly_5[0]), &(gauss_2d_lap_poly_6[0]), &(gauss_2d_lap_poly_7[0]), &(gauss_2d_lap_poly_8[0]), &(gauss_2d_lap_poly_9[0]), &(gauss_2d_lap_poly_10[0]), &(gauss_2d_lap_poly_11[0]), &(gauss_2d_lap_poly_12[0]), &(gauss_2d_lap_poly_13[0]), &(gauss_2d_lap_poly_14[0]), &(gauss_2d_lap_poly_15[0]) };





int main()
{
    for (int n=0; n<=15; n++ )
    {
        printf ( "double gauss_2d_lap_poly_%d [] = { ", n );
        for (int k=0; k<=n-1; k++ )
            printf( "%22.18e, 0.0, ", P[n][n - k] );
            printf( "%22.18e };\n", P[n][0] );
        
    }
    
    cout << endl;
    
    {
        int n = 15;
        printf ( "int gauss_2d_lap_poly_degrees [] = { " );
        for (int k=0; k<=n-1; k++ )
            printf ( "%d, ", 2*k );
            printf ( "%d };\n", 2*n );
        
    }
    
    cout << endl;
    
    {
        int n = 15;
        printf ( "double * gauss_2d_lap_poly [] = { " );
        for (int k=0; k<=n-1; k++ )
            printf ( "&(gauss_2d_lap_poly_%d[0]), ", k );
            printf ( "&(gauss_2d_lap_poly_%d[0]) };\n", n );
        
    }
}

