#include "edginc/coreConfig.hpp"
#include <cmath>
#include <cstdio>

CORE_BEGIN_NAMESPACE

using namespace std;

double SC_InvertCumNormal( double x)

{
    static double  a0=   2.50662823884;
    static double  a1= -18.61500062529;
    static double  a2=  41.39119773534;
    static double  a3= -25.44106049637;

    static double  b0=   1.0;
    static double  b1=  -8.47351093090;
    static double  b2=  23.08336743743;
    static double  b3= -21.06224101826;
    static double  b4=   3.13082909833;

    static double  c[] = {  7.7108870705487895,
                            2.7772013533685169,
                            0.3614964129261002,
                            0.0373418233434554,
                            0.0028297143036967,
                            0.0001625716917922,
                            0.0000080173304740,
                            0.0000003840919865,
                            0.0000000129707170 };

    static double  k1=  0.4179886424926431;
    static double  k2=  4.2454686881376569;

    double d=0.,dd=0.,sv,y,y2,z,z2;

    int m=9;
    int j;

    y=x-0.5;

    if ( fabs(y) <= 0.42 ) {
        y2 = y*y;

        return( y*(a0+y2*(a1+y2*(a2+y2*a3)))/(b0+y2*(b1+y2*(b2+y2*(b3+y2*b4)))) );
    } else {
        y=0.5+fabs(y);
        z = k1*(2.*log(-log(1.-y))-k2);

        z2=2.*z;

        for (j=m-1; j>=1; j--) {
            sv=d;
            d=z2*d-dd+c[j];
            dd=sv;
        }
        y=z*d-dd+0.5*c[0];
        if (x>0.5)
            return(y);
        return(-y);
    }
}




CORE_END_NAMESPACE
