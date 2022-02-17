
#include "edginc/coreConfig.hpp"

#include <cmath>

CORE_BEGIN_NAMESPACE

using namespace std;

/* Returns the cumulative normal density for a given x */

double
SC_Normal(double x)
{
    static double a1    =  0.319381530 ;
    static double a2    = -1.1164195538 ;
    static double a3    = -4.9962391853 ;
    static double a4    = -1.0223286152 ;
    static double a5    = -0.7304159999 ;
    static double Gamma =  0.2316419 ;
    static double oneoversqrttwopi =0.39894228040143300000;

    double k ;
    double derivn ;
    double absx ;
    double auxnormal ;

    derivn = oneoversqrttwopi * exp(-0.5*x*x) ;
    if (x >= 0)
        absx = x ;
    else
        absx = - x ;
    k = 1.0 / (1.0 + Gamma * absx) ;
    auxnormal = 1 - derivn * a1 * k * (1.0 + a2 * k * (1.0 + a3 * k * (1.0 + a4 * k * (1.0 + a5 * k)))) ;

    if (x >= 0)
        return auxnormal ;
    else
        return (1.0 - auxnormal) ;
}

CORE_END_NAMESPACE
