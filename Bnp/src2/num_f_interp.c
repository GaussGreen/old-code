/* ===============================================================
   FILENAME : num_f_interp.c

   PURPOSE:   interpolation functions
   =============================================================== */

#include "utallhdr.h"
#include "math.h"
#include "num_h_interp.h"

/***
        interpolate x, y = f(x) at xt, linearly, or linear on xt*f(xt)
        if method==1, geometric if method == 2, geometric on rt if method==3
        on r*r*t if method=5(by Julia);
        x is assumed to be in increasing order
***/

double interp(double* x, double* y, int l, double xt, int method, double* q)
{
    int    i;
    double yt;

    /***	find times on either side of xt ***/

    i = 0;
    while (x[i] < xt)
    {
        i++;
        if (i >= l)
            break;
    }

    /**	boundary conditions **/

    if (i == 0)
    {
        if (q)
            *q = y[0];
        return y[0];
    }
    if (i == l)
    {
        if (q)
            *q = y[l - 1];
        return y[l - 1];
    }

    /**	linear interpol **/
    if (method == 0)
    {
        yt = y[i - 1] + (xt - x[i - 1]) * (y[i] - y[i - 1]) / (x[i] - x[i - 1]);
    }
    /** 	interpol on rate times time **/
    else if (method == 1)
    {
        yt = y[i - 1] * x[i - 1] +
             (xt - x[i - 1]) * (y[i] * x[i] - y[i - 1] * x[i - 1]) / (x[i] - x[i - 1]);
        yt /= xt;
    }

    /**	geometric interpol i.e. linear interpol on the log of y **/
    else if (method == 2)
    {
        yt = log(y[i - 1]) + (xt - x[i - 1]) * (log(y[i]) - log(y[i - 1])) / (x[i] - x[i - 1]);
        yt = exp(yt);
    }
    /**	linear interpol on (log(y)/x) **/
    else if (method == 3)
    {
        yt = log(y[i - 1]) / x[i - 1] +
             (xt - x[i - 1]) * (log(y[i]) / x[i] - log(y[i - 1]) / x[i - 1]) / (x[i] - x[i - 1]);
        yt = exp(yt * xt);
    }
    /** 	quadratic interpol  on rate*rate times time **/
    else if (method == 5)
    {
        /*
    if(i>2)
    yt=interp_quadratic(xt,x[i-2],y[i-2],x[i-1],y[i-1],x[i],y[i]);
    else yt = y[i-1] + (xt - x[i-1]) * (y[i] - y[i-1])/(x[i]-x[i-1]);
*/

        yt = y[i - 1] * y[i - 1] * x[i - 1] +
             (xt - x[i - 1]) / (x[i] - x[i - 1]) *
                 (y[i] * y[i] * x[i] - y[i - 1] * y[i - 1] * x[i - 1]);
        yt /= xt;
        yt = sqrt(yt);

        /*
                yt = sqrt(y[i-1])*x[i-1]+(xt-x[i-1])/(x[i]-x[i-1])
                     *(sqrt(y[i])*x[i]-sqrt(y[i-1])*x[i-1]);
                yt /= xt;
                yt=yt*yt;
                */
    }

    if (q)
        *q = yt;

    return yt;
}

/////////////////////////////////////////////////////////////////////////////////////////////////
////////////  same function as before in which one can interpolate the columns of matrices
/////////////////////////////////////////////////////////////////////////////////////////////////

int BinarySearch(double* xx, int n, double x)
{
    int ju, jm, jl;

    if (x <= xx[0])
    {
        return 0;
    }
    else if (x >= xx[n - 1])
    {
        return n;
    }

    jl = 0;
    ju = n + 1;

    //	ascnd=(xx[n] > xx[1]);

    while (ju - jl > 1)
    {
        jm = (int)(ju + jl) / 2;
        if (x >= xx[jm - 1])
            jl = jm;
        else
            ju = jm;
    }

    return jl;
}

double interp_columns(double** x, double** y, int l, double xt, int method, double* q, const int n)
{
    int    i;
    double yt;

    /***	find times on either side of xt ***/

    /*
            i = 0;
            while ( x[n][i] < xt )
            {
                    i++;
                    if ( i >= l )
                            break;
            }
    */

    i = BinarySearch(x[n], l, xt);

    /**	boundary conditions **/

    if (i == 0)
    {
        if (q)
            *q = y[n][0];
        return y[n][0];
    }
    if (i == l)
    {
        if (q)
            *q = y[n][l - 1];
        return y[n][l - 1];
    }

    /**	linear interpol **/
    if (method == 0)
    {
        yt = y[n][i - 1] + (xt - x[n][i - 1]) * (y[n][i] - y[n][i - 1]) / (x[n][i] - x[n][i - 1]);
    }
    /** 	interpol on rate times time **/
    else if (method == 1)
    {
        yt = y[n][i - 1] * x[n][i - 1] + (xt - x[n][i - 1]) *
                                             (y[n][i] * x[n][i] - y[n][i - 1] * x[n][i - 1]) /
                                             (x[n][i] - x[n][i - 1]);
        yt /= xt;
    }

    /**	geometric interpol i.e. linear interpol on the log of y **/
    else if (method == 2)
    {
        yt = log(y[n][i - 1]) +
             (xt - x[n][i - 1]) * (log(y[n][i]) - log(y[n][i - 1])) / (x[n][i] - x[n][i - 1]);
        yt = exp(yt);
    }
    /**	linear interpol on (log(y)/x) **/
    else if (method == 3)
    {
        yt = log(y[n][i - 1]) / x[n][i - 1] +
             (xt - x[n][i - 1]) * (log(y[n][i]) / x[n][i] - log(y[n][i - 1]) / x[n][i - 1]) /
                 (x[n][i] - x[n][i - 1]);
        yt = exp(yt * xt);
    }
    /** 	quadratic interpol  on rate*rate times time **/
    else if (method == 5)
    {
        /*
    if(i>2)
    yt=interp_quadratic(xt,x[i-2],y[i-2],x[i-1],y[i-1],x[i],y[i]);
    else yt = y[i-1] + (xt - x[i-1]) * (y[i] - y[i-1])/(x[i]-x[i-1]);
*/

        yt = y[n][i - 1] * y[n][i - 1] * x[n][i - 1] +
             (xt - x[n][i - 1]) / (x[n][i] - x[n][i - 1]) *
                 (y[n][i] * y[n][i] * x[n][i] - y[n][i - 1] * y[n][i - 1] * x[n][i - 1]);
        yt /= xt;
        yt = sqrt(yt);

        /*
                yt = sqrt(y[i-1])*x[i-1]+(xt-x[i-1])/(x[i]-x[i-1])
                     *(sqrt(y[i])*x[i]-sqrt(y[i-1])*x[i-1]);
                yt /= xt;
                yt=yt*yt;
                */
    }

    if (q)
        *q = yt;

    return yt;
}

double lin_interp_2d(
    double x, double y, double* xa, double* ya, double** za, long num_xa, long num_ya)

{
    long   start1 = 0, end1 = 0, start2 = 0, end2 = 0;
    double a, int1, int2, int3;

    if (x <= xa[0])
    {
        start1 = 0;
        end1   = 0;
    }
    else if (x >= xa[num_xa - 1])
    {
        start1 = num_xa - 1;
        end1   = num_xa - 1;
    }
    else
    {
        while (x > xa[end1])
            end1++;
        start1 = end1 - 1;
    }

    if (y <= ya[0])
    {
        start2 = 0;
        end2   = 0;
    }
    else if (y >= ya[num_ya - 1])
    {
        start2 = num_ya - 1;
        end2   = num_ya - 1;
    }
    else
    {
        while (y > ya[end2])
            end2++;
        start2 = end2 - 1;
    }

    a = (xa[start1] != xa[end1]) ? (za[end1][start2] - za[start1][start2]) / (xa[end1] - xa[start1])
                                 : 0.0;
    int1 = za[start1][start2] + a * (x - xa[start1]);

    a    = (xa[start1] != xa[end1]) ? (za[end1][end2] - za[start1][end2]) / (xa[end1] - xa[start1])
                                    : 0.0;
    int2 = za[start1][end2] + a * (x - xa[start1]);

    a = (ya[start2] != ya[end2]) ? (int2 - int1) / (ya[end2] - ya[start2]) : 0.0;

    int3 = int1 + a * (y - ya[start2]);

    return int3;
}

void hunt(double xx[], unsigned long n, double x, unsigned long* jlo)
{
    unsigned long jm, jhi, inc;
    int           ascnd;

    ascnd = (xx[n] >= xx[1]);
    if (*jlo <= 0 || *jlo > n)
    {
        *jlo = 0;
        jhi  = n + 1;
    }
    else
    {
        inc = 1;
        if (x >= xx[*jlo] == ascnd)
        {
            if (*jlo == n)
                return;
            jhi = (*jlo) + 1;
            while (x >= xx[jhi] == ascnd)
            {
                *jlo = jhi;
                inc += inc;
                jhi = (*jlo) + inc;
                if (jhi > n)
                {
                    jhi = n + 1;
                    break;
                }
            }
        }
        else
        {
            if (*jlo == 1)
            {
                *jlo = 0;
                return;
            }
            jhi = (*jlo)--;
            while (x < xx[*jlo] == ascnd)
            {
                jhi = (*jlo);
                inc <<= 1;
                if (inc >= jhi)
                {
                    *jlo = 0;
                    break;
                }
                else
                    *jlo = jhi - inc;
            }
        }
    }
    while (jhi - (*jlo) != 1)
    {
        jm = (jhi + (*jlo)) >> 1;
        if (x >= xx[jm] == ascnd)
            *jlo = jm;
        else
            jhi = jm;
    }
    if (x == xx[n])
        *jlo = n - 1;
    if (x == xx[1])
        *jlo = 1;
}

/*----------------------------------------------------------------------------------*/
// Lagrange Polinomial interpolation

double InterpLagrange1D(double* Xdata, double* Ydata, unsigned int n, double x)

{
    double       y = 0.0, FactNum, FactDen;
    unsigned int i, j;

    ///////////////////////////////////////////////////////////////////////
    ///// step: Computes all the n terms in the series
    ///////////////////////////////////////////////////////////////////////
    for (i = 1; i <= n; i++)
    {
        FactNum = 1.0;
        FactDen = 1.0;
        for (j = 1; j <= n; j++)
        {
            if (j == i) {}
            else
            {
                FactNum *= (x - Xdata[j]);
                FactDen *= (Xdata[i] - Xdata[j]);
            }
        }
        y += FactNum / FactDen * Ydata[i];
    }

    return y;
}

double InterpLagrangeColumns1D(double** Xdata, double** Ydata, unsigned int n, double x, int m)

{
    double       y = 0.0, FactNum, FactDen;
    unsigned int i, j;

    ///////////////////////////////////////////////////////////////////////
    ///// step: Computes all the n terms in the series
    ///////////////////////////////////////////////////////////////////////
    for (i = 0; i < n; i++)
    {
        FactNum = 1.0;
        FactDen = 1.0;
        for (j = 0; j < n; j++)
        {
            if (j == i) {}
            else
            {
                FactNum *= (x - Xdata[m][j]);
                FactDen *= (Xdata[m][i] - Xdata[m][j]);
            }
        }
        y += FactNum / FactDen * Ydata[m][i];
    }

    return y;
}