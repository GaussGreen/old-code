/*******************************************************************************
**
**              Copyright (c) 1993 PARIBAS Capital Markets Group
**
********************************************************************************
**
**      SYSTEM:         SRT     SORT, Fixed Income 2020 Addins
**      SUB_SYSTEM:     OPT     Option Tools
**
**      MODULE NAME:    OPTION_TOOLS
**
**      PURPOSE:        C functions for option pricing
**
**      INCLUDED BY:    prog.c              !!!
**                      prog_fct.c          !!!
**                      gen_2020.c          !!!
**
**      AUTHORS:        Guillaume AMBLARD,
**                      Remy KLAMMERS
**                      Eric AULD
**                      Alexandra POITE
**
**      DATE:           1st October, 1992
**
**      VERSION:        1.1
**
**      DESCRIPTION:    Header for the library of general related option
**                      functions used in the SORT program applications
**
*******************************************************************************/

/*******************************************************************************
**                      Amendment History
*******************************************************************************/

/*******************************************************************************
**
**      AMENDED BY:     Ken N LINTON
**
**      DATE:           3rd May, 1994
**
**      VERSION:        <not applicable>
**
**      REASON:         Restructuring for re-use
**
**      REQUEST NO:     <not applicable>
**
**      DESCRIPTION:    <not applicable>
**
********************************************************************************
**
**      AMENDED BY:     G AMBLARD
**                      E AULD
**                      A POITE
**                      R KLAMMERS
**
**      DATE:           8th October, 1993
**
**      VERSION:        1.0
**
**      REASON:         Development and release
**
**      REQUEST NO:     <not applicable>
**
**      DESCRIPTION:    <not applicable>
**
*******************************************************************************/

/*******************************************************************************
**                      Include Files
*******************************************************************************/

#include "math.h"

#include "utallhdr.h"

/******************************************************************************/

#define SWAP(a, b) \
    tempr = (a);   \
    (a)   = (b);   \
    (b)   = tempr

/*******************************************************************************
 *
 * FUNCTION    	: comb(...)
 *
 * PURPOSE      	: Returns the binomial coefficient:
 *			/n_steps\ * p^(n-ii)(1-p)^ii
 *                       \   n   /
 *
 * DESCRIPTION  	: (Was originally 'static').
 *
 * CALLS		: <none>
 *
 * PARAMETERS   	: n     	- ??
 *              	: n_steps     	- ??
 *              	: p            	- ??
 *              	: ii          	- ??
 *
 * RETURNS      	: prob        	- ??
 *
 *******************************************************************************/

double comb(double n, double n_steps, double p, int ii)
{
    int    i;
    double prob;

    prob = 1.0;

    for (i = 1; i <= n; i++)
    {
        prob *= (n_steps + 1 - i) / i;
        prob *= (i <= ii ? 1 - p : p);
    }

    if (n >= ii)
        prob *= pow(p, n_steps - n);
    else
        prob *= pow(p, n_steps - ii) * pow(1 - p, ii - n);

    if (n == 0)
        prob = pow(p, n_steps - ii) * pow((1 - p), ii);
    if (n == -1)
        prob = 0;

    return (prob);

} /* END comb(...) */

/******************************************************************************/

/* *****************************************************************************
        Private Function Prototypes
   ************************************************************************** */

void realft(double data[], int n, int isign);

void four1(double data[], int nn, int isign);

/*******************************************************************************
 *
 * FUNCTION     	: smooth( 3 )
 *
 * PURPOSE      	: Smoothing module used for:
 *			increasing_strike_opt3(), and
 *			iko_ext3()
 *
 * DESCRIPTION  	: (Was originally 'static')
 *
 * CALLS		: <none>
 *
 * PARAMETERS   	: y[]       	- ??
 *              	: n         	- ??
 *              	: pts   	- number of ??
 *
 * RETURNS      : ??        	- ??
 *
 *******************************************************************************/

double smooth(double y[], int n, double pts)
{
    int nmin;
    int m = 2;
    int mo2;
    int k;
    int j;

    double yn;
    double y1;
    double rn1;
    double fac;
    double cnst;
    double real;

    nmin = n + (int)((2.0 * pts) + 0.5);

    while (m < nmin)
        m *= 2;
    cnst = pts / m, cnst = cnst * cnst;
    y1  = y[1];
    yn  = y[n];
    rn1 = 1.0 / (n - 1);

    for (j = 1; j <= n; j++)
        y[j] += (-rn1 * ((y1 * (n - j)) + (yn * (j - 1))));

    for (j = n + 1; j <= m; j++)
        y[j] = 0.0;

    mo2 = m >> 1;
    realft(y, mo2, 1);
    y[1] /= mo2;
    fac = 1.0;

    for (j = 1; j < mo2; j++)
    {
        k = (2 * j) + 1;

        if (fac)
        {
            if ((fac = (1.0 - cnst * j * j) / mo2) < 0.0)
                fac = 0.0;

            y[k]     = fac * y[k];
            y[k + 1] = fac * y[k + 1];
        }
        else
            y[k + 1] = y[k] = 0.0;
    }

    if ((fac = (1.0 - 0.25 * pts * pts) / mo2) < 0.0)
        fac = 0.0;

    y[2] *= fac;
    realft(y, mo2, -1);

    for (j = 1; j <= n; j++)
        y[j] += rn1 * ((y1 * (n - j)) + (yn * (j - 1)));

    real = y[6];

    return (real);

} /* END smooth(...) */

/******************************************************************************/

/*******************************************************************************
 *
 * FUNCTION     	: realft( 3 )
 *
 * PURPOSE      	: Real FFT function
 *
 * DESCRIPTION  	: ??
 *
 * CALLS		: four1( 3 )
 *
 * PARAMETERS   	: data[]	- ??
 *              	: n         	- ??
 *              	: isign     	- ??
 *
 * RETURNS	: void
 *
 *******************************************************************************/

void realft(double data[], int n, int isign)
{
    int i;
    int i1;
    int i2;
    int i3;
    int i4;
    int n2p3;

    double c1 = 0.5;
    double c2;
    double h1r;
    double h1i;
    double h2r;
    double h2i;

    double wr;
    double wi;
    double wpr;
    double wpi;
    double wtemp;
    double theta;

    theta = SRT_PI / (double)n;

    if (isign == 1)
    {
        c2 = -0.5;
        four1(data, n, 1);
    }
    else
    {
        c2    = 0.5;
        theta = -theta;
    }

    wtemp = sin(0.5 * theta);
    wpr   = -2.0 * wtemp * wtemp;
    wpi   = sin(theta);
    wr    = 1.0 + wpr;
    wi    = wpi;
    n2p3  = 2 * n + 3;

    for (i = 2; i <= n / 2; i++)
    {
        i4       = 1 + (i3 = n2p3 - (i2 = 1 + (i1 = i + i - 1)));
        h1r      = c1 * (data[i1] + data[i3]);
        h1i      = c1 * (data[i2] - data[i4]);
        h2r      = -c2 * (data[i2] + data[i4]);
        h2i      = c2 * (data[i1] - data[i3]);
        data[i1] = h1r + wr * h2r - wi * h2i;
        data[i2] = h1i + wr * h2i + wi * h2r;
        data[i3] = h1r - wr * h2r + wi * h2i;
        data[i4] = -h1i + wr * h2i + wi * h2r;
        wr       = (wtemp = wr) * wpr - wi * wpi + wr;
        wi       = wi * wpr + wtemp * wpi + wi;
    }

    if (isign == 1)
    {
        data[1] = (h1r = data[1]) + data[2];
        data[2] = h1r - data[2];
    }
    else
    {
        data[1] = c1 * ((h1r = data[1]) + data[2]);
        data[2] = c1 * (h1r - data[2]);
        four1(data, n, -1);
    }
} /* END realft(...) */

/******************************************************************************/

/*******************************************************************************
 *
 * FUNCTION     	: four1( 3 )
 *
 * PURPOSE      	: FFT Butterfly ????
 *
 * DESCRIPTION  	: ??
 *
 * CALLS		: SWAP( 2 )
 *
 * PARAMETERS   	: data[]    	- ??
 *              	: nn        	- ??
 *              	: isign     	- ??
 *
 * RETURNS      	: <void>
 *
 *******************************************************************************/

void four1(double data[], int nn, int isign)
{
    int n;
    int mmax;
    int m;
    int j;
    int istep;
    int i;

    double wtemp;
    double wr;
    double wpr;
    double wpi;
    double wi;
    double theta;
    double tempr;
    double tempi;

    n = nn << 1;
    j = 1;

    for (i = 1; i < n; i += 2)
    {
        if (j > i)
        {
            SWAP(data[j], data[i]);
            SWAP(data[j + 1], data[i + 1]);
        }
        m = n >> 1;
        while ((m >= 2) && (j > m))
        {
            j -= m;
            m >>= 1;
        }
        j += m;
    }
    mmax = 2;

    while (n > mmax)
    {
        istep = 2 * mmax;
        theta = 6.28318530717959 / (isign * mmax);
        wtemp = sin(0.5 * theta);
        wpr   = -2.0 * wtemp * wtemp;
        wpi   = sin(theta);
        wr    = 1.0;
        wi    = 0.0;

        for (m = 1; m < mmax; m += 2)
        {
            for (i = m; i <= n; i += istep)
            {
                j           = i + mmax;
                tempr       = wr * data[j] - wi * data[j + 1];
                tempi       = wr * data[j + 1] + wi * data[j];
                data[j]     = data[i] - tempr;
                data[j + 1] = data[i + 1] - tempi;
                data[i] += tempr;
                data[i + 1] += tempi;
            }
            wr = (wtemp = wr) * wpr - wi * wpi + wr;
            wi = wi * wpr + wtemp * wpi + wi;
        }
        mmax = istep;
    }

} /* END 	four1( 3 ) */

/******************************************************************************/

/******************************************************************************/

/*
 computate abscissas and weights for for n point Gauss-Hermite quadrature
 (mostly taken from Numerical Recipes); the reason you would use this is if
 you are integrating something of the form
 f(x) * gauss(x) from -infinity to infinity (something we do a lot);
 you call this function then the answer is sum(f(x[i])*w[i])

 This function is modified from the weight function exp(-x^2);
 this means multiplying the x values by sqrt(2) and the w values by sqrt(1/2*pi)

 E.Auld Jan 96
*/

/*#define GAUHER_EPS 3.0e-14  too small */
#define GAUHER_EPS 1.0e-8
#define GAUHER_PIM4 .7511255444649425
#define GAUHER_MAXIT 10

char* gauher(double* x, double* w, int n)
{
    int    i, its, j, m;
    double p1, p2, p3, pp, z, z1;

    m = (n + 1) / 2;
    for (i = 1; i <= m; i++)
    {
        if (i == 1)
        {
            z = sqrt((double)(2 * n + 1)) - 1.85575 * pow((double)(2 * n + 1), -0.16667);
        }
        else if (i == 2)
        {
            z -= 1.14 * pow((double)n, 0.426) / z;
        }
        else if (i == 3)
        {
            z = 1.86 * z - 0.86 * x[1];
        }
        else if (i == 4)
        {
            z = 1.91 * z - 0.91 * x[2];
        }
        else
        {
            z = 2.0 * z - x[i - 2];
        }
        for (its = 1; its <= GAUHER_MAXIT; its++)
        {
            p1 = GAUHER_PIM4;
            p2 = 0.0;
            for (j = 1; j <= n; j++)
            {
                p3 = p2;
                p2 = p1;
                p1 = z * sqrt(2.0 / j) * p2 - sqrt((double)(j - 1) / j) * p3;
            }
            pp = sqrt((double)2 * n) * p2;
            z1 = z;
            z  = z1 - p1 / pp;
            if (fabs(z - z1) <= GAUHER_EPS)
                break;
        }
        if (its > GAUHER_MAXIT)
            return serror("too many iterations in gauher");
        x[i]         = z;
        x[n + 1 - i] = -z;
        w[i]         = 2.0 / (pp * pp) * INV_SQRT_TWO_PI;
        w[n + 1 - i] = w[i];
    }
    z = sqrt(2.0);
    for (i = 1; i <= n; i++)
        x[i] *= z;
    return NULL;
}
