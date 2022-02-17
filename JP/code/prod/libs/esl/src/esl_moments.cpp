/****************************************************************************/
/*                                                                          */
/*      First and Second Moments of Tree Variables.                         */
/*                                                                          */
/****************************************************************************/
/*                                                                          */
/*      MOMENTS.c                                                           */
/*                                                                          */
/****************************************************************************/

#include <stdlib.h>  
#include <math.h>
#include "esl_moments.h"
#include "esl_macros.h"

/****************************************************************************/
/*                                                                          */
/*    Eq is treated explicitly as a special case of FX.                     */
/*                                                                          */
/*    Currencies are designated as follows:                                 */
/*                                                                          */
/*      Num: the numeraire currency (payment currency).                     */
/*      Ext: a currency other than the numeraire currency.                  */
/*      Dom: the currency of denomination for an equity or FX.              */
/*      For: the foreign currency whose price "FX" represents.              */
/*                                                                          */
/*    Further:                                                              */
/*                                                                          */
/*      Dom0: the domestic currency for fx/equity asset 0.                  */
/*      Dom1: the domestic currency for fx/equity asset 1.                  */
/*      For0: the foreign  currency for fx/equity asset 0 (must be FX).     */
/*      For1: the foreign  currency for fx/equity asset 1 (must be FX).     */
/*                                                                          */
/*    Other processes:                                                      */
/*                                                                          */
/*      Fx:    an FX asset.                                                 */
/*      Eq:    an equity asset.                                             */
/*      QxExt: a background FX process that relates Num to Ext.             */
/*      QxDom: a background FX process that relates Num to Dom              */
/*             for an equity or FX.                                         */
/*      QxFor: a background FX process that relates Num to For              */
/*             for an equity or FX.                                         */
/*                                                                          */
/*    Asset labels are usually translated during calculation to fit         */
/*    into the following order, as in the document.                         */
/*                                                                          */
/*      0:  foreign currency                                                */
/*      1:  domestic currency                                               */
/*      2:  fx                                                              */
/*      3:  eq                                                              */
/*      4:  other currency                                                  */ 
/*                                                                          */
/*      00: QxFor                                                           */
/*      11: QxDom                                                           */
/*      44: QxExt                                                           */
/*                                                                          */
/*    Rij represents correlation between i and j.                           */
/*    RiQ represents correlation between i and ii (in quanto terms).        */
/*                                                                          */
/****************************************************************************/

/*****  xi  *****************************************************************/
/*                                                                          */
/*      Pure math function.                                                 */
/*                                                                          */
/****************************************************************************/

double xi (double b,
           double t)
{
    if (fabs(b) > TINY)
    {
        return (1 - exp(-b * t)) / b;
    }
    else
    {
        return t;
    }
}

/*****  nu  *****************************************************************/
/*                                                                          */
/*      Pure math function.                                                 */
/*                                                                          */
/****************************************************************************/

double nu (double b,
           double t)
{
    if (fabs(b) > TINY)
    {
        return (t - xi(b, t)) / b;
    }
    else
    {
        return 0.5 * t * t;
    }
}

/*****  eta  ****************************************************************/
/*                                                                          */
/*      Pure math function.                                                 */
/*                                                                          */
/****************************************************************************/

double eta (double b0,
            double b1,
            double t)
{
    if (fabs(b1) > TINY)
    {
        return (xi(b0, t) - xi(b0 + b1, t)) / b1;
    }
    else
    {
        return exp(-b0 * t) * nu(-b0, t);
    }
}

/*****  zeta  ***************************************************************/
/*                                                                          */
/*      Pure math function.                                                 */
/*                                                                          */
/****************************************************************************/

double zeta (double b0,
             double b1,
             double t)
{
    if (fabs(b0 + b1) > TINY)
    {
        return (xi(b0, t) - exp(-b0 * t) * xi(b1, t)) / (b0 + b1);
    }
    else
    {
        return nu(b0, t);
    }
}

/*****  tau  ****************************************************************/
/*                                                                          */
/*      Pure math function.                                                 */
/*                                                                          */
/****************************************************************************/

double tau (double b0,
            double b1,
            double t)
{
    if (fabs(b0 * b1) > TINY)
    {
        return (t - xi(b0, t) - xi(b1, t) + xi(b0 + b1, t)) / (b0 * b1);
    }
    else if ((fabs(b0) >  TINY) &&
             (fabs(b1) <= TINY))
    {
        return (0.5 * t * t - exp(-b0 * t) * nu(-b0, t)) / b0;
    }
    else if ((fabs(b0) <= TINY) &&
             (fabs(b1) >  TINY))
    {
        return (0.5 * t * t - exp(-b1 * t) * nu(-b1, t)) / b1;
    }
    else
    {
        return 0.5 * t * t * t;
    }
}

/*****  TreeMoments_H  ******************************************************/
/*                                                                          */
/*      See the document "Moments in Four Dimensions".                      */
/*                                                                          */
/****************************************************************************/

int TreeMoments_H (double *H,              /* (O) */
                   int     NbTP,           /* (I) */
                   double *Time,           /* (I) */
                   double  MrIr0,          /* (I) */
                   double *SvIr0,          /* (I) */
                   double *FwIr1,          /* (I) */
                   double  MrIr1,          /* (I) */
                   double *SvIr1,          /* (I) */
                   double *CrIr0Ir1)       /* (I) */
{
    int    k   = 0;
    double t   = 0;
    double tn  = 0;
    double dt  = 0;

    double r0  = 0;
    double b0  = MrIr0;
    double s0  = 0;
    double xi0 = 0;
    double z0  = 0;

    double r1  = 0;
    double b1  = MrIr1;
    double s1  = 0;
    double xi1 = 0;
    double z1  = 0;
    
    double R01 = 0;
    
    double ha  = 0;
    double hb  = 0;
    double hc  = 0;
    double h   = 0;

    for (k=0; k<NbTP; ++k)
    {
        H[k] = h;

        t   = Time[k];
        tn  = Time[k+1];
        dt  = tn - t;
        s0  = SvIr0[k];
        r1  = FwIr1[k];
        s1  = SvIr1[k];
        R01 = CrIr0Ir1[k];

        /* locals */

        xi0 = xi(b0, dt);
        xi1 = xi(b1, dt);

        z0 = r0 * exp(-b0 * t) * xi0;
        z1 = r1 * exp(-b1 * t) * xi1;

        /* updates */

        hb += z1 * ha;
        ha += R01 * s0 * s1 * exp((b0 + b1) * tn) * xi(b0 + b1, dt);
        hc += R01 * s0 * s1 * r1 * exp(b0 * tn) * zeta(b0, b1, dt);

        h = exp(-b0 * tn) * (hb + hc);
    }

    H[NbTP] = h;

    return SUCCESS;
}

/*****  TreeMoments_II  *****************************************************/
/*                                                                          */
/*      See the document "Moments in Four Dimensions".                      */
/*                                                                          */
/****************************************************************************/

int TreeMoments_II (double *II,             /* (O) */
                    int     NbTP,           /* (I) */
                    double *Time,           /* (I) */
                    double  MrIr0,          /* (I) */
                    double *SvIr0,          /* (I) */
                    double  MrIr1,          /* (I) */
                    double *SvIr1,          /* (I) */
                    double *CrIr0Ir1)       /* (I) */
{
    int    k   = 0;
    double t   = 0;
    double tn  = 0;
    double dt  = 0;

    double b0  = MrIr0;
    double s0  = 0;

    double b1  = MrIr1;
    double s1  = 0;

    double R01 = 0;

    double value = 0;

    for (k=0; k<NbTP; ++k)
    {
        II[k] = value;

        t   = Time[k];
        tn  = Time[k+1];
        dt  = tn - t;

        s0  = SvIr0[k];
        s1  = SvIr1[k];

        R01 = CrIr0Ir1[k];

        value *= exp(-(b0 + b1) * dt);
        value += R01 * s0 * s1 * xi(b0 + b1, dt);
    }

    II[NbTP] = value;

    return SUCCESS;
}
    
/*****  TreeMoments_MM  *****************************************************/
/*                                                                          */
/*      See the document "Moments in Four Dimensions".                      */
/*                                                                          */
/****************************************************************************/

int TreeMoments_MM (double *MM,             /* (O) */
                    int     NbTP,           /* (I) */
                    double *Time,           /* (I) */
                    double *FwIr0,          /* (I) */
                    double  MrIr0,          /* (I) */
                    double *SvIr0,          /* (I) */
                    double *FwIr1,          /* (I) */
                    double  MrIr1,          /* (I) */
                    double *SvIr1,          /* (I) */
                    double *CrIr0Ir1)       /* (I) */
{
    int    k   = 0;
    double t   = 0;
    double tn  = 0;
    double dt  = 0;

    double r0  = 0;
    double b0  = MrIr0;
    double s0  = 0;
    double xi0 = 0;
    double en0 = 0;
    double z0  = 0;

    double r1  = 0;
    double b1  = MrIr1;
    double s1  = 0;
    double xi1 = 0;
    double en1 = 0;
    double z1  = 0;

    double R01  = 0;
    double en01 = 0;
    double xi01 = 0;

    double mma = 0;
    double mmb = 0;
    double mmc = 0;
    double mmd = 0;
    double mme = 0;
    double mmf = 0;
    double mmg = 0;
    double mmh = 0;

    double c = 0;

    double value = 0;

    for (k=0; k<NbTP; ++k)
    {
        MM[k] = value;

        t   = Time[k];
        tn  = Time[k+1];
        dt  = tn - t;

        r0  = FwIr0[k];
        s0  = SvIr0[k];

        r1  = FwIr1[k];
        s1  = SvIr1[k];

        R01 = CrIr0Ir1[k];

        /* locals */

        xi0  = xi(b0, dt);
        xi1  = xi(b1, dt);

        xi01 = xi(b0 + b1, dt);

        en0  = exp(b0 * tn);
        en1  = exp(b1 * tn);
        en01 = exp((b0 + b1) * tn);

        z0   = r0 * exp(-b0 * t) * xi0;
        z1   = r1 * exp(-b1 * t) * xi1;

        /* updates */

        mma += R01 * s0 * s1 * r0 * r1 * tau(b0, b1, dt);

        mmc += z0 * mmb;

        mmb += R01 * s0 * s1 * r1 * en0 * eta(b0, b1, dt);

        mme += z1 * mmd;

        mmd += R01 * s1 * s0 * r0 * en1 * eta(b1, b0, dt);

        mmh += z1 * mmf + z0 * mmg + z0 * z1 * c;

        mmf += z0 * c;

        mmg += z1 * c;

        c   += R01 * s0 * s1 * en01 * xi01;

        value = mma + mmc + mme + mmh;
    }

    MM[NbTP] = value;

    return SUCCESS;
}
    
/*****  TreeMean_IrExt  *****************************************************/
/*                                                                          */
/*      See the document "Moments in Four Dimensions".                      */
/*                                                                          */
/****************************************************************************/

int TreeMean_IrExt (double *Mean,         /* (O) */
                    int     NbTP,         /* (I) */
                    double *Time,         /* (I) */
                    double  MrIrExt,      /* (I) */
                    double *SvIrExt,      /* (I) */
                    double *SvQxExt,      /* (I) */
                    double *CrIrExtQxExt) /* (I) */
{
    int    k   = 0;
    double t   = 0;
    double dt  = 0;

    double b4  = MrIrExt;
    double s4  = 0;
    double xi4 = 0;

    double s44 = 0;

    double R4Q = 0;

    double m   = 0;

    for (k=0; k<NbTP; ++k)
    {
        Mean[k] = m;

        t   = Time[k];
        dt  = Time[k+1] - Time[k];
        s4  = SvIrExt[k];
        s44 = SvQxExt[k];
        R4Q = CrIrExtQxExt[k];

        xi4 = xi(b4, dt);

        m *= exp(-b4 * dt);
        m -= R4Q * s4 * s44 * xi4;
    }

    Mean[NbTP] = m;

    return SUCCESS;
}

/*****  TreeMean_IrNum  *****************************************************/
/*                                                                          */
/*      See the document "Moments in Four Dimensions".                      */
/*                                                                          */
/****************************************************************************/

int TreeMean_IrNum (double *Mean,  /* (O) */
                    int     NbTP)  /* (I) */
{
    int k;

    for (k=0; k<=NbTP; ++k)
    {
        Mean[k] = 0;
    }

    Mean[NbTP] = 0;

    return SUCCESS;
}

/*****  TreeMean_FxExt  *****************************************************/
/*                                                                          */
/*      See the document "Moments in Four Dimensions".                      */
/*                                                                          */
/****************************************************************************/

int TreeMean_FxExt (double *Mean,           /* (O) */
                    int     NbTP,           /* (I) */
                    double *Time,           /* (I) */
                    double *FwIrDom,        /* (I) */
                    double  MrIrDom,        /* (I) */
                    double *SvIrDom,        /* (I) */
                    double *FwIrFor,        /* (I) */
                    double  MrIrFor,        /* (I) */
                    double *SvIrFor,        /* (I) */
                    double *SvQxDom,        /* (I) */
                    double *SvQxFor,        /* (I) */
                    double *CrIrDomQxDom,   /* (I) */
                    double *CrIrForQxFor)   /* (I) */
{
    int    k   = 0;
    double t   = 0;
    double tn  = 0;
    double dt  = 0;

    double r0  = 0;
    double b0  = MrIrFor;
    double s0  = 0;
    double xi0 = 0;
    double z0  = 0;
    double H0  [MAXNBDATE];
    double h0  = 0;
    double n0  = 0;
    double q0a = 0;
    double q0b = 0;
    double q0c = 0;
    double q0  = 0;

    double r1  = 0;
    double b1  = MrIrDom;
    double s1  = 0;
    double xi1 = 0;
    double z1  = 0;
    double H1  [MAXNBDATE];
    double h1  = 0;
    double n1  = 0;
    double q1a = 0;
    double q1b = 0;
    double q1c = 0;
    double q1  = 0;

    double s00 = 0;
    double s11 = 0;

    double R0Q = 0;
    double R1Q = 0;

    double One [MAXNBDATE];

    double m = 0;

    for (k=0; k<MAXNBDATE; ++k)
    {
        One[k] = 1;
    }

    TreeMoments_H(H0,
                  NbTP,
                  Time,
                  MrIrFor,
                  SvIrFor,
                  FwIrFor,
                  MrIrFor,
                  SvIrFor,
                  One);

    TreeMoments_H(H1,
                  NbTP,
                  Time,
                  MrIrDom,
                  SvIrDom,
                  FwIrDom,
                  MrIrDom,
                  SvIrDom,
                  One);

    for (k=0; k<NbTP; ++k)
    {
        Mean[k] = m;

        t   = Time[k];
        tn  = Time[k+1];
        dt  = tn - t;

        r0  = FwIrFor[k];
        s0  = SvIrFor[k];
        h0  = H0[k];

        r1  = FwIrDom[k];
        s1  = SvIrDom[k];
        h1  = H1[k];

        s00 = SvQxFor[k];
        s11 = SvQxDom[k];

        R0Q = CrIrForQxFor[k];
        R1Q = CrIrDomQxDom[k];

        /* locals */

        xi0 = xi(b0, dt);
        xi1 = xi(b1, dt);

        z0 = r0 * exp(-b0 * t) * xi0;
        z1 = r1 * exp(-b1 * t) * xi1;

        /* q0 updates */

        q0b += z0 * q0a;
        q0a += R0Q * s0 * s00 * exp(b0 * tn) * xi0;
        q0c += r0 * R0Q * s0 * s00 * nu(b0, dt);

        q0 = q0b + q0c;

        /* q1 updates */

        q1b += z1 * q1a;
        q1a += R1Q * s1 * s11 * exp(b1 * tn) * xi1;
        q1c += r1 * R1Q * s1 * s11 * nu(b1, dt);

        q1 = q1b + q1c;

        /* n0 updates */

        n0 += r0 * h0 * dt;

        /* n1 updates */

        n1 += r1 * h1 * dt;

        /* m */
        
        m = (n1 - q1) -
            (n0 - q0);
    }

    Mean[NbTP] = m;

    return SUCCESS;
}

/*****  TreeMean_FxNum  *****************************************************/
/*                                                                          */
/*      See the document "Moments in Four Dimensions".                      */
/*                                                                          */
/****************************************************************************/

int TreeMean_FxNum (double *Mean,           /* (O) */
                    int     NbTP,           /* (I) */
                    double *Time,           /* (I) */
                    double *FwIrDom,        /* (I) */
                    double  MrIrDom,        /* (I) */
                    double *SvIrDom,        /* (I) */
                    double *FwIrFor,        /* (I) */
                    double  MrIrFor,        /* (I) */
                    double *SvIrFor,        /* (I) */
                    double *SvFx,           /* (I) */
                    double *CrIrForFx)      /* (I) */
{
    double  Zero[MAXNBDATE];

    double *SvQxDom = Zero;
    double *SvQxFor = SvFx;

    double *CrIrDomQxDom = Zero;
    double *CrIrForQxFor = CrIrForFx;

    int k;

    for (k=0; k<MAXNBDATE; ++k)
    {
        Zero[k] = 0;
    }

    TreeMean_FxExt (Mean,
                    NbTP,
                    Time,
                    FwIrDom,
                    MrIrDom,
                    SvIrDom,
                    FwIrFor,
                    MrIrFor,
                    SvIrFor,
                    SvQxDom,
                    SvQxFor,
                    CrIrDomQxDom,
                    CrIrForQxFor);

    return SUCCESS;
}

/*****  TreeMean_EqExt  *****************************************************/
/*                                                                          */
/*      See the document "Moments in Four Dimensions".                      */
/*                                                                          */
/****************************************************************************/

int TreeMean_EqExt (double *Mean,           /* (O) */
                    int     NbTP,           /* (I) */
                    double *Time,           /* (I) */
                    double *FwIrDom,        /* (I) */
                    double  MrIrDom,        /* (I) */
                    double *SvIrDom,        /* (I) */
                    double *SvQxDom,        /* (I) */
                    double *CrIrDomQxDom)   /* (I) */
{
    double  Zero[MAXNBDATE];

    double *FwIrFor = Zero;
    double  MrIrFor = 0;
    double *SvIrFor = Zero;

    double *SvQxFor = Zero;

    double *CrIrForQxFor = Zero;

    int k;

    for (k=0; k<MAXNBDATE; ++k)
    {
        Zero[k] = 0;
    }

    TreeMean_FxExt (Mean,
                    NbTP,
                    Time,
                    FwIrDom,
                    MrIrDom,
                    SvIrDom,
                    FwIrFor,
                    MrIrFor,
                    SvIrFor,
                    SvQxDom,
                    SvQxFor,
                    CrIrDomQxDom,
                    CrIrForQxFor);

    return SUCCESS;
}

/*****  TreeMean_EqNum  *****************************************************/
/*                                                                          */
/*      See the document "Moments in Four Dimensions".                      */
/*                                                                          */
/****************************************************************************/

int TreeMean_EqNum (double *Mean,           /* (O) */
                    int     NbTP,           /* (I) */
                    double *Time,           /* (I) */
                    double *FwIrDom,        /* (I) */
                    double  MrIrDom,        /* (I) */
                    double *SvIrDom)        /* (I) */
{
    double  Zero[MAXNBDATE];

    double *SvQxDom = Zero;

    double *CrIrDomQxDom = Zero;

    int k;

    for (k=0; k<MAXNBDATE; ++k)
    {
        Zero[k] = 0;
    }

    TreeMean_EqExt (Mean,
                    NbTP,
                    Time,
                    FwIrDom,
                    MrIrDom,
                    SvIrDom,
                    SvQxDom,
                    CrIrDomQxDom);

    return SUCCESS;
}

/*****  TreeCovariance_IrIr  ************************************************/
/*                                                                          */
/*      See the document "Moments in Four Dimensions".                      */
/*                                                                          */
/****************************************************************************/

int TreeCovariance_IrIr (double *Covariance,     /* (O) */
                         int     NbTP,           /* (I) */
                         double *Time,           /* (I) */
                         double  MrIr0,          /* (I) */
                         double *SvIr0,          /* (I) */
                         double  MrIr1,          /* (I) */
                         double *SvIr1,          /* (I) */
                         double *CrIr0Ir1)       /* (I) */
{
    TreeMoments_II(Covariance,
                   NbTP,
                   Time,
                   MrIr0,
                   SvIr0,
                   MrIr1,
                   SvIr1,
                   CrIr0Ir1);
    
    return SUCCESS;
}

/*****  TreeCovariance_IrFx  ************************************************/
/*                                                                          */
/*      See the document "Moments in Four Dimensions".                      */
/*                                                                          */
/****************************************************************************/

int TreeCovariance_IrFx (double *Covariance,     /* (O) */
                         int     NbTP,           /* (I) */
                         double *Time,           /* (I) */
                         double  MrIrExt,        /* (I) */
                         double *SvIrExt,        /* (I) */
                         double *FwIrFor,        /* (I) */
                         double  MrIrFor,        /* (I) */
                         double *SvIrFor,        /* (I) */
                         double *FwIrDom,        /* (I) */
                         double  MrIrDom,        /* (I) */
                         double *SvIrDom,        /* (I) */
                         double *SvFx,           /* (I) */
                         double *CrIrExtIrFor,   /* (I) */
                         double *CrIrExtIrDom,   /* (I) */
                         double *CrIrExtFx)      /* (I) */
{
    int k;

    double H40  [MAXNBDATE];
    double H41  [MAXNBDATE];

    double II42 [MAXNBDATE];

    TreeMoments_H(H40,
                  NbTP,
                  Time,
                  MrIrExt,
                  SvIrExt,
                  FwIrFor,
                  MrIrFor,
                  SvIrFor,
                  CrIrExtIrFor);

    TreeMoments_H(H41,
                  NbTP,
                  Time,
                  MrIrExt,
                  SvIrExt,
                  FwIrDom,
                  MrIrDom,
                  SvIrDom,
                  CrIrExtIrDom);

    TreeMoments_II(II42,
                   NbTP,
                   Time,
                   MrIrExt,
                   SvIrExt,
                   0,
                   SvFx,
                   CrIrExtFx);

    for (k=0; k<=NbTP; ++k)
    {
        Covariance[k] = H41[k] - H40[k] + II42[k];
    }
    
    return SUCCESS;
}

/*****  TreeCovariance_IrEq  ************************************************/
/*                                                                          */
/*      See the document "Moments in Four Dimensions".                      */
/*                                                                          */
/****************************************************************************/

int TreeCovariance_IrEq (double *Covariance,     /* (O) */
                         int     NbTP,           /* (I) */
                         double *Time,           /* (I) */
                         double  MrIrExt,        /* (I) */
                         double *SvIrExt,        /* (I) */
                         double *FwIrDom,        /* (I) */
                         double  MrIrDom,        /* (I) */
                         double *SvIrDom,        /* (I) */
                         double *SvEq,           /* (I) */
                         double *CrIrExtIrDom,   /* (I) */
                         double *CrIrExtEq)      /* (I) */
{
    double  Zero[MAXNBDATE];

    double *FwIrFor = Zero;
    double  MrIrFor = 0; 
    double *SvIrFor = Zero;

    double *SvFx    = SvEq;

    double *CrIrExtIrFor = Zero;
    double *CrIrExtFx    = CrIrExtEq;

    int k;

    for (k=0; k<MAXNBDATE; ++k)
    {
        Zero[k] = 0;
    }

    TreeCovariance_IrFx (Covariance,
                         NbTP,
                         Time,
                         MrIrExt,
                         SvIrExt,
                         FwIrFor,
                         MrIrFor,
                         SvIrFor,
                         FwIrDom,
                         MrIrDom,
                         SvIrDom,
                         SvFx,   
                         CrIrExtIrFor,
                         CrIrExtIrDom,
                         CrIrExtFx);

    return SUCCESS;
}

/*****  TreeCovariance_FxFx  ************************************************/
/*                                                                          */
/*      See the document "Moments in Four Dimensions".                      */
/*                                                                          */
/****************************************************************************/

int TreeCovariance_FxFx (double *Covariance,     /* (O) */
                         int     NbTP,           /* (I) */
                         double *Time,           /* (I) */
                         double *FwIrFor0,       /* (I) */
                         double  MrIrFor0,       /* (I) */
                         double *SvIrFor0,       /* (I) */
                         double *FwIrDom0,       /* (I) */
                         double  MrIrDom0,       /* (I) */
                         double *SvIrDom0,       /* (I) */
                         double *FwIrFor1,       /* (I) */
                         double  MrIrFor1,       /* (I) */
                         double *SvIrFor1,       /* (I) */
                         double *FwIrDom1,       /* (I) */
                         double  MrIrDom1,       /* (I) */
                         double *SvIrDom1,       /* (I) */
                         double *SvFx0,          /* (I) */
                         double *SvFx1,          /* (I) */
                         double *CrIrFor0IrFor1, /* (I) */
                         double *CrIrFor0IrDom1, /* (I) */
                         double *CrIrFor0Fx1,    /* (I) */
                         double *CrIrDom0IrFor1, /* (I) */
                         double *CrIrDom0IrDom1, /* (I) */
                         double *CrIrDom0Fx1,    /* (I) */
                         double *CrIrFor1Fx0,    /* (I) */
                         double *CrIrDom1Fx0,    /* (I) */
                         double *CrFx0Fx1)       /* (I) */
{
    int k;

    /* S0 = Fx0    */
    /* D0 = IrDom0 */
    /* F0 = IrFor0 */

    /* S1 = Fx1    */
    /* D1 = IrDom1 */
    /* F1 = IrFor1 */

    double HS0D1  [MAXNBDATE];
    double HS0F1  [MAXNBDATE];
    double HS1D0  [MAXNBDATE];
    double HS1F0  [MAXNBDATE];

    double IIS0S1 [MAXNBDATE];

    double MMD0D1 [MAXNBDATE];
    double MMD0F1 [MAXNBDATE];
    double MMF0D1 [MAXNBDATE];
    double MMF0F1 [MAXNBDATE];

    TreeMoments_H(HS0D1,
                  NbTP,
                  Time,
                  0,
                  SvFx0,
                  FwIrDom1,
                  MrIrDom1,
                  SvIrDom1,
                  CrIrDom1Fx0);

    TreeMoments_H(HS0F1,
                  NbTP,
                  Time,
                  0,
                  SvFx0,
                  FwIrFor1,
                  MrIrFor1,
                  SvIrFor1,
                  CrIrFor1Fx0);

    TreeMoments_H(HS1D0,
                  NbTP,
                  Time,
                  0,
                  SvFx1,
                  FwIrDom0,
                  MrIrDom0,
                  SvIrDom0,
                  CrIrDom0Fx1);

    TreeMoments_H(HS1F0,
                  NbTP,
                  Time,
                  0,
                  SvFx1,
                  FwIrFor0,
                  MrIrFor0,
                  SvIrFor0,
                  CrIrFor0Fx1);

    TreeMoments_II(IIS0S1,
                   NbTP,
                   Time,
                   0,
                   SvFx0,
                   0,
                   SvFx1,
                   CrFx0Fx1);

    TreeMoments_MM(MMD0D1,
                   NbTP,
                   Time,
                   FwIrDom0,
                   MrIrDom0,
                   SvIrDom0,
                   FwIrDom1,
                   MrIrDom1,
                   SvIrDom1,
                   CrIrDom0IrDom1);

    TreeMoments_MM(MMD0F1,
                   NbTP,
                   Time,
                   FwIrDom0,
                   MrIrDom0,
                   SvIrDom0,
                   FwIrFor1,
                   MrIrFor1,
                   SvIrFor1,
                   CrIrDom0IrFor1);

    TreeMoments_MM(MMF0D1,
                   NbTP,
                   Time,
                   FwIrFor0,
                   MrIrFor0,
                   SvIrFor0,
                   FwIrDom1,
                   MrIrDom1,
                   SvIrDom1,
                   CrIrFor0IrDom1);

    TreeMoments_MM(MMF0F1,
                   NbTP,
                   Time,
                   FwIrFor0,
                   MrIrFor0,
                   SvIrFor0,
                   FwIrFor1,
                   MrIrFor1,
                   SvIrFor1,
                   CrIrFor0IrFor1);

    for (k=0; k<=NbTP; ++k)
    {
        Covariance[k] =   MMD0D1 [k]
                        - MMD0F1 [k]
                        - MMF0D1 [k]
                        + MMF0F1 [k]
                        + HS0D1  [k]
                        - HS0F1  [k]
                        + HS1D0  [k]
                        - HS1F0  [k]
                        + IIS0S1 [k];
    } 
    
    return SUCCESS;
}

/*****  TreeCovariance_FxEq  ************************************************/
/*                                                                          */
/*      See the document "Moments in Four Dimensions".                      */
/*                                                                          */
/****************************************************************************/

int TreeCovariance_FxEq (double *Covariance,     /* (O) */
                         int     NbTP,           /* (I) */
                         double *Time,           /* (I) */
                         double *FwIrFor0,       /* (I) */
                         double  MrIrFor0,       /* (I) */
                         double *SvIrFor0,       /* (I) */
                         double *FwIrDom0,       /* (I) */
                         double  MrIrDom0,       /* (I) */
                         double *SvIrDom0,       /* (I) */
                         double *FwIrDom1,       /* (I) */
                         double  MrIrDom1,       /* (I) */
                         double *SvIrDom1,       /* (I) */
                         double *SvFx,           /* (I) */
                         double *SvEq,           /* (I) */
                         double *CrIrFor0IrDom1, /* (I) */
                         double *CrIrFor0Eq,     /* (I) */
                         double *CrIrDom0IrDom1, /* (I) */
                         double *CrIrDom0Eq,     /* (I) */
                         double *CrIrDom1Fx,     /* (I) */
                         double *CrFxEq)         /* (I) */
{
    double  Zero[MAXNBDATE];

    double *FwIrFor1 = Zero;
    double  MrIrFor1 = 0;
    double *SvIrFor1 = Zero;

    double *SvFx0    = SvFx;
    double *SvFx1    = SvEq;

    double *CrIrFor0IrFor1 = Zero;
    double *CrIrFor0Fx1    = CrIrFor0Eq;
    double *CrIrDom0IrFor1 = Zero;
    double *CrIrDom0Fx1    = CrIrDom0Eq;
    double *CrIrFor1Fx0    = Zero;
    double *CrIrDom1Fx0    = CrIrDom1Fx;
    double *CrFx0Fx1       = CrFxEq;

    int k;

    for (k=0; k<MAXNBDATE; ++k)
    {
        Zero[k] = 0;
    }

    TreeCovariance_FxFx (Covariance,
                         NbTP,
                         Time,
                         FwIrFor0,
                         MrIrFor0,
                         SvIrFor0,
                         FwIrDom0,
                         MrIrDom0,
                         SvIrDom0,
                         FwIrFor1,
                         MrIrFor1,
                         SvIrFor1,
                         FwIrDom1,
                         MrIrDom1,
                         SvIrDom1,
                         SvFx0,
                         SvFx1,
                         CrIrFor0IrFor1,
                         CrIrFor0IrDom1,
                         CrIrFor0Fx1,
                         CrIrDom0IrFor1,
                         CrIrDom0IrDom1,
                         CrIrDom0Fx1,
                         CrIrFor1Fx0,
                         CrIrDom1Fx0,
                         CrFx0Fx1);

    return SUCCESS;
}

/*****  TreeCovariance_EqEq  ************************************************/
/*                                                                          */
/*      See the document "Moments in Four Dimensions".                      */
/*                                                                          */
/****************************************************************************/

int TreeCovariance_EqEq (double *Covariance,     /* (O) */
                         int     NbTP,           /* (I) */
                         double *Time,           /* (I) */
                         double *FwIrDom0,       /* (I) */
                         double  MrIrDom0,       /* (I) */
                         double *SvIrDom0,       /* (I) */
                         double *FwIrDom1,       /* (I) */
                         double  MrIrDom1,       /* (I) */
                         double *SvIrDom1,       /* (I) */
                         double *SvEq0,          /* (I) */
                         double *SvEq1,          /* (I) */
                         double *CrIrDom0IrDom1, /* (I) */
                         double *CrIrDom0Eq1,    /* (I) */
                         double *CrIrDom1Eq0,    /* (I) */
                         double *CrEq0Eq1)       /* (I) */
{
    double Zero[MAXNBDATE];

    double *FwIrFor0 = Zero;
    double  MrIrFor0 = 0;
    double *SvIrFor0 = Zero;

    double *FwIrFor1 = Zero;
    double  MrIrFor1 = 0;
    double *SvIrFor1 = Zero;

    double *SvFx0    = SvEq0;
    double *SvFx1    = SvEq1;

    double *CrIrFor0IrFor1 = Zero;
    double *CrIrFor0IrDom1 = Zero;
    double *CrIrFor0Fx1    = Zero;
    double *CrIrDom0IrFor1 = Zero;
    double *CrIrDom0Fx1    = CrIrDom0Eq1;
    double *CrIrFor1Fx0    = Zero;
    double *CrIrDom1Fx0    = CrIrDom1Eq0;
    double *CrFx0Fx1       = CrEq0Eq1;

    int k;

    for (k=0; k<MAXNBDATE; ++k)
    {
        Zero[k] = 0;
    }

    TreeCovariance_FxFx (Covariance,
                         NbTP,
                         Time,
                         FwIrFor0,
                         MrIrFor0,
                         SvIrFor0,
                         FwIrDom0,
                         MrIrDom0,
                         SvIrDom0,
                         FwIrFor1,
                         MrIrFor1,
                         SvIrFor1,
                         FwIrDom1,
                         MrIrDom1,
                         SvIrDom1,
                         SvFx0,   
                         SvFx1,   
                         CrIrFor0IrFor1,
                         CrIrFor0IrDom1,
                         CrIrFor0Fx1,   
                         CrIrDom0IrFor1,
                         CrIrDom0IrDom1,
                         CrIrDom0Fx1,   
                         CrIrFor1Fx0,   
                         CrIrDom1Fx0,   
                         CrFx0Fx1);
    
    return SUCCESS;
}
