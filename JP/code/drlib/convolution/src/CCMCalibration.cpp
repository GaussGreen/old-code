//----------------------------------------------------------------------------
//
//   Group       : Derivatives Research
//
//   Filename    : CCMCalibration.cpp
//
//   Description : 
//
//   Date        : Aug 2004
//
//
//----------------------------------------------------------------------------
#include "edginc/config.hpp"
#include "edginc/CCMCalibration.hpp"
#include "edginc/Format.hpp"
#include "edginc/mathlib.hpp"

DRLIB_BEGIN_NAMESPACE

typedef double (*TObjectFunc)(double, void*);

static double ccmRootFindBrent(           
   TObjectFunc funcd,                   /* (I) Function to call */
   void       *data,                    /* (I) Data to pass into funcd */
   double      boundLo,                 /* (I) Lower bound on legal X */
   double      boundHi,                 /* (I) Upper bound on legal X */
   int         numIterations,           /* (I) Maximum number of iterations */
   double      guess,                   /* (I) Initial guess */
   double      initialXStep,            /* (I) Size of step in x */
   double      initialFDeriv,           /* (I) Initial derivative or 0*/
   double      xacc,                    /* (I) X accuracy tolerance */
   double      facc);                   /* (I) Function accuracy tolerance */

/* ---------------------------------------------------------------------------
 * Recovery Calibration
 */
/* loss given default payoff calculation utility */
static double lgdCalc(double x,double k1,double k2) 
{
    ASSERT(k2>=k1);
    if (x<=k1) return k1;
    if (x>=k2) return k2;
    return x;
}

/* struct with the parameters needed for recovery calibration */
struct SolverData
{
    double survival;
    double pdep;
    double pind;
    double pgauss;
    double l1;
    double l2;
    double lbar;
    double W1ind;
    double indep;
    double betaDef;
    double betaRec;
    double qM;
};

/* Objective function, find fK such that expected loss is matched */
static double recoveryFuncToSolve(double fK, void *data){
    const char routine[] = "recoveryFuncToSolve";
    const struct SolverData *d  = (const struct SolverData *)data;
    double fP = N1InverseBetter(d->pgauss);
    double K  = N1(fK);  
    double c  = d->betaDef * d->betaRec;

    if (d->qM != 0.) {
        throw ModelException(routine, 
                             "cannot have non 0 qm and correlated recovery");
    }

    double tmp = N2(fP, fK, c);
    double tmp1  = K - tmp * d->pind;
    tmp1 *= (1. - d->indep);
    tmp1 += d->indep * d->W1ind * (1. - d->pind*d->pgauss);   
    tmp1  = d->l2 * (1. - d->pind*d->pgauss) - (d->l2 - d->l1)*tmp1;
    
    return ((d->pdep - d->survival) / d->pdep * d->lbar - tmp1);
}

/** Caluclates the weight on LGD1 as a function of M and recovery
 * parameters */
double CCMCalibration::recoveryWeightCalc(
    double  W1ind,          /* (I) */
    double  T1,             /* (I) */
    bool    isT1Calibrated, /* (I) */
    double  a,              /* (I) */
    double  beta_r,         /* (I) */
    double  qM,             /* (I) */
    double  M){             /* (I) */
    if (isT1Calibrated) {
        if (qM != 0.) {
            throw ModelException("CCMCalibration::recoveryWeightCalc",
                                 "qM is not zero");
        }
        if (beta_r == 1.0) {
            if (T1 < M) {
                return a*W1ind;
            } else if (T1 > M) {
                return a*W1ind + (1-a);
            } else {
                throw ModelException(
                    "CCMCalibration::recoveryWeightCalc",
                    "Normal distribution input is undefined");
            }
        } else {
            return (a*W1ind + (1-a) * N1((T1-beta_r*M)/sqrt(1.-beta_r*beta_r)));
        }
    } else {
        return W1ind;
    }
}


#define TINY 1e-12

/* Calibration of recovery parameters for one name */
void CCMCalibration::recoveryModelCalibrate(
    const CCMConvolution::NameParam&  b, /* (I) info on 1 name: survival, loss, correlation */
    double                        pind,   //(I) indep survival
    double                        pgauss, //(I) gaussian copula survival
    bool                          isDisc, //(I) do we calibrate to discretised loss 0=false
    CCMConvolution::LgdParamSP&   p){     //(O) calibrated loss param 
    const char routine[] = "recoveryModelCalibrate";
    try{
        p->isT1Calibrated = false;
        p->T1 = 0.;

        double RC = b.R*b.cataRecFactor;
        double RNC;
        if (fabs(b.survival-b.pdep) > TINY){
            RNC = (b.R*(1.-b.survival)-RC*(1.-b.pdep))
                /((1.-b.survival)-(1.-b.pdep));    
        } else {
            RNC = RC;
        }

        p->lc = b.ntl * lgdCalc(b.lgdNotional- b.R*b.cataRecFactor, 
                                 b.lgdFloor, b.lgdCap);

        if (fabs(b.ntl) < TINY) {
            p->LGD1  = p->LGD2  = 0;
            p->LGD1d = p->LGD2d = 0;
            p->W1ind = 1.;
            return;
        }

        /* 1. Desired Recovery levels */
        double R1 = (1. - b.lambda) * RNC;
        double R2 = RNC + b.lambda*(1. - RNC);
    
        /* 2. Rounding */
        /* ensure abs(LGD1)>=abs(LGD2) */
        double l1 = b.ntl*lgdCalc(b.lgdNotional-R1,
                                  b.lgdFloor,
                                  b.lgdCap);
        double l2 = b.ntl*lgdCalc(b.lgdNotional-R2,
                                  b.lgdFloor,
                                  b.lgdCap);
        p->LGD1 = (long)(b.ntl>0 ? ceil(l1-TINY)  : -ceil(-l1-TINY));
        p->LGD2 = (long)(b.ntl>0 ? floor(l2+TINY) : -floor(-l2+TINY));
        p->LGD1d = l1;
        p->LGD2d = l2;

#ifdef ASSERT_CHECK
        if (fabs(l2-l1) < TINY && fabs(l1-p->LGD1) < TINY){
            ASSERT(p->LGD1 == p->LGD2);
        }
#endif

        if (p->LGD1 == p->LGD2) {
            p->W1ind = 1.;
            return;
        }

        double Lbar = lgdCalc(b.lgdNotional-RNC,b.lgdFloor,
                              b.lgdCap)*b.ntl;
    
        p->W1ind = (Lbar - p->LGD2) / (p->LGD1-p->LGD2);
        ASSERT(fabs(p->W1ind-.5) <= .5);
        if (pgauss > 1.-3e-16){
            return;
        }
        /* performance: do not calibrate if beta_r = 0 */
        if (b.betaRec==0.){
            return;
        }
        /* 3. Calibration of T1 */
        {
            struct SolverData data;
            data.survival = b.survival;
            data.pdep     = b.pdep;
            data.pind     = pind;
            data.pgauss   = pgauss;
            data.l1       = isDisc ? p->LGD1 : p->LGD1d;
            data.l2       = isDisc ? p->LGD2 : p->LGD2d;
            data.lbar     = Lbar;
            data.W1ind    = p->W1ind;
            data.indep    = b.indep;
            data.betaDef  = b.beta;
            data.betaRec  = b.betaRec;
            data.qM       = b.qM;
            p->T1 = ccmRootFindBrent(           
                recoveryFuncToSolve,    /* (I) Function to call              */
                &data,                   /* (I) Data to pass into funcd      */
                -10.,                   /* (I) Lower bound on legal X        */
                10.,                    /* (I) Upper bound on legal X        */
                100,                    /* (I) Maximum number of iterations  */
                0.,                     /* (I) Initial guess                 */
                0.01,                   /* (I) Size of step in x             */
                0.,                     /* (I) Initial derivative or 0       */
                1e-7,                   /* (I) X accuracy tolerance on K     */
                1e-10);                 /* (I) Function accuracy tolerance   */
            /* check we got the right weights in the 0 correlation case */
            /* ASSERT(b.beta_r!=0. || fabs(ccmNormalCum(p.T1) - p.W1ind)<1.e-8);*/
            p->isT1Calibrated = true;
        }
    } catch (exception& e){
        throw ModelException(e, routine);
    }
}


/* ---------------------------------------------------------------------- 
CREATED BY :    Julia Tolpin

CONTAINS   :    RootFindBrent

Copyright 1995 J.P. Morgan & Co. Incorporated.   All rights reserved.
-------------------------------------------------------------------------  */

#undef ABS
#define ABS(a) ((a)>0 ? (a) : -(a))

static void secantMethod 
(TObjectFunc funcd,        /* (I) Name of function to call */
 void       *data,         /* (I) Data to pass into funcd */
 int        numIterations, /* (I) Maximum number of iterations */
 double     xacc,          /* (I) X accuracy */
 double     facc,          /* (I) F accuracy */
 double     boundLo,       /* (I) Lower bound */
 double     boundHi,       /* (I) Upper bound */
 double     *xPoints,      /* (I/O) Array of x values */
 double     *yPoints,      /* (I/O) Array of y values */
 bool   *foundIt,      /* (O) If solution was found */
 bool   *bracketed,    /* (O) If root was bracketed */
 double     *solution);    /* (O) Root of equation */

static double brentMethod 
(TObjectFunc funcd,         /* (I) Name of function to call */
 void       *data,         /* (I) Data to pass into funcd */
 int        numIterations, /* (I) Maximum number of iterations*/
 double     xacc,          /* (I) X accuracy tolerance */
 double     facc,          /* (I) Function accuracy tolerance */
 double     *xPoints,      /* (I) Array of x values */
 double     *yPoints);     /* (I) Array of y values */

/** the point of this function is that it catches exceptions and then
    propagates an exception with a more useful message */
static double callFunction
(TObjectFunc funcd,        /* (I) Name of function to call */
 double      x,            /* (I) Function argument */  
 void       *data);        /* (I) Data to pass into funcd */

static char routine[]="ccmRootFindBrent"; 


#define SWITCH(a,b)  \
{                    \
    double temp = a; \
    a = b;           \
    b = temp;        \
}               

/* 
 * This is the percentage of x=(boundHi - boundLo) 
 * used to set x[2] if x[2] is outside 
 * [boundLo, boundHi]. 
 */

#define ONE_PERCENT 0.01

#undef REQUIRE
#define REQUIRE(cond) do { if (!(cond))\
{\
    throw ModelException(routine, Format::toString (\
         "Required condition (%s) fails!", #cond));\
}} while(0)


/*---------------------------------------------------------------------
Function:   ccmRootFindBrent

Created by: Julia Tolpin

Description: Finds the root of f(x) =  0 using a combination of 
             secant,bisection and an inverse quadratic interpolation method.
-------------------------------------------------------------------------  */
double ccmRootFindBrent(           
   TObjectFunc funcd,                   /* (I) Function to call */
   void       *data,                    /* (I) Data to pass into funcd */
   double      boundLo,                 /* (I) Lower bound on legal X */
   double      boundHi,                 /* (I) Upper bound on legal X */
   int         numIterations,           /* (I) Maximum number of iterations */
   double      guess,                   /* (I) Initial guess */
   double      initialXStep,            /* (I) Size of step in x */
   double      initialFDeriv,           /* (I) Initial derivative or 0*/
   double      xacc,                    /* (I) X accuracy tolerance */
   double      facc)                    /* (I) Function accuracy tolerance */
{
   double      fLo;                    /* Function evaluated at boundLo */
   double      fHi;                    /* Function evaluated at boundHi */
   bool    bracketed;              /* If root was bracketed by secant */
   bool    foundIt;                /* If root was found by secant */
   double      xPoints[3];             /* Array of x values */
   double      yPoints[3];             /* Array of y values */
   double      boundSpread;            /* Diff between the hi and lo bounds */ 
   double      solution;

   xPoints[0] = guess;

   REQUIRE(boundLo < boundHi);
   REQUIRE(guess >= boundLo);
   REQUIRE(guess <= boundHi);

   yPoints[0] = callFunction (funcd, xPoints[0], data);

   /* Check if guess is the root (use bounds for x-accuracy) */
   if (yPoints[0] == 0.0 || 
       (ABS(yPoints[0]) <= facc && (ABS(boundLo-xPoints[0]) <= xacc ||
                                    ABS(boundHi-xPoints[0]) <= xacc)))
   {
       solution = xPoints[0];
       return solution;
   }

   /* If the initialXStep is 0, set it to ONE_PERCENT of 
    *  of (boundHi - BoundLo). 
    */
   boundSpread  = boundHi - boundLo; 
   if ( initialXStep == 0.0 ) 
   {
       initialXStep = ONE_PERCENT * boundSpread; 
   }
 
   /* Take a step of the size passed in, if the derivative is not
      passed in */
   if (initialFDeriv == 0)
   {
       xPoints[2] = xPoints[0] + initialXStep;
   }
   else
   {
   /* If initial derivative is known, use Newton's Method 
    * it to find the next point.
    */
       xPoints[2] = xPoints[0] - (yPoints[0])/initialFDeriv;
   }

   /* Now check to make sure that xPoints[2] 
    * is within the hi-lo bounds. If it isn't, then adjust it 
    * so that it is. 
    */
   if ( xPoints[2] < boundLo || xPoints[2] > boundHi ) 
   {
       /* Switch the direction of the step */ 
       xPoints[2] = xPoints[0] - initialXStep; 
       if ( xPoints[2] < boundLo )
       {
           /* xPoints[2] is still too small, so we make 
            * it boundLo
            */ 
           xPoints[2] = boundLo; 
       }
       
       if ( xPoints[2] > boundHi )
       {
           /* xPoints[2] is too large, then set it 
            * to boundHi. 
            */ 
           xPoints[2] = boundHi; 
       }

       if ( xPoints[2] == xPoints[0] )
       {
           /* We cannot have xPoints[0] and 
            * xPoints[2] be the same. 
            */ 
           if ( xPoints[2] == boundLo ) 
           {
               xPoints[2] = boundLo + ONE_PERCENT * boundSpread; 
           }
           else 
           {
               xPoints[2] = boundHi - ONE_PERCENT * boundSpread; 
           }
       }
   }
       
   /* Finally, try to call (*funcd) with xPoints[2], to make 
    * that the function can return a value at that point. 
    */
   yPoints[2] = callFunction (funcd, xPoints[2], data);

   /* Check if the new point meets the tolerance requirements */
   if (yPoints[2] == 0.0 || 
       (ABS(yPoints[2]) <= facc  && ABS(xPoints[2]-xPoints[0]) <= xacc))
   {
       solution = xPoints[2];
       return solution;
   }

   /* Call secant method to find the root, or to get a 
      third point, so that two of the three points bracket the root. */
   
   secantMethod (funcd, data, numIterations,
                 xacc, facc, boundLo, boundHi,
                 xPoints, yPoints, 
                 &foundIt, &bracketed, &solution);

   if (foundIt)
       return solution;

   if (! bracketed)
   {
       /* Root was not bracketed, now try at the bounds 
        */
       fLo = callFunction (funcd, boundLo, data);
       if (fLo == 0.0 || 
           (ABS(fLo) <= facc && ABS(boundLo - xPoints[0]) <= xacc))
       {
           solution = boundLo; 
           return solution;
       }
   
       /* If these points bracket the root, assign points so that
          xPoints[0] < xPoints[2] */
       if (yPoints[0]*fLo < 0)
       {
           xPoints[2] = xPoints[0];
           xPoints[0] = boundLo;
           
           yPoints[2] = yPoints[0];
           yPoints[0] = fLo;
           
       }
       else 
       {
           /* Root is still not bracketed, so try at the upper bound now. */
           fHi = callFunction (funcd, boundHi, data);
           if (fHi == 0.0 || 
               (ABS(fHi) <= facc && ABS(boundHi-xPoints[0]) <= xacc)) 
           {
               solution = boundHi; 
               return solution;
           }
           
           /* If points bracket the root, assign xPoints[2] to boundHi */
           if (yPoints[0]*fHi < 0)
           {
               xPoints[2] = boundHi;
               yPoints[2] = fHi;
           }
           else 
           {
               /* Root could not be bracketed at the bounds. */
               throw ModelException (routine, Format::toString (
                   "Function values (%2.6e,%2.6e) at bounds "
                   "(%2.6e, %2.6e) imply no root exists.",
                   fLo, fHi, boundLo, boundHi));
           }
       }        
       
       /* xPoints[0] and xPoints[2] bracket the root, but we need third
          point to do Brent method. Take the midpoint. */
       xPoints[1] = 0.5*(xPoints[0]+xPoints[2]);
       yPoints[1] = callFunction (funcd, xPoints[1], data);
       if (yPoints[1] == 0.0 || 
           (ABS(yPoints[1]) <= facc && ABS(xPoints[1] - xPoints[0]) <= xacc)) 
       {
           solution = xPoints[1]; 
           return solution;
       }
   }

   /* xPoints is an array of three points, two of which bracket the root.
      Call brent Method now, to find the root */
   solution = brentMethod (funcd, data, numIterations, xacc, facc,
                           xPoints, yPoints);
   return solution;
}




/* ----------------------------------------------------------------------
FUNCTION    :   brentMethod
CREATED BY  :   Julia Tolpin
DESCRIPTION :   Finds the root using a combination of inverse quadratic
                method and bisection method.
---------------------------------------------------------------------- */
static double brentMethod 
(TObjectFunc  funcd,        /* (I) Function to evaluate */
 void        *data,         /* (I) Data to pass into funcd */
 int         numIterations, /* (I) Maximum number of iterations*/
 double      xacc,          /* (I) X accuracy tolerance */
 double      facc,          /* (I) Function accuracy tolerance */
 double      *xPoints,      /* (I) Array of x values */
 double      *yPoints)      /* (I) Array of y values */
{
    int         j;                      /* Index */
    double      ratio;                  /* (x3-x1)/(x2-x1) */
    double      x31;                    /* x3-x1*/
    double      x21;                    /* x2-x1*/
    double      f21;                    /* f2-f1 */
    double      f31;                    /* f3-f1 */
    double      f32;                    /* f3-f2 */
    double      xm;                     /* New point found using Brent method*/
    double      fm;                     /* f(xm) */

    double      x1 = xPoints[0];        /* xN short hand for xPoints[n] */
    double      x2 = xPoints[1];
    double      x3 = xPoints[2];

    double      f1 = yPoints[0];
    double      f2 = yPoints[1];
    double      f3 = yPoints[2];

    double      solution;

    for (j=1; j<=numIterations; j++) 
    {
        /* Always want to be sure that f1 and f2 have opposite signs,
         * and f2 and f3 have the same sign.
         */
        if (f2*f1>0.0)
        {   
            SWITCH(x1,x3);
            SWITCH(f1,f3);
        }
        f21 = f2-f1;
        f32 = f3-f2;
        f31 = f3-f1;
        x21 = x2-x1;
        x31 = x3-x1;
        /* Check whether this is suitable for interpolation. When checking
         * for f21,f31,f32 = 0, we don't use IS_ALMOST_ZERO for efficiency
         * reasons. If the objective function has been truncated to a 
         * certain number of digits of accuracy, f21,f31,or f32 could be
         * (were in one case) zero. In this case we need to protect against
         * division by zero. So we use bisection instead of brent.
         */
        ratio = (x3-x1)/(x2-x1);
        if (f3*f31<ratio*f2*f21 || f21 == 0. || f31 == 0. || f32 == 0.) 
        {
            /* This is not suitable, do bisection 
             */
            x3 = x2;
            f3 = f2; 

        }
        else 
        {
            xm = x1 - (f1/f21)*x21 + ((f1*f2)/(f31*f32))*x31 - 
                ((f1*f2)/(f21*f32))*x21;
            fm = callFunction(funcd, xm, data);
            if (fm == 0.0 || (ABS(fm) <= facc && ABS(xm-x1) <= xacc))
            {
                solution=xm;
                return solution;
            }
            /* If the new point and point1 bracket the root,
               replace point3 with the new point */
            if (fm*f1<0.0)
            {
                x3=xm;
                f3=fm;
            }
            /* If the root is not bracketed, replace point1 with new point,
               and point3 with point2 */
            else
            {
                x1=xm;
                f1=fm;
                x3=x2;
                f3=f2;
            }
        }             
        x2 = 0.5*(x1+x3); 
        f2 = callFunction (funcd, x2, data);
        if (f2 == 0.0 || (ABS(f2) <= facc && ABS(x2 - x1) <= xacc))
        {
            solution = x2;
            return solution;
        }
    }

    throw ModelException (routine, "Maximum number of iterations exceeded");
}    
       


/* -----------------------------------------------------------------------
FUNCTION    :     secantMethod
CREATED BY  :     Julia Tolpin
DESCRIPTION :     Uses the secant method until one of the following things
                  happens. On the right, the setting of foundIt, bracketed,
                  and the return status is given. 

                                        foundIt  bracketed  status
                                        -------  ---------  ------
    1) the root is found                true     true       SUCCESS
    2) root is bracketed                false    true       SUCCESS
    3) max # steps reached              false    false      SUCCESS
    4) next x not in [boundLo, boundHi] false    false      SUCCESS
    5) Objective func returns FAILURE   false    false      FAILURE
 
    Note that the routine only returns FAILURE if the objective 
    function returns FAILURE. 
---------------------------------------------------------------------- */
static void secantMethod 
(TObjectFunc funcd,        /* (I) Function to evaluate */
 void       *data,         /* (I) data to pass into funcd */
 int         numIterations, /* (I) Maximum number of iterations */
 double      xacc,          /* (I) Accuracy of x */
 double      facc,          /* (I) Accuracy of f */
 double      boundLo,       /* (I) Lower bound */
 double      boundHi,       /* (I) Upper bound */
 double     *xPoints,      /* (I/O) Array of x points */
 double     *yPoints,      /* (I/O) Array of y points */
 bool       *foundIt,      /* (O) If solution was found */
 bool       *bracketed,    /* (O) if root was bracketed*/
 double     *solution)     /* (O) Root of function */
{
    int           j=numIterations;      /* Index */
    double        dx;                   /* Delta x used for secant */       

    *foundIt = false;           /* Until solution is found. */
    *bracketed = false;         /* Until bracketed */

    while (j--)
    {
        /* Swap points so that yPoints[0] is smaller than yPoints[2] */
        if (ABS(yPoints[0])>ABS(yPoints[2]))
        {   
            SWITCH(xPoints[0],xPoints[2]);
            SWITCH(yPoints[0],yPoints[2]);
        }

        /* Make sure that you do not divide by a very small value */
        if (ABS(yPoints[0]-yPoints[2]) <= facc)
        {
            if (yPoints[0] - yPoints[2] > 0)
            {
                dx = -yPoints[0]*(xPoints[0] - xPoints[2])/facc;
            }
            else
            {
                dx = yPoints[0]*(xPoints[0] - xPoints[2])/facc;
            }
        }
        else
        {
            dx= (xPoints[2]-xPoints[0])* yPoints[0]/(yPoints[0]-yPoints[2]);
        }
        xPoints[1] = xPoints[0] + dx;

        /* Make sure that the point is within bounds 
         */
        if (xPoints[1] < boundLo || xPoints[1] > boundHi)
            return; /* Not bracketed, not found */

        yPoints[1] = callFunction (funcd, xPoints[1], data);
        if (yPoints[1] == 0.0 || 
            (ABS(yPoints[1]) <= facc && ABS(xPoints[1] - xPoints[0]) <= xacc))
        {
            *solution  = xPoints[1];
            *foundIt   = true;
            *bracketed = true;
            return;      /* Found, bracketed */
        }

        if ((yPoints[0] < 0 && yPoints[1] < 0 && yPoints[2] < 0) ||
            (yPoints[0] > 0 && yPoints[1] > 0 && yPoints[2] > 0))
        {
            /* Swap points so that yPoints[0] is always smallest 
             */
            if (ABS(yPoints[0]) > ABS(yPoints[1]))
            {
                xPoints[2] = xPoints[0];
                yPoints[2] = yPoints[0];
                xPoints[0] = xPoints[1];
                yPoints[0] = yPoints[1];
            }
            else 
            {
                xPoints[2] = xPoints[1];
                yPoints[2] = yPoints[1];
            }
            continue;
        }
        else
        { 
            /* Root was bracketed. 
             * Swap points so that yPoints[0]*yPoints[2] < 0 
             */
            if (yPoints[0]*yPoints[2] > 0)
            {
                /* Make sure that you swap so that 
                 * xPoints[0]<xPoints[1]<xPoints[2] 
                 */
                if (xPoints[1] < xPoints[0])
                {
                    SWITCH(xPoints[0], xPoints[1]);
                    SWITCH(yPoints[0], yPoints[1]);
                }
                else
                {
                    SWITCH (xPoints[1], xPoints[2]);
                    SWITCH (yPoints[1], yPoints[2]);
                }
            }
            /* Root was bracketed, but not found.
             */
            *bracketed = true;
            return;
        }
    } /* while */
    
    /* Root not bracketed or found.
     */
    return;
}

static double callFunction
(TObjectFunc funcd,        /* (I) Name of function to call */
 double      x,            /* (I) Function argument */  
 void       *data)         /* (I) Data to pass into funcd */
{
    double y;
    try
    {
        y = (*funcd)(x,data);
    } catch (exception &e) {
        throw ModelException (e, routine, Format::toString (
            "Objective function failed at %2.6e", x));
    }
    return y;
}

DRLIB_END_NAMESPACE


