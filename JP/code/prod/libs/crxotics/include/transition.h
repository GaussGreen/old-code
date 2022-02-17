/*********************************************************************************
 * Transition matrix for defaults
 *
 ********************************************************************************/

#ifndef _TRANSITION_H_
#define _TRANSITION_H_

#ifdef __cplusplus
extern "C" {
#endif

#include "tcurve.h"
#include "gtobf.h"
#include "datelist.h"
#include "ldate.h"
#include "convert.h"
#include "ctype.h"
    
#define MAX_NB 131                        /* make it larger than 125            */

#ifndef TOL
#define TOL   1e-4
#endif

#ifndef SMALL
#define SMALL   1e-14
#endif

#define MAT_INTERP GTO_FLAT_FORWARDS

typedef struct{
    int size;
    int size2;
    double ptr[MAX_NB * MAX_NB];
}Mat;

typedef struct{
    TDate               today;            /* base date                          */
    int                 numName;          /* size of the matrix                 */
    int                 numTotal;         /* total number of names              */
    TCurve              *cCurve;          /* catastrophe curve                  */
    TDateInterval       step;             /* time step                          */
    double              recovery;         /* recovery                           */

    double              *loss;            /* loss                               */
    TCurve              **pCurves;        /* effective curves                   */
}EffCurves;
    
typedef struct{
    int     numDate;
    TDate   cutoffDate;                      /* protection cutoff date */
    TDate   payDate[MAX_NB];
    double  payoff[MAX_NB][MAX_NB];
}PayOff;

typedef struct{
    double  price[MAX_NB];
    double  prob[MAX_NB];
    double  fwdPrice[MAX_NB];
}PriceOutPut;

typedef struct{
    int     nbVol;
    TDate   volDate[MAX_NB];
    TDate   swapSt[MAX_NB];
    TDate   swapMat[MAX_NB];
    double  vol[MAX_NB];
    char    volType[MAX_NB];
}VolInput;

typedef struct{
    double  alpha;                        /* alpha                              */
    double  fac[MAX_NB];                  /* conditional factor                 */
    TDate   facDate[MAX_NB];              /* date corresponding to each factor  */
    int     numDate;                      /* number of fac dates                */
    
    double  sum1[MAX_NB];                 /* zero moment of exponential tail    */
    double  sum2[MAX_NB];                 /* first moment of exponential tail   */
    double  tail[MAX_NB];                 /* tail component                     */
    
    double  kFac;                         /* rotation factor                    */
    double  decay;                        /* decay factor                       */
    int     numName;                      /* matrix dimension                   */
    int     numTotal;                     /* total number of names              */    
    int     numT;                         /* number of factors in time dimension*/
    TDate   baseDate;                     /* base date                          */
    TDateInterval step;                   /* time step used to calibrate        */
    TCurve  *portCurve;                   /* portfolio spread curve             */
}DTMParm;    

typedef struct{
    DTMParm parms[15];
}DTMParms;
    
typedef struct{
    int     numDate;                      /* number of target dates             */        
    int     numStrike;                    /* number of target loss              */    
    TDate   maturityDate[MAX_NB];         /* maturity date                      */
    int     index[MAX_NB];                /* index in loss dimension            */
    double  strikes[MAX_NB];              /* upper sizes                        */
    double  loss[MAX_NB][MAX_NB];         /* expected loss                      */

    double  lossDist[MAX_NB][MAX_NB];     /* calib. loss dist at target dates   */
    double  recovery;
    int     numTotal;
    
}TargetLoss;
    
/*********************************************************************************
 * default transition matrix from startDate to endDate
 *
 ********************************************************************************/
int  TransMatBK(
    Mat             *mat,                 /* (O) transition matrix              */
    double          *Gamma,               /* (O) gamma                          */    
    int             numName,              /* (I) num of tranchelets             */
    TCurve          **pCurve,             /* (I) eff curves for tranchelets     */
    TDate           startDate,            /* (I) start Date                     */
    TDate           endDate,              /* (I) end Date                       */
    double          slope,                /* (I) conditional slope              */
    double          beta,                 /* (I) correlation                    */
    long            curveInterp,          /* (I) curve Interp method            */
    TDateInterval   step,                 /* (I) date step                      */
    int             cName);               /*                                    */
    
/*********************************************************************************
 * default transition matrix from startDate to endDate
 *
 ********************************************************************************/
int  TransMatExp(
    Mat             *mat,                 /* (O) transition matrix              */
    TDate           startDate,            /* (I) start Date                     */
    TDate           endDate,              /* (I) end Date                       */
    EffCurves       *effCurves);          /* (I) effective curves               */

/*********************************************************************************
 * default transition matrix from startDate to endDate
 *
 ********************************************************************************/
int  TransMat(
    Mat             *mat,                 /* (O) transition matrix              */
    TDate           startDate,            /* (I) start Date                     */
    TDate           endDate,              /* (I) end Date                       */
    EffCurves       *effCurves);          /* (I) effective curves               */

/*********************************************************************************
 * default transition matrix kernel from startDate to endDate
 *
 ********************************************************************************/
int  TransMatExpKernel(
    Mat             *mat,                 /* (O) transition matrix              */
    TDate           startDate,            /* (I) start Date                     */
    TDate           endDate,              /* (I) end Date                       */
    EffCurves       *effCurves);          /* (I) effective curves               */

/*********************************************************************************
 * Calculate residual expected loss
 *
 ********************************************************************************/
int ResidualLoss(
    double          *loss,                /* (O) residual expected loss         */
    TDate           date,                 /* (I) a specific date to observe     */
    EffCurves       *effCurves);          /* (I) effective curves               */

/*********************************************************************************
 * Calibrate DTM to a given loss distritubion
 * 
 *
 ********************************************************************************/
int  DTMCalib(
    DTMParms        *dtmParm,             /* (O) DTM param structure            */
    TDate           today,                /* (I) today                          */
    TargetLoss      *target,              /* (I) target tranche structure       */
    TCurve          *portCurve,           /* (I) average spread                 */
    double          beta,                 /* (I) CM correlation                 */
    double          kFac,                 /* (I) tail decay                     */
    EffCurves       *effCurves);          /* (I) effective curves               */    
    
/*********************************************************************************
 * Calc DTM given dtmParm
 * 
 *
 ********************************************************************************/
int  DTMCalc(
    Mat             *mat,                 /* (O) transition matrix              */
    DTMParms        *dtmParm,             /* (O) DTM param structure            */        
    TDate           startDate,            /* (I) today                          */
    TDate           endDate);             /* (I) target tranche mat date        */
    
#include "util.h"
#include "component.h"
#include "lapackinterface.h"
#include "smile.h"
    
#ifdef __cplusplus
}
#endif

#endif
