/*********************************************************************************
 * pricer based on matrix
 *
 ********************************************************************************/

#ifndef _COMPONENT_H_
#define _COMPONENT_H_

#ifdef __cplusplus
extern "C" {
#endif

#include "transition.h"

typedef enum{
    UNCOND,
    COND_NO_DEFAULT,        
}SettleConv;
    
/* snapshot of the vector at a certain time                                     */
typedef struct{
    double         spd[MAX_NB];           /* par spread                         */
    double         annuity[MAX_NB];       /* annuity                            */
    double         prot[MAX_NB];          /* protection                         */
    double         prob[MAX_NB];          /* prob vec from L=0,t=0->L=i,t=time  */
}MatState;

/*********************************************************************************
 * Set Tranche Payoff Vector
 *
 ********************************************************************************/
int SetTranchePayoff(
    double          *payoff,              /* (O) tranche contingent payoff vec  */
    double          lowerSize,            /* (I) lower strike                   */
    double          upperSize,            /* (I) upper strike                   */
    double          notional,             /* (I) notional of the payoff         */
    EffCurves       *effCurves);          /* (I) effective curves               */
    
/*********************************************************************************
 * protection leg
 *
 ********************************************************************************/
int MatrixProtection(
    double          *price,               /* (O) price of the payoff matrix     */
    TDate           valueDate,            /* (I) value date                     */
    TDate           today,                /* (I) today                          */
    TDate           startDate,            /* (I) protection start date          */
    TDate           endDate,              /* (I) protection end date            */
    double          *prot,                /* (I) protection payoff vector       */
    EffCurves       *effCurves,           /* (I) effective curves               */    
    TCurve          *irCurve);            /* (I) ir curve                       */

/*********************************************************************************
 * annuity pricer based on transition matrix
 *
 ********************************************************************************/
int  MatrixAnnuity(
    double          *price,               /* (O) price of the payoff matrix     */
    TDate           valueDate,            /* (I) value date                     */
    TDate           today,                /* (I) today                          */
    TDate           issueDate,            /* (I) issueDate                      */
    TDate           maturityDate,         /* (I) issueDate                      */
    double          cpnRate,              /* (I) coupon rate                    */
    TDateInterval   cpnInterval,          /* (I) coupon interval                */
    long            DCC,                  /* (I) coupon DCC                     */
    double          lowerSize,            /* (I) lower strike                   */
    double          upperSize,            /* (I) upper strike                   */
    double          notional,             /* (I) notional                       */
    int             payPrincipal,         /* (I) =1, pay principal at maturity  */
    EffCurves       *effCurves,           /* (I) effective curves               */    
    TCurve          *irCurve);            /* (I) ir curve                       */

/*********************************************************************************
 * fwd par spread based on transition matrix
 *
 ********************************************************************************/
int  MatrixParSpread(
    MatState        *res,                 /* (O) mat state                      */
    TDate           valueDate,            /* (I) value date                     */
    TDate           today,                /* (I) today                          */
    TDate           issueDate,            /* (I) issueDate                      */
    TDate           maturityDate,         /* (I) maturityDate                   */
    TDateInterval   cpnInterval,          /* (I) coupon interval                */
    long            DCC,                  /* (I) coupon DCC                     */
    double          lowerSize,            /* (I) lower strike                   */
    double          upperSize,            /* (I) upper strike                   */
    EffCurves       *effCurves,           /* (I) effective curves               */    
    TCurve          *irCurve);            /* (I) ir curve                       */

/*********************************************************************************
 * pricer based on transition matrix
 *
 ********************************************************************************/
int  MatrixPrice(
    double          *price,               /* (O) price of the payoff matrix     */
    TDate           valueDate,            /* (I) value date                     */
    TDate           today,                /* (I) today                          */
    PayOff          *pPayoff,             /* (I) payoff matrix                  */
    EffCurves       *effCurves,           /* (I) effective curves               */    
    TCurve          *irCurve);            /* (I) ir curve                       */


/*********************************************************************************
 * fwd par spread distribution based on transition matrix
 *
 ********************************************************************************/
int  FwdSpdDistribution(
    MatState        *res,                 /* (O) mat state                      */
    TDate           fwdDate,              /* (I) fwd date                       */
    TDate           today,                /* (I) today                          */
    TDate           issueDate,            /* (I) issueDate                      */
    TDate           maturityDate,         /* (I) issueDate                      */
    TDateInterval   cpnInterval,          /* (I) coupon interval                */
    long            DCC,                  /* (I) coupon DCC                     */
    double          lowerSize,            /* (I) lower strike                   */
    double          upperSize,            /* (I) upper strike                   */
    EffCurves       *effCurves,           /* (I) effective curves               */    
    TCurve          *irCurve);            /* (I) ir curve                       */

/*********************************************************************************
 * Insert payoff at a certain pos
 *
 ********************************************************************************/
int PayoffInsert(
    PayOff    *pPayoff,                   /* (I/O) payoff                       */
    int       pos,                        /* (I)   pos to insert into           */
    TDate     date,                       /* (I)   payoff date                  */
    int       numName);                   /* (I)   num of tranchelets           */

/*********************************************************************************
 * extend payoff to contain exercise dates
 *
 ********************************************************************************/
int ExtendPayoff(
    PayOff     *pPayoffNew,               /* (O) extended payoff                */
    int        *exerFlag,                 /* (O) exer flag, 0: not exer, 1: yes */
    PayOff     *pPayoffOld,               /* (I) original payoff                */
    TDate      *pDate,                    /* (I) exercise dates                 */
    int         numDate,                  /* (I) num. of exer. dates            */
    int         numName);                 /* (I) num of names                   */

#ifdef __cplusplus
}
#endif

#endif
