/****************************************************************************/
/*      Declaration file for structure of data.                             */
/****************************************************************************/
/*      FIX123.h                                                            */
/****************************************************************************/

/*
$Header$
*/

/* Use safe multiple inclusion for Jerry Cohen's stream libraries */

#ifndef _fix123_h
#define _fix123_h

#include "esl.h"

/* TREE DATA STRUCTURE */
typedef struct _TREE_DATA
{
    /* Critical dates */
    CRIT_DATE   *CritDate[NBCRITDATE];     /* Critical dates description     */
    char        CritType[NBCRITDATE];      /* Type of critical dates         */
    int         *TPtype[NBCRITDATE];       /* Critical type of current TP    */
    int         NbZeros[NBCRITDATE];       /* Number of zeros in zero bank   */

    /* Express DEV dates */ 
    int         NbEDevDates;               /* Nb of dates for express DEV    */               
    long        *EDevDate;                 /* Dates where express DEV req    */
    double      **EDevStPrice;             /* Corresp state prices           */

    /* Time points */
    int         NbTP;                      /* Nb of time points in the tree  */
    long        *TPDate;                   /* Date of each time point        */
    double      *Length;                   /* Length of time steps (ACT/365) */
    int         Ppy;                       /* Number of time point per year  */
    int         JumpPpy;                   /* Used to calculate CET jump size*/
    int         NbDailyPts;                /* Number of daily points         */


    /* Zero curve */
    double      *ZeroCoupon[3];            /* Zero price at time point       */
    double      *ZeroRate[3];              /* Zero rate at time point        */
    double      *FwdRate[3];               /* One period forward at TP       */

    /* Internal assigment of zero curve */
    int         CvDiff;                    /* Diffused curve                 */
    int         CvIdx1;                    /* First index curve              */
    int         CvIdx2;                    /* Second index curve             */
    int         CvDisc;                    /* Discount curve                 */ 

    /* Model and volatility */
    int         NbFactor;                  /* Number of factors              */
    double      *Aweight[6];               /* Weights at time point          */

    /* Tree geometry */
    int         NbSigmaMax;                /* Nb of std dev to cut the tree  */
    int         Width[3];                  /* Ellipsoid width                */
    int         HalfWidth[3];              /* Ellipsoid half width           */
    double      *ZCenter;                  /* Center of the tree in X-space  */
    double      *LengthJ;                  /* Time step for jump size        */
    int         *Top1, *Bottom1;           /* Limits of 1D tree              */
    int         **Top2, **Bottom2;         /* Limits of 2D tree              */
    int         ***Top3, ***Bottom3;       /* Limits of 3D tree              */
    int         *OutTop1, *OutBottom1;     /* Outer limits                   */
    int         **OutTop2, **OutBottom2;   
    int         ***OutTop3, ***OutBottom3; 

} FIX3_TREE_DATA;



/* DEV DATA STRUCTURE */
typedef struct _DEV_DATA
{
    /* Node branching */
    int     *Shift1;
    int     *Shift2;
    int     *Shift3;

    /* Discount factors */
    double  *Discount[3];

    /* Probabilities */
    double  *pu, *p0, *pd;      /* 1D */
    double  *quu, *qu0,	*qud;   /* 2D */
    double  *q0u, *q00,	*q0d;
    double  *qdu, *qd0,	*qdd;
    double  *ru, *r0, *rd;      /* 3D */

    double  *NewPrice;          /* Auxiliary slice */

} FIX3_DEV_DATA;


#endif /* _fix123_h */
