/****************************************************************************/
/*      Declaration file for structure of data.                             */
/****************************************************************************/
/*      CUPSMODL.h                                                          */
/****************************************************************************/


/*
$Header$
*/

#include "esl.h"

#ifndef _CUPSMODL_H
#define _CUPSMODL_H


#define         MAXNBPOINTS     1000
#define         MAXNBEQDIV      15000        /* Maximum number of stock dividends                  */
#define         MAXMAT          60           /* Maximum swap maturity 	                           */
#define         MAXCORRELATION  .95          /* Absolute value of maximum allowable correlation    */
#define         NBSPLINE        1000          /* Number of SPLINE points in grid interpolation      */

/* Tree types supported by the HYB3 library  */
#define         TTYPE_1IR         -1         /* One interest rate, only for CET purposes */
#define         TTYPE_2IR          0         /* Two interest rates, one in CUPS mode               */
#define         TTYPE_FX2IR        1         /* Two interest rates (one in CUPS mode) plus FX      */
#define         TTYPE_EQD2IR       2         /* Two IR's (one CUPS) plus a domestic equity         */
#define         TTYPE_EQF2IR       3         /* Two IR's (one CUPS) plus a foreign equity          */
#define         TTYPE_EQC2IR       4         /* Two IR's (one CUPS) plus a composite equity        */
#define         TTYPE_EQ1IR        5         /* One IR plus an equity                              */
#define         TTYPE_EQDFX2IR     6
#define         TTYPE_EQFFX2IR     7
#define         TTYPE_1IR2F        8         /* One IR with two factor, only for CET               */
#define         TTYPE_2IR2F1D      9         /* Two IR, 2 factors for. in CUPS mode, 1 factor Dom. */

/* Discounting modes suported at run-time */
#define         DISC_1D_NOCUPS       1       /* DEV in 1-D using the simple foreign IR, no CUPS    */
#define         DISC_2D_CUPS         2       /* DEV in 2-D using the currency protected foreign IR */
#define         DISC_3D_CUPS         3       /* DEV in 3-D using the currency protected foreign IR */
#define         DISC_2D_NOCUPS       4       /* DEV in 2-D using the IR on the 1st dimension       */
#define         DISC_2D_1IR2F_NOCUPS 5       /* DEV in 2-D using the for. IR with 2 fact., no CUPS */
#define         DISC_3D_2IR2F1D_CUPS 6       /* DEV in 3-D with domestic IR on dim. 3, */
                                             /*     for 2 factor CUPS foreign IR                   */

/* FX Bootstrap modes */
#define         FX_NO_FAILURE_ALLOWED   0L
#define         FX_USE_LAST_LEVEL       1L
#define         FX_CONSTANT_SPOT_VOL    2L
#define         FXSPOT_CUTOFF_R   10.0        /* Minimum level (ratio) of FX spot vol for tree      */

/* EQ Bootstrap modes */
#define         EQ_NO_FAILURE_ALLOWED   0L
#define         EQ_USE_LAST_LEVEL       1L
#define         EQ_CONSTANT_SPOT_VOL    2L
#define         EQSPOT_CUTOFF_R   10.0        /* Minimum level (ratio) of FX spot vol for tree      */

#define         FOR               0          /* arrary indices or multi currency arrays            */
#define         DOM               1
#define         CORR_NOT_AVAIL    -999.0     /* Value used tby Kapital to indicate N/A correlation */
#define         DR_SHRT_MIN     -32768
#define         DR_SHRT_MAX      32767
#define         A2CUTOFF          1E-4 
#define         MAX_FX_ITER           2      /* Number iterations for FX drift search               */
#define         DOMZERO_TOL           0.00005
#define         FWDFX_TOL             0.0003  /* previously one less zero*/
#define         FWDEQ_TOL             0.0001  /* this is a FRACTIONAL tolerance */
#define         INNER_TREE_SIZE       1.0    /* */
#define         IR_MINSIZE            3.0
#define         STDDEV_FACTOR         0.3    /* Coeff used to determine IR NbStdev                  */

#define         NBFXBARR          50         /* Specific to the binary turbo suite     */
                                             /* (number of FX levels in barr spec.)     */
#define         BARRIER_ADJ     0.58259716   /* Barrier adjustment for continuous knock-out freq   */



/* MARKET VOLATILITY CURVE INPUT DATA STRUCTURE */

typedef struct                  /* STOCK INPUT DATA STRUCTURE                */
{
    long ValueDate;             /* eq value date                             */

    char SettleType;            /* Settlement type: 'F'ixed or 'R'olling     */
    char FwdType[MAXNBEQDIV];   /* Dividend / forward type                   */

    int  NbFwd;                 /* Number of dividend / forward prices       */
    int  NbSettle;              /* Number of settlement dates(<0 if rolling) */
    int  NbVol;                 /* Number of points in volatility curve      */
    int  NbBorrow;              /* Number of points in borrowing cost curve  */
    int  NbCredit;              /* Number of points in credit spread curve   */
    
    long  FwdDate[MAXNBEQDIV];    /* Dividends / forwards dates              */	
    long  LastTrading[MAXNBDATE]; /* Last trading date of settlement period  */
    long  SettleDate[MAXNBDATE];  /* Corresponding settlement date           */
    long  VolDate[MAXNBDATE];     /* Volatility dates                        */
    long  BorrowDate[MAXNBDATE];  /* Borrowing curve dates                   */
    long  CreditDate[MAXNBDATE];  /* Credit spread curve dates               */

    long    InpSpotVolDate[MAXNBDATE+1];
    int     NbInpSpotVol;

    double  Spot;                 /* Spot value of the stock                 */
    double  Fwd[MAXNBEQDIV];      /* Dividend / forward prices               */
    double  Vol[MAXNBDATE+1];     /* Stock volatility curve                  */
    double  InpSpotVol[MAXNBDATE+1];
    double  Rho[3][MAXNBDATE+1];  /* Corr's Eq x IrFor, Eq x IrDom, Eq x FX  */
    double  Borrow[MAXNBDATE];    /* Borrowing cost zero rate (ann ACT/365F) */
    double  Credit[MAXNBDATE];    /* Credit spread (bp spread to Libor curve)*/

    /* Calibration flags */
    int     EqCutOffFlag;            /* True if cutoff allowed when calib fails */
    int     EqCutOffLast;            /* True if cutoff at last spotvol level    */
    int     EqBootStrapMode;         /*                                         */
    double  EqCutOffLevel;           /* Index calibration flag                  */

    /* Input Smile Parameters */
    double  a1[MAXNBDATE];           /* Skew parameters               */
    double  a2[MAXNBDATE];           /* convexity parameter           */
    double  a3[MAXNBDATE];           /* Controls Max Slope of Bps vol */
    long    SmileDate[MAXNBDATE];    /* dates for time dependent smile term-structure*/
    int     NbSmilePt;    
} EQ_DATA;




typedef struct                       /*    FX INPUT DATA STRUCTURE              */
{                                                 
    /* The FX spot as of today */    
    double  Spot; 
        
    /* Dates */                 
    long    Today;                   /* Today's date                            */
    long    ValueDate;               /* Output: fx value date                   */
    long    VolDate[MAXNBDATE+1];    /* Benchmark option expiration dates       */
    long    InpSpotVolDate[MAXNBDATE+1];
    int     SpotDays;                /* Nb of spot days (today to value date)   */
    int     NbVol;                   /* Nb of points in fx volatility curve     */
    int     NbInpSpotVol;
                 
    /* Vols and correlations */
    double  FxVol[MAXNBDATE+1];      /* Benchmark option volatilities           */
    double  InpSpotVol[MAXNBDATE+1];
    double  Rho[3][MAXNBDATE+1];     /* Correlations (IRvsIR, IRvsFX, IRvsFX)   */

    /* Calibration flags */
    int     FxCutOffFlag;            /* True if cutoff allowed when calib fails */
    int     FxCutOffLast;            /* True if cutoff at last spotvol level    */
    int     FxBootStrapMode;         /*                                         */
    double  FxCutOffLevel;           /* Index calibration flag                  */

    /* Input Smile Parameters */
    double  a1[MAXNBDATE];           /* Skew parameters               */
    double  a2[MAXNBDATE];           /* convexity parameter           */
    double  a3[MAXNBDATE];           /* Controls Max Slope of Bps vol */
    long    SmileDate[MAXNBDATE];    /* dates for time dependent smile term-structure*/
    int     NbSmilePt;    
} FX_DATA;




/*************************/ 
/*     SMILE CACHING     */
/*************************/
typedef struct 
{
    /* Forward FX or Equity mapping function */
    /* the whole setup is not very elegant neither very save but it's the best     */
    /* that can be done using old fashioned C: needs rewriting using better design */

    double  **X;               /* for tabulated K^(-1) map: normalised FX values   */
    double  **K;               /* for tabulated K^(-1) map: corresponding K values */
    double  **SPL;             /* for tabulated K^(-1) map: corresponding 2nd deriv*/
    double  **SPL_Inv;         /* for tabulated K map: corresponding 2nd deriv     */
    int     nbCachePts;       
    int     *nbPtSPL;

    /* for tabulated gDash and kDashtimesX function for FX smile drift */
    double  **gd, **gd_SPL;    /* like K & SPL: tabulated gDash-function        */
    double  **kdX, **kdX_SPL;  /* like K & SPL: tabulated kDashtimesX-function  */
    int     isCached;           /* flag to indicate if there is a cache present  */

} HYB3_FXSMILE_CACHE;



/*************************/ 
/*     TREE MODELLING    */
/*************************/
typedef struct                  /* TREE BUILDING DATA STRUCTURE              */
{   

    /* Time points */
    int      NbTP;                     /* Total number of nodes              */
    long    *TPDate;                   /* Date of each node                  */
    int      Ppy;                      /* Number of period per year          */
    double  *Length;                   /* Length of time steps (ACT/365)     */
    double  *LengthJ;                  /* Length of steps for jump size calc */


    /* Critical dates */    
    CRIT_DATE   *CritDate[NBCRITDATE]; /* Critical dates description         */  
    char         CritType[NBCRITDATE]; /* Type of critical dates             */
    int         *TPType[NBCRITDATE];   /* Critical type of current node      */
    int          NbZeros[NBCRITDATE];  /* Nb of zeros carried back in the    */
                                       /* tree in relation to this CritDate  */

    /* Express DEV dates */ 
    int         NbEDevDates;               /* Nb of dates for express DEV    */               
    long        *EDevDate;                 /* Dates where express DEV req    */
    TSLICE      *EDevStPrice;              /* Corresp state prices           */

    
    /* Deterministic zero curve */
    double  *ZeroCoupon[2][3];  /* ZC mat at each node 2 CCY's X 3 TYPES     */
    double  *ZeroRate[2][3];    /* Z rate at each node 2 CCY's X 3 TYPES     */
    double  *FwdRate[2][3];     /* Fwd rate from one period to the following */
                                /* NB: CCY's are [0]=foreign, [1]= domestic  */

    /* Internal assigment of zero curve */
    int      CvDiff[2];         /* Diffused curve                                 */
    int      CvIdx1[2];         /* First index curve                              */
    int      CvIdx2[2];         /* Second index curve                             */
    int      CvDisc[2];         /* Discount curve                                 */ 

                                
    /* IR Volatility */
    double  *SpotVol[2];        /* Instant interest rate vol at each node         */
    double  *IrAweight[2][3];   /* Instant IR A weight at each node (2 fact. mode */
    double  *Rho[6];            /* Instant correlation between all processes      */
    char     Index[2][MAXINDEX];/* Index to calibrate                             */


    /* Equity */
    double  *FwdEq;             /* Forward stock at each node in the tree         */
    double  *EqMidNode;         /* Center of the equity tree                      */
    double  *NodeSettleTime;    /* Time from current node to eq settlement        */
    long    *NodeSettleDate;    /* Equity settlmt date corresp to curr node       */
    double  *SpotEqVol;         /* Instantaneous volatility of stock              */
    double  *EqVol;             /* Equity option volatility at each node          */

    /* FX */                    
    double  *FwdFx;             /* Deterministic Forward FX values alongside the tree    */
    double  *FxMidNode;         /* Centre of the FX tree                                 */
    double  *SpotFxVol;         /* Instantaneous volatility of the FX                    */
    double  *FxVol;             /* FX volatility (average) alongside the tree            */

    /* Smile (FX or Equity) */
    double  *A1;                /* Skew parameter (as input by user ,[on specific dates])*/
    double  *A2;                /* convexity parameter (as input by user)                */
    double  *A3;                /* scale parameter (as input by user)                    */
    double  *A1C;               /* A1 coeffs on the complete TimeLine                    */ 
    double  *A2C;               /* A2 coeffs on the complete TimeLine                    */
    double  *A3C;               /* A3 coeffs on the complete TimeLine                    */


    /* caching structure for tabulated KdashX, K and g functions */
    HYB3_FXSMILE_CACHE  FXsmileCache;

    /* Forward FX or Equity mapping function */
    int     *SmileIndex;        /* Index for time dependent smile coeffs         */

    long    *tMin;
    long    *tMax;

    int      CalcCheckSlices;   /* calculate slices to check Fwd FX in lattice.c */

    
    /* Model */   
    int      TreeType;          /* One of:TTYPE_2IR,_FX2IR,_EQD2IR,_EQF2IR   */
    double  *Aweight[10];        /* Orthogonal weights at each node           */
    double  *DriftCUPS[2];         /* Drift adjustment to the foreign zero bond */
    double  *IrZCenter[2];      /* Centre of IR dimensions in the tree       */

    
    /* Tree geometry */
    int      NbSigmaMax;        /* Nb of std devs used when cutting the tree */
    double   NbSigmaDBL;         /* Nb of std devs used when cutting the tree */
    double   NbIRSigmaMax;      /* Nb of std devs used for IR dims           */
    int      Width[4];          /* Maximum number of nodes in each dimension */
    int      HalfWidth[4];      /* Minimum number of nodes in each dimension */
    int     *Top1;              /* Limits of tree in its 1st dimension       */
    int     *Bottom1;           /* Limits of tree in its 1st dimension       */
    int    **Top2;              /* Limits of tree in its 2nd dimension       */
    int    **Bottom2;           /* Limits of tree in its 2nd dimension       */
    int   ***Top3;              /* Limits of tree in its 3rd dimension       */
    int   ***Bottom3;           /* Limits of tree in its 3rd dimension       */
    int  ****Top4;              /* Limits of tree in its 4th dimension       */
    int  ****Bottom4;           /* Limits of tree in its 4th dimension       */

    /* These dimensions are obsolete and not currently used*/
    int     *OutTop1;           /* Outer limits of tree in its 1st dimension */
    int     *OutBottom1;        /* Outer limits of tree in its 1st dimension */
    int    **OutTop2;           /* Outer limits of tree in its 2nd dimension */
    int    **OutBottom2;        /* Outer limits of tree in its 2nd dimension */
    int   ***OutTop3;           /* Outer limits of tree in its 3rd dimension */
    int   ***OutBottom3;        /* Outer limits of tree in its 3rd dimension */
    int  ****OutTop4;           /* Outer limits of tree in its 4th dimension */
    int  ****OutBottom4;        /* Outer limits of tree in its 4th dimension */

    /* We only do FwdFX forced calibration ( i.e search) inside the Inner Tree */
    int      InnerTreeWidth[4];          /* Maximum number of nodes in each dimension */
    int      InnerTreeHalfWidth[4];      /* Minimum number of nodes in each dimension */

    int     *InnerTreeTop1;              /* Limits of tree in its 1st dimension       */
    int     *InnerTreeBottom1;           /* Limits of tree in its 1st dimension       */
    int    **InnerTreeTop2;              /* Limits of tree in its 2nd dimension       */
    int    **InnerTreeBottom2;           /* Limits of tree in its 2nd dimension       */
    int   ***InnerTreeTop3;              /* Limits of tree in its 3rd dimension       */
    int   ***InnerTreeBottom3;           /* Limits of tree in its 3rd dimension       */
    int  ****InnerTreeTop4;              /* Limits of tree in its 4th dimension       */
    int  ****InnerTreeBottom4;           /* Limits of tree in its 4th dimension       */

    /* Arrays containing the information for the moment matching procedure  */
    int    FxMomentMatching;             /* launches the CUPS drift adjustment calculation */

    /* temporary hyb4 support */
    int      xT;
    long     xDate;
    int      xWidth[4];          /* Maximum number of nodes in each dimension */
    int      xHalfWidth[4];      /* Minimum number of nodes in each dimension */

} HYB3_TREE_DATA;


typedef struct                        /* DEV DATA STRUCTURE */
{                            
     /* Node shifts */
    int    *Shift1;                   /* Node shift in 1st factor (IR1 CUPS) */
    int    *Shift2;                   /* Node shift in 2nd factor (IR2)      */
    int    *Shift3;                   /* Node shift in 3rd factor (Eq or FX) */
    int    *Shift4;                   /* Node shift in 1st factor, NO CUPS   */
    int    *Shift5;                   /* Node shift in 2st dim, (2nd factor IR1, NO CUPS)   */
    
    TSLICE     FxSpot;                /* Spot FX rate at each node           */
    TSLICE     NextFxSpot;
    TSLICE     EqSpot;                /* Equity spot price at each node      */
    TSLICE     NextEqSpot;
    TSLICE     EqFwd;                 /* Equity fwd price(NodeSettleTime fwd)*/

    /* auxiliary variables to speed up the local smile calculations */
    TSLICE     gDash;                 /* needed for Ito drift calculations   */
    TSLICE     kDashTimesX;           /* needed for Ito drift calculations   */
    TSLICE     kVar;                  /* for calculating deltaK=K_t(FX)-K_t-1(FX) */

  
    TSLICE     Discount_1D[3];        /* I, C and R discount for foreign ccy */
    TSLICE     Discount_2D[3];        /* I, C and R discount for domestic ccy*/
    TSLICE     Discount_3D[3];        /* 3 dim. dom. disc. for 2IR-2F1D mode  */


    TSLICE      Aux1D;                /* Auxiliary slices used in Hyb3_Dev()      */
    TSLICE      Aux2D;
    TSLICE      Aux3D;

    TSLICE      DomZero;              /* Slices used to calculate check quantities*/
    TSLICE      FwdFX;

    TPROB_0 *p;
    TPROB_0 *s;
    TPROB_0 *q;
    TPROB_0 *t;
    TPROB_0 *r;
             
    TSLICE   quu,                     /* 2-D probabilities (domestic ccy)    */
             qu0,
             qud,
             q0u,
             q00,
             q0d,
             qdu,
             qd0,
             qdd;
             
    TSLICE   t2,
             t3;
             
} HYB3_DEV_DATA;



/*************************/ 
/*   MARKET ENVIRONMENT  */
/*************************/

/* these structures are based on the format of the market data in the */
/* various *.dat files                                                */

typedef struct FXVOLATILITY_DATA
{
    long    ValueDate;
    double  FXSpotRate;
    char    BaseVolFreq;

    int     NbBaseVols;
    long    BaseVolDates[MAXNBDATE];
    double  BaseVols[MAXNBDATE];

    int     NbSpotVols;
    long    SpotVolDates[MAXNBDATE];
    double  SpotVols[MAXNBDATE];

} FXVOLATILITY_DATA;

typedef struct EQVOLATILITY_DATA
{
    long    ValueDate;
    double  EQSpotRate;

    int     NbBaseVols;
    long    BaseVolDates[MAXNBDATE];
    double  BaseVols[MAXNBDATE];

    int     NbSpotVols;
    long    SpotVolDates[MAXNBDATE];
    double  SpotVols[MAXNBDATE];

} EQVOLATILITY_DATA;


typedef struct CORRELATION_DATA
{
    double CorrIR;
    double CorrForIRFX;
    double CorrDomIRFX;

} CORRELATION_DATA;


typedef struct FXSMILE_DATA
{
    int     NbSmileDates;
    long    SmileDates[MAXNBDATE];
    double  A1[MAXNBDATE];
    double  A2[MAXNBDATE];
    double  A3[MAXNBDATE];

} FXSMILE_DATA;


typedef struct MODELPARAMETERS_DATA
{
    int    nbFactors;

    /* one factor */
    double OneFactorMR;
    double OneFactorVol;
    int    OneFactorPPY;

    /* two factor */
    double TwoFactorMR1;
    double TwoFactorMR2;
    double TwoFactorVol1;
    double TwoFactorVol2;
    double TwoFactorCorr;
    int    TwoFactorPPY;

    /* three factor */
    double ThreeFactorMR1;
    double ThreeFactorMR2;
    double ThreeFactorMR3;
    double ThreeFactorVol1;
    double ThreeFactorVol2;
    double ThreeFactorVol3;
    double ThreeFactorCorr12;
    double ThreeFactorCorr13;
    double ThreeFactorCorr23;
    int    ThreeFactorPPY;

    /* smile parameters */
    double  QLeft;
    double  QRight;
    double  FwdShift;
    int     CetNbIter;

} MODELPARAMETERS_DATA;


/* hyb3 market data structure is a container for the various market structures  */
typedef struct HYB3_MARKET_DATA
{
    SWAPVOL_DATA  FSwapVol;
    SWAPVOL_DATA  DSwapVol;
    BASEVOL_DATA  FBaseVol;
    BASEVOL_DATA  DBaseVol;

    int NbFZeroCurves, NbDZeroCurves;
    T_CURVE       FZeroCurves[3];
    T_CURVE       DZeroCurves[3];

    FXVOLATILITY_DATA FXVolatility;
    CORRELATION_DATA  Correlation;
    FXSMILE_DATA      FXSmile;

    MODELPARAMETERS_DATA  FModelParameters;
    MODELPARAMETERS_DATA  DModelParameters;

} HYB3_MARKET_DATA;

#endif /* CUPSMODL_H */
