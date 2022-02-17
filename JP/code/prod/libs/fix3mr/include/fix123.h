/****************************************************************************/
/*      Declaration file for structure of data.                             */
/****************************************************************************/
/*      FIX123.h                                                            */
/****************************************************************************/

/*
$Header: /nasdev/export2/home/drdev/cvsadmin/cvs/libs/fix3/include/fix123.h,v 1.27 2004/07/27 15:38:51 hjorge Exp $
*/

/* Use safe multiple inclusion for Jerry Cohen's stream libraries */

#ifndef _fix123_h
#define _fix123_h



#define  TRUE          1
#define  FALSE         0
#define  SUCCESS       0
#define  FAILURE       -1
#define  MAXBUFF       250
#define  PI            3.141592653
#define  BIG           exp (15.)
#define  TINY          1E-10
#define  ERROR         1E-7           /* Error tolerance                */
#define  BARRIER_TOL   1E-4           /* Error tolerance for barriers   */
#define  JUMPCOEFF     3.0            /* Jump size coefficient          */
#define  NBCRITDATE    50             /* Number of critical date arrays */
#define  ZbkEVENT      (NBCRITDATE-3) /* the first zerobank event type  */

#define  MAXNBDATE     500    /* Max nb of elements in input date array */
#define  MAXNBSTATES   200    /* Max nb of state variables              */
#define  MAXNBEVCURVES 5      /* Max nb of event curves in EVENT_LIST   */
#define  QCUTOFF       1E-4   /* Normal model for |q|<QCUTOFF           */

#define   MAX_ITERATIONS    500
#define   MAX_INST          NBCRITDATE

/* different mode choices for GetDLOffset function in date.c */
#define  CbkEXACT   0
#define  CbkLOWER   -1
#define  CbkHIGHER  1


/* Types for memory allocation */
#define     CRITDATE        0
#define     INT             1
#define     INT_PTR         2
#define     INT_D_PTR       3
#define     LONG            4
#define     LONG_PTR        5
#define     LONG_D_PTR      6
#define     DOUBLE          7
#define     DOUBLE_PTR      8
#define     DOUBLE_D_PTR    9
#define     CHAR            10
#define     CHAR_PTR        11

/* Rate types */
#define CASH_RATE 'C'
#define SWAP_RATE 'S'


/* Macros */
#ifndef NEAR_INT
#define NEAR_INT(a) (int) ((a) + ((a) >= 0. ? 0.5 : -0.5));
#endif

#ifndef MAX
#define MAX(a,b) ((a) > (b) ? (a) : (b))
#endif

#ifndef MIN
#define MIN(a,b) ((a) < (b) ? (a) : (b))
#endif

#ifndef ABS
#define ABS(a) ((a) > 0 ? (a) : -(a))
#endif

#ifndef SIGN
#define SIGN(a) ((a) >= 0. ? 1 : -1)
#endif

#ifndef SQUARE
#define SQUARE(a) ((a)*(a))
#endif

#ifndef CUBE
#define CUBE(a) ((a)*(a)*(a))
#endif

#ifndef COLLAR
#define COLLAR(amt,cap,flr) \
                MAX(MIN((amt),(cap)),(flr))
#endif


#ifndef ACC_FN
#define ACC_FN(rate,dcf,isSIMPLE)            \
              ( (isSIMPLE) ?                 \
                ((rate)*(dcf)) :             \
                (pow((1+(rate)),(dcf))-1.0)  \
              )
#endif



#ifndef ACC_FN_DIFF  /* equals ACC_FN(r1,dcf,type) - ACC_FN(r2,dcf,type) */
#define ACC_FN_DIFF(r1,r2,dcf,isSIMPLE)                    \
              ( (isSIMPLE) ?                               \
                (((r1)-(r2))*(dcf)) :                      \
                (pow((1+(r1)),(dcf))-pow((1+(r2)),(dcf)))  \
              )
#endif


#ifndef IS_EQUAL
#define IS_EQUAL(m,n) (fabs((m)-(n)) < TINY)
#endif

#ifndef DR_REALLOC
#define DR_REALLOC(dl,s) (((dl)==NULL) ? malloc((s)) : realloc((dl),(s)))
#endif

#ifndef IS_Q
#define IS_Q(q) (fabs((q)) > QCUTOFF)
#endif

#ifndef SIZEOFARRAY
#define SIZEOFARRAY(a)	(sizeof(a)/sizeof(a[0]))
#endif


/* 
 *	Input structures.
 */


/* ZERO CURVE DATA STRUCTURE */  
typedef struct _T_CURVE
{
    /* Base date information */
    long    Today;                  /* Today's date */
    int     SpotDays;               /* Spot days    */
    long    ValueDate;              /* Value date   */

    /* Underlying yield curve conventions */
    char    SwapFreq;               /* Benchmark swap frequency            */
    char    SwapDCC[4];             /* Benchmark swap day count convention */
    int     MMB;                    /* Money market basis (360 or 365)     */

    /* Zero curve */
    int     NbZero;                 /* Number of zeros                   */
    long    ZeroDate[MAXNBDATE];    /* Zero maturity dates               */
    double  Zero[MAXNBDATE];        /* Zero rates (Annual ACT/365 basis) */

} T_CURVE;



/* MARKET VOLATILITY CURVE INPUT DATA STRUCTURE */
typedef struct _MKTVOL_DATA
{
    /* Volatility curve */
    long    BaseDate;               /* Volatility curve base date      */
    int     NbVol;                  /* Nb of points in vol curve       */
    int     NbCetVol;               /* Number of vols used in CET      */
    long    VolDate[MAXNBDATE];     /* Volatility dates                */
    double  Vol[MAXNBDATE];         /* Vol curve                       */
    int     VolUsed[MAXNBDATE];     /* TRUE if vol used in calibration */
    int     CetNbIter;              /* Number of iterations in Cet     */
    double  CetVolShift;            /* Vol shift for initial CET input */
    double  Aweight[6][MAXNBDATE];  /* Output: orthogonal weights      */

    /* Forward swap information */
    char    Freq;                   /* Frequency of underlying rate    */
    char    DCC;                    /* Day count convention            */
    long    SwapSt[MAXNBDATE];      /* Underlying swap start           */
    long    SwapMat[MAXNBDATE];     /* Underlying swap maturity for Idx[0] */
    long    SwapMat2[MAXNBDATE];    /* Underlying swap maturity for Idx[1] */

    /* Model parameters */
    double  QLeft;                  /* Left Q mapping coeff            */
    double  QRight;                 /* Right Q mapping coeff           */
    double  FwdShift;               /* Shift mapping coeff             */
    double  Alpha[3];               /* Relative size factors           */
    double  Beta[3];                /* Mean reversions                 */
    double  Rho[3];                 /* Correlation                     */

    /* Backbone parameters */
    double  Bbq;                    /* Backbone weight: 1-Logn, 0-Norm */
    double  VolNorm;                /* Total normal model vol          */
    double  VolLogn;                /* Total lognormal model vol       */

    /* Calibration flags */
    int     SkipFlag;               /* Skip calibration failure points */
    int     CalibFlag;              /* Index calibration flag          */
    int     FilterSpotVolFlag;      /* TRUE = filter vol, else FALSE   */ 
    char    SmoothingFlag;          /* Cet smoothing (Y or N)          */
    char    TraceFlag;              /* Print CET.prn or not            */

    /* Initial vol guess for CET */
	int     NbCetVolShift;
    long    CetVolShiftExp[MAXNBDATE];
	long    CetVolShiftDate[MAXNBDATE];
    double  CetVolShiftVal[MAXNBDATE];
    double  CetVolShiftArray[MAXNBDATE];
	double  TargetVol[MAXNBDATE];

    /* Time-dependent mr */
    int     NbMr;                   /* number of mr dates              */
	long    MrExp[MAXNBDATE];
    long    MrDate[MAXNBDATE];      /* mean-reversion dates            */
    double  MrInput[MAXNBDATE];     /* mean-reversion term structure   */ 
    double  MrVNFM;                 /* input VNFM mean-reversion       */
    double  MrGuess;                /* mr guess                        */
    double  CetMrShift;             /* mr shift from VNFM for cet      */
    int     MrCalibFlag;            /* calibration flag                */
	int     NbMrCet;
    
	double  BetaTD[3][MAXNBDATE];   /* time-dependent mean reversion   */
    double  BetaBmk[3][MAXNBDATE];  /* mean-reversion on intvals 
	                                   between benchmarks              */  
    int     NbBmkMr;                /* number of mr on bmk dates and beyond */
    long    BmkDate[MAXNBDATE];
    
    long    RatioExp[MAXNBDATE];    /* expiries at which vol ratios are matched */
	long    RatioDate[MAXNBDATE];   /* dates for which vol ratios are 
                                       matched  */  
    double  VNFMRatio[MAXNBDATE];   /* corresponding VNFM ratios */
    int     RatioFlag[MAXNBDATE];   /* TRUE if vol ratio needs to be 
                                       matched at current bmark */
    int     LastRatioIdx;
    double  VolRatio[MAXNBDATE];
    double  VolRatioInput[MAXNBDATE];

} MKTVOL_DATA;

/* New from ESL */
/** Interpolation Type
*/
typedef enum 
{
    LINEAR_INTERP,
    STAIRCASE_INTERP

} ESL_INTERP_TYPE;



/* 
*	Output structures.
*/

  
/* CRITICAL DATE DATA STRUCTURE */  
typedef struct _CRIT_DATE
{
    int     Type;        /* Critical date type  */
    long    CritDate;    /* Critical date       */
    long    SuppDate[3]; /* Supplementary dates */
    double  Value[5];    /* Associated values   */

} CRIT_DATE;



/* EVENT LIST DATA STRUCTURE */
typedef struct _EVENT_LIST
{
    int     NbEntries;             /* Length of Dates and Curves             */
    int     NbCurves;              /* Nb of curves used out of MAXNBEVCURVES */
    long    *Dates;                /* Event dates (e.g. knock out, exercise) */
    double  *Curve[MAXNBEVCURVES]; /* Curves characterising event            */

} EVENT_LIST;



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

    /* Time-dependent mr */
    long    MrDate[MAXNBDATE];             /* mean-reversion dates            */
    double  *BetaTD[3];                    /* time-dependent mean reversion   */
    int     NbMrInt;                       /* Number of MR intervals          */
    long    MrSwitchT[MAX_INST + 10];     /* Array of TP when MR changes     */
    double  *BetaInt[3];                   /* MR on each MR intval            */

} TREE_DATA;



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

} DEV_DATA;



/* OUTPUT DATA STRUCTURE */
#define OPT_OUT_DATA_SIZE   128

typedef struct _PROB_OUT_DATA
{
    double       TotalEventProb;                 /**< Total probability */
    double       ExpEventTime;                   /**< Expected time to event */
    double       StdEventTime;                   /**< Time to event standard deviation */
    double       Fugit;                          /**< Fugit */
    int          Count;                          /**< size of the event schedule */

    long         EventDate[OPT_OUT_DATA_SIZE];   /**< Events dates */
    double       EventProb[OPT_OUT_DATA_SIZE];   /**< Events probabilities */

} PROB_OUT_DATA;

typedef enum
{
    EXER_EVENT,
    KO_EVENT,
    ESL_NB_EVENT
} ESL_EVENT;

typedef struct _OPT_OUT_DATA
{
    double         Option;                       /**< Output price            */
    double         Price  [10];                  /**< Additional outputs - backward compatibility */
    
    PROB_OUT_DATA  prob_out_data[ESL_NB_EVENT];  /**< struct stroring statistic *
							  *  0 is for the Call proba
							  *  1 is for the KO proba    */

} OPT_OUT_DATA;


/* CLAIM BANK STRUCTURE */
typedef struct _CLAIM_BANK
{
	/* element information */
    double **Slices;           /* array of slices in the bank             */
    long    *EvDates;          /* date the payoff is known for each slice */
    long    *ErDates;          /* slice is obsolete for all t < this date */

	/* bank information */
    int      TotNbSlices;
    int      NbActiveSlices;   /* nb of active slices                     */

    /* internal parameters */
    double  *auxSlice;         /* extra slice for interp etc              */
    long     MaxErDate;        /* max earliest use date                   */
    int      Locked;           /* if TRUE, disable all activities but Add */

} CLAIM_BANK;




/* CET (Calibration enhancement tool) STRUCTURES */
typedef struct _CET_OUT_DATA         /* CET  output with current mkt vols    */
{                  

    double   Annuity[MAX_INST];      /* Discounted annuities for the b'marks */
    double   Annuity2[MAX_INST];     /* Discounted annuities for the b'marks */
    double   ParYield[MAX_INST];     /* Par yields for the benchmarks        */
    double   ParYield2[MAX_INST];    /* Par yields for the benchmarks (2nd idx) */

    double   BSPrice[MAX_INST];      /* B & Scholes prices for the benchmarks*/
    double   BSVega[MAX_INST];       /* B&S vegas for the benchmarks         */
    double   TreePrice[MAX_INST];    /* Tree prices for the benchmarks       */
    double   TreePricePrev[MAX_INST];/* Tree prices for previous iteration   */
    double   TreePrice2[MAX_INST];    /* Tree prices for the benchmarks: Idx[1]       */
    double   TreePricePrev2[MAX_INST];/* Tree prices for previous iteration: Idx[1]   */

    double   PriceDiff[MAX_INST];    /* Differences between B&S and tree     */
    double   PriceDiffInVega[MAX_INST];

    double   Vol[MAX_INST];          /* Previous iteration vol               */
    double   FilterVol[MAX_INST];    /* Previous iteration filtered vol      */

    double   FOVega[MAX_INST];       /* Price/Vol approximation of vega      */

} CET_OUT_DATA;
/* NB: The size of this array is consistent  with the limitaiton imposed */
/* by the critical date structure. The number of chosen benchmark inst's */
/* must be less than MAX_INST.                                           */


extern int  ZeroInterpTypeFlag;    /* 0=Linear Zero Cpn; 1=Flat Fwd */


/* IR MC simulation structure for state variable */
typedef struct
{
    int    NbPaths;                      /* Nb of paths                      */

    long   BaseDate;                     /* Simulation base date             */
    long   CurrDate;                     /* Current date                     */

    /* Current OU states of a swap yield */
    double X[MAXNBSTATES];               /* Current OU states                */
    long   SwapSt;                       /* Current swap start               */
    long   SwapMat;                      /* Current swap maturity            */
    char   Freq;                         /* Current swap frequency           */
    char   DCC;                          /* Current swap DCC                 */

    long   Seed;                         /* Random number seed               */
} IR_SIM;


#endif /* _fix123_h */
