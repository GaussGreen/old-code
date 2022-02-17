/****************************************************************************/
/*      Declaration file for structure of data.                             */
/****************************************************************************/
/*      TMX123.h                                                            */
/****************************************************************************/

/*
$Header$
*/

#ifndef _bmx123_h
#define _bmx123_h



/* just for now, as we don't have the multifactor case so there are
 * loads of unreferenced variables */
#pragma warning(disable: 4101)

#define  TRUE          1
#define  FALSE         0
#define  SUCCESS       0
#define  FAILURE       -1
#define  MAXBUFF       250
#define  PI            3.141592653
#define  INVSQRT2PI    0.398942280401433   
#define  BIG           exp (15.)
#define  TINY          1E-10
#define  ERROR         1E-7           /* Error tolerance                */
#define  BARRIER_TOL   1E-4           /* Error tolerance for barriers   */
#define  JUMPCOEFF     3.0            /* Jump size coefficient          */
#define  NBCRITDATE    50             /* Number of critical date arrays */
#define  ZbkEVENT      (NBCRITDATE-3) /* the first zerobank event type  */
#define  NMREVENT      (NBCRITDATE-4) /* Numeraire date event type      */

#define  MAXNBDATE     500    /* Max nb of elements in input date array */
#define  MAXNBSTATES   200    /* Max nb of state variables              */
#define  MAXNBEVCURVES 5      /* Max nb of event curves in EVENT_LIST   */
#define  QCUTOFF       1E-4   /* Normal model for |q|<QCUTOFF           */
#define  NMR_CUTOFF    1E-16  /* Low numeraire cutoff (> 0)             */

#define   MAX_ITERATIONS    10
#define   MAX_INST          NBCRITDATE
#define   SIGMACORR         0.001
#define   CORR_ERROR        0.01
#define   DELTA_CORR        0.05
#define   NB_CORREL_MAX     20
#define   NUM_MAP_DIR1      1.0
#define   NUM_MAP_DIR2      0.10

 




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
#define SIZEOFARRAY(a)  (sizeof(a)/sizeof(a[0]))
#endif


/* 
 *  Input structures.
 */


/* ZERO CURVE DATA STRUCTURE */  
typedef struct                                                                
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


#define    MAXNBEXPIRY   30
#define    MAXNBFWDMAT   20
#define    NBVOLPARS     9          /* sig,q,vvol,bbV,bbR,dL,tl,dR,tR */
#define    NBSTRIKE      5

/* MARKET VOLATILITY SMILE MATRIX INPUT DATA STRUCTURE */
typedef struct
{
    /* Volatility smile matrix coordinates */
    long    VolDate[MAXNBDATE];     /* Basevol: expiration dates       */
    long    Expiry[MAXNBDATE];      /* Swapvol: expiration month       */
    long    FwdMat[MAXNBDATE];      /* Swapvol: tenor month            */

    long    NbExpiry;               /* Swapvol: number of expiries     */
    long    NbFwdMat;               /* Swapvol: number of tenors       */

    /* Volatility smile matrix data */
    double  Vol[NBVOLPARS][MAXNBEXPIRY][MAXNBFWDMAT]; /* vol and smile */

    /* Volatility smile market attributes */
    char    Freq;                   /* Rate frequency in vol quotes    */

} MKTSMILE_DATA;

#define  AA_FWD  0
#define  AA_ATM  1
#define  AA_VOL  2
#define  AA_VAR  3
#define  AA_ALL  4
#define  OU_FWD  0
#define  OU_ATM  1
#define  OU_NMR  2
#define  OU_VAR  3
#define  OU_ALL  4

/* TREE VOLATILITY AND SMILE CURVE DATA STRUCTURE */
typedef struct
{
    /* Volatility and smile curve */
    long    TermDate;               /* Terminal date of TMX            */
    long    BaseDate;               /* Volatility curve base date      */
    int     NbVol;                  /* Nb of points in vol curve       */
    long    VolDate[MAXNBDATE];     /* Volatility dates                */
    double  Vol[NBVOLPARS][MAXNBDATE]; /* Vol curve                    */
    int     VolUsed[MAXNBDATE];     /* TRUE if vol used in calibration */
    long    SwapSt[MAXNBDATE];      /* Underlying swap start           */
    long    SwapMat[MAXNBDATE];     /* Underlying swap maturity        */
    int     SwapMatMos[MAXNBDATE];  /* Swap maturity in months         */
    char    Freq;                   /* Frequency of underlying rate    */
    char    DCC;                    /* Day count convention            */
    double  Aweight[6][MAXNBDATE];  /* Output: orthogonal weights      */

    /* Volatility info at numeraire dates */
    int     NmrStatFlag;            /* Nmr statistics True/False       */
    int     NbNmr;                  /* Nb numeraire dates              */
    long    NmrDate[MAXNBDATE];     /* Numeraire vol dates             */
    long    NmrSwapSt[MAXNBDATE];   /* Numeraire swap start date       */
    long    NmrSwapMat[MAXNBDATE];  /* Numeraire swap maturity date    */
    double  NmrVol[MAXNBDATE][NBVOLPARS]; /* Numeraire volatilities    */
    double  Strike[MAXNBDATE][NBSTRIKE]; /* Strike for Nmr stat        */
    double  VanlPr[MAXNBDATE][NBSTRIKE]; /* Vanila price for Nmr stat  */
    double  TreePr[MAXNBDATE][NBSTRIKE]; /* Tree price for Nmr stat    */
    double  AASta[MAXNBDATE][5];     /* Tree price stat with AA factor */
    double  OUSta[MAXNBDATE][5];     /* Tree price stat                */
    double  TreeSta[MAXNBDATE][6];   /* Other tree price stat          */

    /* Zero curve information as of ValueDate */
    double  ParYield0[MAXNBDATE];   /* Forward swap yield              */
    double  AtmPrice0[MAXNBDATE];   /* Forward ATM price               */
    double  Annuity0 [MAXNBDATE];   /* Annuity as of ValueDate         */
    double  StZero0  [MAXNBDATE];   /* Zero to nmr date                */
    double  TermZero0;              /* Zero to terminal date           */

    /* Model parameters */
    double  Alpha[3];               /* Relative size factors           */
    double  Beta[3];                /* Mean reversions                 */
    double  Rho[3];                 /* Correlation                     */
    double  OWVol[NBVOLPARS];       /* Overwrite vol smile             */
    double  SmileFact;              /* Smile impact on spot vol        */

    /* Tail definitions and numerical configuration */
    double  NbSigmaMQ;              /* Number of std, MQ normal cutoff */
    double  NckMQ;                  /* MultiQ configuration key        */
    double  DeltaNMQ;               /* MultiQ negative delta           */
    double  TauNMQ;                 /* MultiQ negative tau             */
    double  NbSigmaBinL;            /* Nb std for bin map, low strikes */
    double  NbSigmaBinR;            /* Nb std for bin map, hi  strikes */

    /* Backbone parameters */
    double  Bbq;                    /* Backbone weight: 1-Logn, 0-Norm */
    double  VolNorm;                /* Total normal model vol          */
    double  VolLogn;                /* Total lognormal model vol       */

    /* Calibration flags */
    int     SkipFlag;               /* Skip calibration failure points */
    int     CalibFlag;              /* Index calibration flag          */
    int     FilterSpotVolFlag;      /* TRUE = filter vol, else FALSE   */ 
    char    SmoothingFlag;          /* Nmr_Calc smoothing (Y or N)     */

    int     Trace;                  /* Detail level for debug info     */

} MKTVOL_DATA;



/* 
*   Output structures.
*/

  
/* CRITICAL DATE DATA STRUCTURE */  
typedef struct                                                                  
{
    int     Type;        /* Critical date type  */
    long    CritDate;    /* Critical date       */
    long    SuppDate[3]; /* Supplementary dates */
    double  Value[5];    /* Associated values   */

} CRIT_DATE;


/* EVENT LIST DATA STRUCTURE */
typedef struct                                                                  
{
    int     NbEntries;             /* Length of Dates and Curves             */
    int     NbCurves;              /* Nb of curves used out of MAXNBEVCURVES */
    long    *Dates;                /* Event dates (e.g. knock out, exercise) */
    double  *Curve[MAXNBEVCURVES]; /* Curves characterising event            */

} EVENT_LIST;



/* TREE DATA STRUCTURE */
typedef struct                                                                  
{
    /* Critical dates */
    CRIT_DATE   *CritDate[NBCRITDATE];     /* Critical dates description     */
    char        CritType[NBCRITDATE];      /* Type of critical dates         */
    int         *TPtype[NBCRITDATE];       /* Critical type of current TP    */
    int         NbZeros[NBCRITDATE];       /* Number of zeros in zero bank   */
    long        LastProdDate;              /* Last product date              */

    /* Express DEV dates (numeraire + product critical dates) */ 
    int         NbEDevDates;               /* Nb of dates for express DEV    */               
    long        *EDevDate;                 /* Dates where express DEV req    */
    double      **EDevStPrice;             /* Corresp state prices           */

    /* Numeraire (benchmark) dates */
    int         NbNmr;                     /* Nb of numeraire dates          */
    long        *NmrDate;                  /* Numeraire dates                */
    double      **NmrInv;                  /* Inverse numeraire prices       */


    /* Annuities & last zeroes for numeraire interpolation */
    double      **Ann1,  **Ann0;           /* Consecutive calib annuities    */
    double      **Zero1, **Zero0;          /* Consective last zeroes         */
    int         NmrInterpOn;               /* Numeraire interpolation active */

    /* Yield mapping function */
    int         *YLMin;                    /* Lowest level (index) for Y     */
    int         *YLMax;                    /* Highest level (index) for Y    */
    double      **Yield;                   /* Mapped yield at NMR dates      */
    double      **YCumP;                   /* Cum prob level in TN measure   */
    double      **ZAInv;                   /* Lower bound for yield on level */

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
    double      *TermZero[3];              /* Zero to terminal point         */

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

    /* New arguments added by BM */

    /* Discretized Gaussian values for mapping */
    double      **XValues;                 /* First Gaussian slice values               */
    double      **YValues;                 /* Second Gaussian slice values              */

    double      Map_dir1[2];               /* Direction vector for first yield mapping  */
    double      Map_dir2[2];               /* Direction vector for second yield mapping */  
    double      NmrMap_dir[2];             /* Numeraire mapping direction               */

    long        *MappSize;                 /* Mapping size for each Nmr date            */
    
    
    long        NbCorrel;                      /* Number of input correlations      */
    double      Biv_Expiry[NB_CORREL_MAX];     /* Expiry in yrs                     */
    double      Biv_Correl1[NB_CORREL_MAX];    /* Correlation array 1               */
    double      Biv_Correl2[NB_CORREL_MAX];    /* Correlation array 2               */



} TREE_DATA;



/* DEV DATA STRUCTURE */
typedef struct                                                                  
{

    /* Node branching */
    int     *Shift1;
    int     *Shift2;
    int     *Shift3;

    /* Currency or numeraire denominated quantities flags */
    int     NmrToCcy;           /* Numeraire used at current TP or not  */
    int     CcyToNmr;           /* Numeraire used at previous TP or not */

    /* Inverse numeraire; include rescaled values for basis curves  */
    double  *NmrInv[3];         /* Inverse numeraire at current TP      */ 
    double  *NmrInvLag[3];      /* Inverse numeraire at previous TP     */

    /* Probabilities */
    double  *pu, *p0, *pd;      /* 1D */
    double  *quu, *qu0, *qud;   /* 2D */
    double  *q0u, *q00, *q0d;
    double  *qdu, *qd0, *qdd;
    double  *ru, *r0, *rd;      /* 3D */

    double  *NewPrice;          /* Auxiliary slice */

} DEV_DATA;



/* OUTPUT DATA STRUCTURE */
typedef struct                                                                  
{
    double  Option;     /* Output price          */
    double  Price[10];  /* Additional outputs    */

} OPT_OUT_DATA;


/* CLAIM BANK STRUCTURE */
typedef struct                                                                  
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


extern int ZeroInterpTypeFlag; /* 0=Linear Zero Cpn; 1=Flat Fwd */


#endif /* _tmx123_h */
