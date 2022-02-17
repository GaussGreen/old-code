/****************************************************************************/
/*      Declaration file for structure of data.                             */
/****************************************************************************/
/*      TMX123.h                                                            */
/****************************************************************************/

#ifndef _tmx123_h
#define _tmx123_h

/* just for now, as we don't have the multifactor case so there are
 * loads of unreferenced variables */
//#pragma warning(disable: 4101)

#define  TRUE              1
#define  FALSE             0
#define  SUCCESS           0
#define  FAILURE          -1
#define  MAXBUFF           250
#define  PI                3.141592653
#define  INVSQRT2PI        0.398942280401433   
#define  BIG               exp (15.)
#define  TINY              1E-10
#define  ERROR             1E-7           /* Error tolerance                */
#define  BARRIER_TOL       1E-4           /* Error tolerance for barriers   */
#define  JUMPCOEFF         3.0            /* Jump size coefficient          */
#define  NBCRITDATE        50             /* Number of critical date arrays */
#define  ZbkEVENT         (NBCRITDATE-3) /* the first zerobank event type  */
#define  NMREVENT         (NBCRITDATE-4) /* Numeraire date event type      */

#define  MAXNBDATE         500    /* Max nb of elements in input date array */
#define  MAXNBSTATES       200    /* Max nb of state variables              */
#define  MAXNBEVCURVES     5      /* Max nb of event curves in EVENT_LIST   */
#define  QCUTOFF           1E-4   /* Normal model for |q|<QCUTOFF           */
#define  NMR_CUTOFF        1E-16  /* Low numeraire cutoff (> 0)             */

#define  MAX_ITERATIONS    50
#define  MAX_INST          NBCRITDATE
#define  DEFAULT_CETNBITER 20
#define  Nb_Daily_Pts      0              /* Number of daily points         */



/* different mode choices for GetDLOffset function in date.c */
#define  CbkEXACT          0
#define  CbkLOWER         -1
#define  CbkHIGHER         1

/* Types for memory allocation */
#define  CRITDATE          0
#define  INT               1
#define  INT_PTR           2
#define  INT_D_PTR         3
#define  LONG              4
#define  LONG_PTR          5
#define  LONG_D_PTR        6
#define  DOUBLE            7
#define  DOUBLE_PTR        8
#define  DOUBLE_D_PTR      9
#define  CHAR              10
#define  CHAR_PTR          11
#define  MIN_SMILE_EXP_MTHS 6

/* Rate types */
#define  CASH_RATE        'C'
#define  SWAP_RATE        'S'


/* Macros */
#ifndef  NEAR_INT
#define  NEAR_INT(a) (int) ((a) + ((a) >= 0. ? 0.5 : -0.5));
#endif

#ifndef  MAX
#define  MAX(a,b) ((a) > (b) ? (a) : (b))
#endif

#ifndef  MIN
#define  MIN(a,b) ((a) < (b) ? (a) : (b))
#endif

#ifndef  ABS
#define  ABS(a) ((a) > 0 ? (a) : -(a))
#endif

#ifndef  SIGN
#define  SIGN(a) ((a) >= 0. ? 1 : -1)
#endif

#ifndef  SQUARE
#define  SQUARE(a) ((a)*(a))
#endif

#ifndef  CUBE
#define  CUBE(a) ((a)*(a)*(a))
#endif

#ifndef  COLLAR
#define  COLLAR(amt,cap,flr) \
                MAX(MIN((amt),(cap)),(flr))
#endif


#ifndef  ACC_FN
#define  ACC_FN(rate,dcf,isSIMPLE)            \
              ( (isSIMPLE) ?                 \
                ((rate)*(dcf)) :             \
                (pow((1+(rate)),(dcf))-1.0)  \
              )
#endif


#ifndef  ACC_FN_DIFF  /* equals ACC_FN(r1,dcf,type) - ACC_FN(r2,dcf,type) */
#define  ACC_FN_DIFF(r1,r2,dcf,isSIMPLE)                    \
               ( (isSIMPLE) ?                               \
                 (((r1)-(r2))*(dcf)) :                      \
                 (pow((1+(r1)),(dcf))-pow((1+(r2)),(dcf)))  \
               )
#endif


#ifndef  IS_EQUAL
#define  IS_EQUAL(m,n) (fabs((m)-(n)) < TINY)
#endif

#ifndef  DR_REALLOC
#define  DR_REALLOC(dl,s) (((dl)==NULL) ? malloc((s)) : realloc((dl),(s)))
#endif

#ifndef  IS_Q
#define  IS_Q(q) (fabs((q)) > QCUTOFF)
#endif

#ifndef  SIZEOFARRAY
#define  SIZEOFARRAY(a)  (sizeof(a)/sizeof(a[0]))
#endif


/** Interpolation Type
*/
typedef enum 
{
    LINEAR_INTERP,
    STAIRCASE_INTERP

} INTERP_TYPE;


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


#define  MAXNBEXPIRY       30
#define  MAXNBFWDMAT       20
#define  NBVOLPARS         9          /* sig,q,vvol,bbV,bbR,dL,tl,dR,tR */
#define  NBSTRIKE          5
#define  MIDSTRIKE         2

/* Delta points for Smile definition */
extern const double DSTRIKE[];


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

#define  NBSTATS           5
#define  AA_FWD            0
#define  AA_ATM            1
#define  AA_VOL            2
#define  AA_VAR            3
#define  AA_ALL            4
#define  OU_FWD            0
#define  OU_ATM            1
#define  OU_NMR            2
#define  OU_VAR            3
#define  OU_ALL            4

/* TREE VOLATILITY AND SMILE CURVE DATA STRUCTURE */
typedef struct
{
    /* Volatility and smile curve */
    long    TermDate;                        /* Terminal date of TMX             */
    long    BaseDate;                        /* Volatility curve base date       */
    int     NbVol;                           /* Nb of points in vol curve        */
    long    VolDate[MAXNBDATE];              /* Volatility dates                 */
    double  Vol[NBVOLPARS][MAXNBDATE];       /* Vol curve                        */
    int     VolUsed[MAXNBDATE];              /* TRUE if vol used in calibration  */
    double  ParYield0[MAXNBDATE];            /* Forward swap yield               */
    double  Annuity0 [MAXNBDATE];             /* Annuity as of ValueDate         */
    
    double  SwapStrk[NBSTRIKE][MAXNBDATE];   /* Swap Strikes                     */
    double  SwapVnVol[NBSTRIKE][MAXNBDATE];  /* Swap Vanilla Vols                */
    double  SwapTreeVol[NBSTRIKE][MAXNBDATE];/* Swap Vanila Vols                 */    
    long    SwapSt[MAXNBDATE];               /* Underlying swap start            */
    long    SwapMat[MAXNBDATE];              /* Underlying swap maturity         */
    int     SwapMatMos[MAXNBDATE];           /* Swap maturity in months          */
    char    Freq;                            /* Frequency of underlying rate     */
    char    DCC;                             /* Day count convention             */
    double  Aweight[6][MAXNBDATE];           /* Output: orthogonal weights       */
    int     CetNbIter;                       /* Number of iterations in Cet      */
    int     SmlLiqDate[MAXNBDATE];           /* if TRUE, used in smile fitting   */
    
    /* Volatility info at numeraire dates */
    int     NbNmr;                               /* Nb numeraire dates               */
    long    NmrDate[MAXNBDATE];                  /* Numeraire vol dates              */
    double  NmrLibVol[NBVOLPARS][MAXNBDATE];     /* Libor SV-MQ                      */
    
    
    /* Zero curve information as of ValueDate */
    double  AtmLiborPrice0[MAXNBDATE];/* Forward Libor ATM price       */
    double  TermZero0;              /* Zero to terminal date           */
    double  ParLibor0[MAXNBDATE];   /* Forward libor                   */
    double  AnnuityLibor0[MAXNBDATE];/* Libor Annuity as of ValueDate  */

    /* Model parameters */
    double  Alpha[3];               /* Relative size factors           */
    double  Beta[3];                /* Mean reversions                 */
    double  Rho[3];                 /* Correlation                     */
    double  OWVol[NBVOLPARS];       /* Overwrite vol smile             */

    /* MultiQ numerical configuration */
    double  NckMQ;                  /* MultiQ configuration key        */
    double  NbSigmaMQ;              /* Normal Cut-off                  */
    
    /* Backbone parameters */
    double  Bbq;                    /* Backbone weight: 1-Logn, 0-Norm */
 
    /* Calibration flags   */
    int     CalibSmileFlag;         /* TRUE if sml calibration required*/
    int     SkipFlag;               /* Skip Vol Flag                   */
    char    SmoothingFlag;          /* Cet smoothing (Y or N)          */

    /* Control NMR output  */
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

    /* Numeraire variables */
    int         NbNmr;                     /* Nb of numeraire dates          */
    double      **NmrInv;                  /* Inverse numeraire prices       */

    /* Annuities & last zeroes for numeraire interpolation */
    double      *LastZero;                 /* Zero to Nxt Critical date      */
    int         NmrInterpOn;               /* Numeraire interpolation active */
    long        NxtCritDate;               /* Next Critical date             */

    /* Yield mapping function */
    double      **Libor;                   /* Mapped libor at NMR dates      */
    
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
