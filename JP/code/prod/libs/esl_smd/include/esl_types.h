#ifndef ESL_TYPES_DOT_H
#define ESL_TYPES_DOT_H

#include "esl_macros.h"

#ifdef  __cplusplus
extern "C" {
#endif

/** Coupon payment stub types. Describe how coupon payments are made if they 
    occure on non coupon payment dates.
*/
typedef enum
{
    ESL_STUB_NONE  = 'N',  /**< Receive no coupon when selling. Receive full 
                                coupon on payment date when buying. This is 
                                dirty price of a bond.  
                            */
    ESL_STUB_BOND  = 'B',  /**< Receive partial coupon with accrual period from
                                accrual start to current date when selling. Pay 
                                partial coupon with accrual period from accrual 
                                start to current date, receive full coupon on 
                                payment date. Coupon has to be known in advance.
                                This is clean price of a bond.        
                            */
    ESL_STUB_SIMPLE= 'S',
    ESL_STUB_SWAP  = 'S',  /**< Receive partial coupon with accrual period from
                                current date to accrual end when buying.
                                _SWAP is the usual term for this but some analytics
                                may use _SIMPLE, so the enum is 
                            */
    
    ESL_STUB_PAR   = 'P'   /**< Receive par rate from current date to payment 
                                date and accrued over the remaining period.  
                                This will ensure that the floating leg is priced
                                at par.                       
                            */

} ESL_STUB_TYPE;
 


/** Stub location. Partial coupon payment results from non integral number 
    for coupon periods per tenor. While determining coupon payment dates we
    can start from the maturity date going backward. This will result in 
    shorter first coupon period or the 'front' stub. Alternatively, we start 
    from the inception date and move forward to maturity. This will result in 
    shorter last coupon period or the 'back' stub.  'Long' stub refers to
    a short stub plus one regular coupon interval.
*/
typedef enum
{
    ESL_STUB_SHORT_FRONT = 0,
    ESL_STUB_SHORT_BACK,
    ESL_STUB_LONG_FRONT,
    ESL_STUB_LONG_BACK

} ESL_STUB_LOC;


/** Day count convention
*/
typedef enum
{
    ESL_DCC_30_360 = '3',
    ESL_DCC_ACT_360= '0',
    ESL_DCC_ACT_365= '5',
    ESL_DCC_ACT_ACT= 'A'

} ESL_DCC;


/** Frequency type. Set of well defined frequencies. For other intervals
    one must use DateInterval type.
*/
typedef enum
{
    ESL_FREQ_ANNUAL      = 'A',
    ESL_FREQ_SEMI_ANNUAL = 'S',
    ESL_FREQ_QUARTERLY   = 'Q',
    ESL_FREQ_IMM         = 'I',
    ESL_FREQ_MONTHLY     = 'M',
    ESL_FREQ_WEEKLY      = 'W',
    ESL_FREQ_DAILY       = 'D'

} ESL_FREQ;


/** Money market basis
*/
typedef enum
{
    ESL_MMB_ACT = 0,
    ESL_MMB_360,
    ESL_MMB_365

} ESL_MMB;


/** Option type
*/
typedef enum
{
    ESL_OPT_CALL = 0,
    ESL_OPT_PUT,
    ESL_OPT_CALLABLE,
    ESL_OPT_PUTTABLE,
    ESL_OPT_DIGITAL,
    ESL_OPT_REDEEMER,
    ESL_OPT_STRADDLE,
    ESL_OPT_STRANGLE,
    ESL_OPT_RISKREVERSAL,
    ESL_OPT_COLLAR

} ESL_OPT_TYPE;


/** Knock in/out type.
*/
typedef enum
{
    ESL_KNOCK_IN,
    ESL_KNOCK_OUT

} ESL_KNOCK_TYPE;


/** Knock in/out range type.
*/
typedef enum
{
    ESL_KNOCK_INSIDE,
    ESL_KNOCK_OUTSIDE

} ESL_KNOCK_RANGE;


/** Long or short 
*/
typedef enum
{
    ESL_LONG = 0,
    ESL_SHORT

} ESL_LONGSHORT;

/** Simple or compound rate
*/
typedef enum
{
    SIMPLE_RATE = 'S',
    COMPOUND_RATE = 'C'

} ESL_RATECONV;

/** Reset in advance or arrears
*/
typedef enum
{
    RESET_ADVANCE,
    RESET_ARREARS

} ESL_RESET_TYPE;

/** types for Dev Kit convert usage 
*/
typedef enum
{
    DEVKIT_INT,
    DEVKIT_LONG,
    DEVKIT_DOUBLE,
    DEVKIT_CHAR,
    DEVKIT_TDATE

} DEVKIT_TYPE;

/** Interpolation Type
*/
typedef enum 
{
    LINEAR_INTERP,
    STAIRCASE_INTERP

} ESL_INTERP_TYPE;

/* 
 *  Input structures.
 */

/**    THE TREE SLICE                         */
typedef  double *  TSLICE;      

/** ZERO CURVE DATA STRUCTURE */  
typedef struct _T_CURVE
{
    /* Base date information */
    long    Today;                  /**< Today's date */
    int     SpotDays;               /**< Spot days    */
    long    ValueDate;              /**< Value date   */

    /* Underlying yield curve conventions */
    char    SwapFreq;               /**< Benchmark swap frequency            */
    char    SwapDCC[4];             /**< Benchmark swap day count convention */
    int     MMB;                    /**< Money market basis (360 or 365)     */

    /* Zero curve */
    int     NbZero;                 /**< Number of zeros                   */
    long    ZeroDate[MAXNBDATE];    /**< Zero maturity dates               */
    double  Zero[MAXNBDATE];        /**< Zero rates (Annual ACT/365 basis) */

} T_CURVE;


/** MARKET VOLATILITY CURVE INPUT DATA STRUCTURE */
typedef struct _MKTVOL_DATA
{
    /* Volatility curve */
    long    BaseDate;               /**< Volatility curve base date      */
    int     NbVol;                  /**< Nb of points in vol curve       */
    int     NbCetVol;               /**< Number of vols used in CET      */
    long    VolDate[MAXNBDATE];     /**< Volatility dates                */
    double  Vol[MAXNBDATE];         /**< Vol curve                       */
    int     VolUsed[MAXNBDATE];     /**< TRUE if vol used in calibration */
    int     CetNbIter;              /**< Number of iterations in Cet     */
    double  CetVegaError;           /**< Vega error in Cet               */
    double  Aweight[6][MAXNBDATE];  /**< Output: orthogonal weights      */

    /* Forward swap information */
    char    Freq;                   /**< Frequency of underlying rate    */
    char    DCC;                    /**< Day count convention            */
    long    SwapSt[MAXNBDATE];      /**< Underlying swap start           */
    long    SwapMat[MAXNBDATE];     /**< Underlying swap maturity        */

    /* Model parameters */
    double  QLeft;                  /**< Left Q mapping coeff            */
    double  QRight;                 /**< Right Q mapping coeff           */
    double  FwdShift;               /**< Shift mapping coeff             */
    double  Alpha[3];               /**< Relative size factors           */
    double  Beta[3];                /**< Mean reversions                 */
    double  Rho[3];                 /**< Correlation                     */

    /* Backbone parameters */
    double  Bbq;                    /**< Backbone weight: 1-Logn, 0-Norm */
    double  VolNorm;                /**< Total normal model vol          */
    double  VolLogn;                /**< Total lognormal model vol       */

    /* Calibration flags */
    int     SkipFlag;               /**< Skip calibration failure points */
    int     CalibFlag;              /**< Index calibration flag          */
    int     FilterSpotVolFlag;      /**< TRUE = filter vol, else FALSE   */ 
    char    SmoothingFlag;          /**< Cet smoothing (Y or N)          */
    char    TraceFlag;              /**< Print CET.prn or not            */

	/* Spread smile parameters */
	double  Afac;                   /**< Skew level                      */
	double  Bfac;                   /**< Smile level                     */
	double  Cfac;                   /**< Intensity                       */
	double  Dfac;                   /**< Decay                           */


} MKTVOL_DATA;

/** BASE VOLATILITY DATA STRUCTURE
 *  based on the format of the market data in the various *.dat files */
typedef struct BASEVOL_DATA
{
    char        Frequency;
    int         NbVols;
    long        VolDates[MAXNBDATE];
    double      Vols[MAXNBDATE];

} BASEVOL_DATA;

/** SWAP VOLATILITY DATA STRUCTURE
 *  based on the format of the market data in the various *.dat files */
typedef struct SWAPVOL_DATA
{
    int NbSwaptionExpiries;
    int NbSwapTenors;

    long SwaptionExpiries[NBSWAPTION];  /* in months */
    long SwapTenors[NBSWAPTION];        /* in months */

    double VolMatrix[NBSWAPTION][NBSWAPTION];

} SWAPVOL_DATA;


/** CLAIM BANK STRUCTURE */
typedef struct _CLAIM_BANK
{
	/* element information */
    TSLICE  *Slices;           /**< array of slices in the bank             */
    long    *EvDates;          /**< date the payoff is known for each slice */
    long    *ErDates;          /**< slice is obsolete for all t < this date */

	/* bank information */
    int      TotNbSlices;
    int      NbActiveSlices;   /**< nb of active slices                     */

    /* internal parameters */
    TSLICE   auxSlice;         /**< extra slice for interp etc              */
    long     MaxErDate;        /**< max earliest use date                   */
    int      Locked;           /**< if TRUE, disable all activities but Add */

} CLAIM_BANK;


/** alib date structure */
typedef long int TDATE;

/** EVENT LIST DATA STRUCTURE */
typedef struct _EVENT_LIST
{
    /** Length of Dates and Curves             */
    int     NbEntries;             
    /** Nb of curves used out of MAXNBEVCURVES */
    int     NbCurves;              
    /** Event dates (e.g. knock out, exercise) */
    long    *Dates;                
    /** Curves characterising event            */
    double  *Curve[MAXNBEVCURVES]; 

} EVENT_LIST;


/** CRITICAL DATE DATA STRUCTURE */  
typedef struct _CRIT_DATE
{
    /** Critical date type  */
    int     Type;        

    /** Critical date       */
    long    CritDate;    
    
    /** Supplementary dates */
    long    SuppDate[3]; 
    
    /** Associated values   */
    double  Value[5];    

} CRIT_DATE;

/** 
*   Output structures.
*/

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
    
    int            prob_calc[ESL_NB_EVENT];      /**< Flag setting if the specific event has been populated */
    PROB_OUT_DATA  prob_out_data[ESL_NB_EVENT];  /**< Struct stroring statistic *
                                                    *  0 is for the Call proba  *
                                                    *  1 is for the KO proba    */

} OPT_OUT_DATA;


/** CET (Calibration enhancement tool) STRUCTURES */
typedef struct _CET_OUT_DATA         /* CET  output with current mkt vols    */
{                  
    double   Annuity[MAX_INST];      /**< Discounted annuities for the b'marks */
    double   ParYield[MAX_INST];     /**< Par yields for the benchmarks        */

    double   BSPrice[MAX_INST];      /**< B & Scholes prices for the benchmarks*/
    double   BSVega[MAX_INST];       /**< B&S vegas for the benchmarks         */
    double   TreePrice[MAX_INST];    /**< Tree prices for the benchmarks       */
    double   TreePricePrev[MAX_INST];/**< Tree prices for previous iteration   */

    double   PriceDiff[MAX_INST];    /**< Differences between B&S and tree     */
    double   PriceDiffInVega[MAX_INST];

    double   VolPrev[MAX_INST];          /**< Previous iteration vol               */
    double   FilterVol[MAX_INST];    /**< Previous iteration filtered vol      */

    double   FOVega[MAX_INST];       /**< Price/Vol approximation of vega      */

} CET_OUT_DATA;
/*  NB: The size of this array is consistent  with the limitaiton imposed *
  * by the critical date structure. The number of chosen benchmark inst's *
  * must be less than MAX_INST.                                           */


/** IR MC simulation structure for state variable */
typedef struct
{
    int    NbPaths;                      /**< Nb of paths                      */

    long   BaseDate;                     /**< Simulation base date             */
    long   CurrDate;                     /**< Current date                     */

    /* Current OU states of a swap yield */
    double X[MAXNBSTATES];               /**< Current OU states                */
    long   SwapSt;                       /**< Current swap start               */
    long   SwapMat;                      /**< Current swap maturity            */
    char   Freq;                         /**< Current swap frequency           */
    char   DCC;                          /**< Current swap DCC                 */

    long   Seed;                         /**< Random number seed               */

} IR_SIM;




#ifdef  __cplusplus
}
#endif


#endif

