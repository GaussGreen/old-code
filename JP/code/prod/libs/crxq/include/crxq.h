/******************************************************************************
 * Module:      CRXQ
 * Submodule:
 * File:        crxq.h
 * Function:    Multi-q mapping and forward measure functionality to enable univariate
 *              and bivariate quasi pricing across credit and rates.
 * Author:      Credit QRD
 * Revision:    $Header: $
 *****************************************************************************/

#ifndef  CRXQ_H
#define  CRXQ_H

#define  INVSQRT2PI             0.398942280401433   /* 1/sqrt(2*pi) */
#define  INF                    1E+10

/* types for memory allocation */
#define  INT                    1
#define  INT_PTR                2
#define  INT_D_PTR              3
#define  LONG                   4
#define  LONG_PTR               5
#define  LONG_D_PTR             6
#define  DOUBLE                 7
#define  DOUBLE_PTR             8
#define  DOUBLE_D_PTR           9
#define  CHAR                   10
#define  CHAR_PTR               11

/* Approximated MQ constants */
#define  CRXQ_NBQ                 34    /* Maximum number of L or R q-values */
#define  CRXQ_Q_SHIFT             1E-6
#define  CRXQ_CUM_NUM_STDEV       10.        /* max allowed arg for NormCum */

/* FA constants */
#define  CRXQ_1D_INT_PTS          500
#define  CRXQ_1D_NUM_STDEV        5
#define  CRXQ_1D_STEP             (2 * (double)CRXQ_1D_NUM_STDEV / ((double)CRXQ_1D_INT_PTS - 1))
#define  CRXQ_BETA_SHIFT          1E-4
#define  CRXQ_FA_BDRY             1E-2  /* cutoff below this * fwdRate */ 
#define  CRXQ_QGUESS_LIMIT        10 

/* parameter ranges (rate, vol, vol of vol, skew) */
#define  CRXQ_MIN_FWD             1E-6
#define  CRXQ_MIN_VOL             1E-8       /* Vol below this is considered to be 0 */
#define  CRXQ_MAX_VOL             2.50
#define  CRXQ_MIN_VOL_CALIB       5E-3
#define  CRXQ_S_STEPS             100
#define  CRXQ_S_DELTA             1E-5
#define  CRXQ_MQ_RESN             1E-14      /* min df/dq resolution in NR  */

/* PA constants */
#define  CRXQ_TAIL_PROB           0.0005
#define  CRXQ_QSOLV_ITER          120         /* maximum number of iterations in q bootstrapping */
#define  CRXQ_QSOLV_TOL           1E6 * DBL_EPSILON    /* tolerance for convergence in q bootstapping */

/* Bivariate corellation limits */
#define CRXQ_MAX_CORR  0.999999
#define CRXQ_MIN_CORR -0.999999

/* Bivariate epsilon limits */
#define CRXQ_MIN_EPS   1e-6

/* macros */

#ifndef  COLLAR
#define  COLLAR(a,b,x) (MIN((b), MAX((a), (x))))
#endif   

#ifndef  SHIFT_ZERO
#define  SHIFT_ZERO(a) ((fabs((a)) < TINY) ? TINY : (a))
#endif   

#ifndef  SHIFT_EPS
#define  SHIFT_EPS(a) ((fabs((a)) < 1e-6) ? 1e-6 : (a))
#endif   

/******************************************************************************
 * CRXQDATA
 * Structure that encapsulates Multi-Q distribution parameters.
 * Computes distribution of forward as a mapping of 
 *     - Gaussian z ~ N(0,1). 
 *     - z -> x = muMQ + sigMQ * sqrt(optExpy) * z 
 *     - x -> H(x) 
 *     - H(x) -> Y = fwdRate * ( volC (H(x) - 1) + fwdC), 
 *     .
 *  where fwdRate is the expected forward rate, K, C are calibrated constants 
 *  so that K( C E[H(x)] - 1) + 1 = 1, so that fwdRate is the expectation.
 *  The other constraint is that the ATM option prices correctly. Given that one
 *  has chosen the q and delta coefficients, there are still two more degrees
 *  of freedom. For credit variables, we will tend to choose muMQ=0 and 
 *  sigmaMQ = sigATM * sqrt(optExpy) and let C and K be calculated. With 
 *  rates variables, the practice has been to choose K=C=1 and to calibrate
 *  muMQ and sigmaMQ so that the forward and option prices match.
 *  H(x) = ki + di * (exp(-qi(x-xi)) - 1)/qi is a sectional q-mapping
 *  where q=0 corresponds to local Normal shape, and q=1 corresponds to
 *  local Log-normal shape. We choose kR0=kL0=1 and dL0=dR0=1.
 ******************************************************************************/
typedef struct {

    /* Keep records of some of the individual input params for convenience   */
    double  optExpy;    /* Time to option expiry in years                   */
    double  sigATM;     /* BS annualized vol for ATM option                 */
    double  optATM;     /* BS ATM option price                              */

    /* Calibrated Multi-Q parameters for left and right sides. */
    /* xL[0]=xR[0] = 0; kR[0]=kL[0] = 1; dR[0]=kR[0]=1; */
    long    nbQL;        /*  Number of qs on the left side                  */
    long    nbQR;        /* Number of qs on the right side                 */
    double  qL[CRXQ_NBQ+1]; /* Left-side q-values.                             */
    double  kL[CRXQ_NBQ+1]; /* Left-side strike boundaries                     */
    double  dL[CRXQ_NBQ+1]; /* Left-side slope-parameters                      */
    double  xL[CRXQ_NBQ+1]; /* Left-side gaussian-space strike boundaries      */
    double  qR[CRXQ_NBQ+1]; /* Right-side q-values.                            */
    double  kR[CRXQ_NBQ+1]; /* Right-side strike boundaries                    */
    double  dR[CRXQ_NBQ+1]; /* Right-side slope-parameters                     */
    double  xR[CRXQ_NBQ+1]; /* Right-side gaussian-space strike boundaries     */

    double  fwdRate;     /* Expected forward rate                           */
    double  sigMQ;       /* Volatility of gaussian driver variable          */
    double  muMQ;        /* Mean of gaussian driver variable                */
    double  K;           /* Affine scaling constant determined by ATM price */
    double  C;           /* 1/E[H(x)] scaling constant                      */
    long    calcFwd;     /* Tell if should calc conditional fwd in bivar  */
} CRXQDATA;

#define CRXQ_CREDIT_SPREAD 0
#define CRXQ_INTEREST_RATE 1

/******************************************************************************
 * CRXQFADATA
 * FA (forward adjusted) parameters. This has been expanded from the Q3 version
 * to allow for both credit spreads and interest rates. The annuity is 
 * calculated differently for each.
 *****************************************************************************/
typedef struct {
    long   rateType; /* CRXQ_CREDIT_SPREAD or CRXQ_INTEREST_RATE              */
    double freqSwap; /* Frequency of the rate/CDS spread                   */

    /* DELAY PARAMETERS */
    double freqDel;  /* Frequency of delay rate                             */
    double matDel;   /* Delay interval in years                             */
    double alphaDel; /* Annuity alpha for delay zero rate                   */
    double powerDel; /* VNFM power for delay calc                           */

    /* ANNUITY RATE CONVEXITY PARAMETERS */
    double freqAnn;  /* Frequency of annuity zero rate/cds spread           */
    double matAnn;   /* Maturity (tenor) of rate/CDS spread                 */
    double alphaAnn; /* Annuity alpha for ann zero rate                     */
    double powerAnn; /* VNFM power for annuity calc                         */

    /* Extras added for credit spreads (so can calculate risky annuity)      */
    double recovery; /* Assumed recovery rate*/
    double ryAnnuity;/* Expected forward risky annuity*/
    double IRate;    /* fixed IR zero rate from exercise to maturity        */
    double accrFactor; /* ratio between accr and no accr fwd par spreads    */

    CRXQDATA  *mq;             /* (O) Underlying annuity rate distribution MQ*/ 
} CRXQFADATA;

/******************************************************************************
 * POSSIBLE INTRUMENT TYPES - USE AS optType IN CRXQPAYOFF
 *****************************************************************************/
#define  CRXQ_ADJ_FWD             0
#define  CRXQ_CALL                1
#define  CRXQ_PUT                 2
#define  CRXQ_IN                  3
#define  CRXQ_OUT                 4

#define  CRXQ_VNL                 0
#define  CRXQ_BIN                 10
#define  CRXQ_ANN                 100   
#define  CRXQ_TEC                 200
#define  CRXQ_VEGA                1000
#define  CRXQ_POW                 1100
#define  CRXQ_QUAD                1110

/* bivariate quasi option types */
#define  CRXQ_CALL_SUM            301
#define  CRXQ_PUT_SUM             302

#define  CRXQ_CALL_PROD           331     
#define  CRXQ_PUT_PROD            332     
#define  CRXQ_CALL_PERC           341     
#define  CRXQ_PUT_PERC            342     
#define  CRXQ_CALL_PERC_WGT       343     
#define  CRXQ_PUT_PERC_WGT        344     

#define  CRXQ_FLR_W_FLR           305
#define  CRXQ_FLR_W_CAP           306
#define  CRXQ_CAP_W_FLR           307
#define  CRXQ_CAP_W_CAP           308

#define  CRXQ_FLR_W_FLR_EMBED     355
#define  CRXQ_FLR_W_CAP_EMBED     356
#define  CRXQ_CAP_W_FLR_EMBED     357
#define  CRXQ_CAP_W_CAP_EMBED     358

#define  CRXQ_IN_BIRIB            313
#define  CRXQ_OUT_BIRIB           314
#define  CRXQ_IN_SPDRIB           323
#define  CRXQ_OUT_SPDRIB          324

#define  CRXQ_IN_BIRIB_EPS        413
#define  CRXQ_OUT_BIRIB_EPS       414
#define  CRXQ_IN_SPDRIB_EPS       423
#define  CRXQ_OUT_SPDRIB_EPS      424

#define  CRXQ_IN_AND_IN           361
#define  CRXQ_IN_OR_IN            362
#define  CRXQ_IN_AND_OUT          363
#define  CRXQ_IN_OR_OUT           364
#define  CRXQ_OUT_AND_IN          365
#define  CRXQ_OUT_OR_IN           366
#define  CRXQ_OUT_AND_OUT         367
#define  CRXQ_OUT_OR_OUT          368

#define  CRXQ_IN_AND_IN_EPS       461
#define  CRXQ_IN_OR_IN_EPS        462
#define  CRXQ_IN_AND_OUT_EPS      463
#define  CRXQ_IN_OR_OUT_EPS       464
#define  CRXQ_OUT_AND_IN_EPS      465
#define  CRXQ_OUT_OR_IN_EPS       466
#define  CRXQ_OUT_AND_OUT_EPS     467
#define  CRXQ_OUT_OR_OUT_EPS      468

#define  CRXQ_COMPLEX_SPD         390

#define  CRXQ_JOINT_FWD           -999

#define  CRXQ_BS_SPRD_CALL        -1 
#define  CRXQ_BS_SPRD_PUT         -2
#define  CRXQ_BS_PERC_CALL        -3
#define  CRXQ_BS_PERC_PUT         -4
#define  CRXQ_BS_SPRD_RATE_CALL   -11
#define  CRXQ_BS_SPRD_RATE_PUT    -12
#define  CRXQ_BS_PERC_RATE_CALL   -13
#define  CRXQ_BS_PERC_RATE_PUT    -14
#define  CRXQ_BS_PERC_SPRD_CALL   -21
#define  CRXQ_BS_PERC_SPRD_PUT    -22

/* basis options calibration type */
#define  CRXQ_BS_CAL_FWD           1
#define  CRXQ_BS_CAL_VOL           2

/* arithmetic on option types */
#define  CRXQ_PAY_TYPE(x)         ((x)/10 * 10)
#define  CRXQ_COP_TYPE(x)         ((x) - CRXQ_PAY_TYPE(x))
#define  CRXQ_COP(cop)            (((cop) == 1) ? CRXQ_CALL : CRXQ_PUT)
#define  CRXQ_COP_TO_COP(CRXQ_cop)  (((CRXQ_cop) == CRXQ_CALL) ? 1. : -1.)

/* instrument description: settlement type */
#define  CRXQ_CASH_SETL           1
#define  CRXQ_PHYS_SETL           0      

/* instrument description: payoff parameters */
#define  CRXQ_MAX_PAY_PARAMS      10
#define  CRXQ_MAX_V_PARAMS         2

/******************************************************************************
 * CRXQPAYOFF
 * Structure representing the type of option - univariate or bivariate
 *****************************************************************************/
typedef struct {
    double  params[CRXQ_MAX_PAY_PARAMS];  /* Scalar parameters              */
    double  *vParam[CRXQ_MAX_V_PARAMS];   /* Vector parameters              */
    long    vParamLen[CRXQ_MAX_V_PARAMS]; /* Vector parameter lengths       */
    CRXQDATA *mq[2];                      /* Measure data (poss. bivariate) */
    double  strike;                       /* Strike                         */
    double  cop;                          /* (1,-1) = (Call,Put)            */
    double  corr;                         /* Correlation (if bivariate)     */
    long    optType;                      /* See instruments above          */
} CRXQPAYOFF;

/******************************************************************************
 * CRXQFPAYOFF
 * Payoff function signature
 *****************************************************************************/
typedef int CRXQFPAYOFF(
    CRXQPAYOFF *prm, 
    double *pt, 
    double smooth,
    double *payoff
    );

/* instrument description: payoff smoothing */
#define  CRXQ_SMOOTH_FACTOR       2.    

/*-----------------------------------------------------------------------------
 * SOURCE FILE: CRXQ.c
 *---------------------------------------------------------------------------*/

/******************************************************************************
 *  CRXQBSQPricer
 *  Prices a vailla call or a put using a single q. 
 *  Returns SUCCESS or FAILURE
 *****************************************************************************/
int CRXQBSQPricer(
    double  Y,                  /* Fwd yield                     */
    double  K,                  /* Strike                        */
    double  T,                  /* Option expiration in years    */
    double  s,                  /* Annualized volatility         */
    double  Q,                  /* Q weight                      */
    long    I,                  /* Instrument)                   */
    double *P);                 /* Price & Vega                  */


/******************************************************************************
 *  CRXQPricer
 *  Prices a vailla call or a put using the multi-q mapping defined by mq. 
 *  You must make sure that mq has been initialized.
 *  Returns SUCCESS or FAILURE
 *****************************************************************************/
int CRXQPricer (
    const CRXQDATA *mq,      /* (I) Initialized Q-distribution parameters    */
    long           optType,  /* (I) Type of option (CRXQ_CALL or CRXQ_PUT only)*/
    double         strike,   /* (I) Option strike price                      */
    double         *price    /* (O) Option price, if returned SUCCESS        */
    );

/******************************************************************************
 *  CRXQBinPricer
 *  Prices a binary call or put 
 *  You must make sure that mq has been initialized.
 *  Returns SUCCESS or FAILURE
 *****************************************************************************/
int CRXQBinPricer(
    const CRXQDATA *mq,      /* (I) Initialized Q-distribution parameters    */
    long           optType,  /* (I) Type of option (CRXQ_CALL or CRXQ_PUT only)*/
    double         strike,   /* (I) Option strike price                      */
    double         *price    /* (O) Option price, if returned SUCCESS        */
    );

/******************************************************************************
 *  CRXQLevPricer
 *  Prices a leveraged call or put: (L*Y-K)+ or (K-LY)+
 *  You must make sure that mq has been initialized.
 *  Returns SUCESS or FAILURE
 *****************************************************************************/
int CRXQLevPricer(
    const CRXQDATA *mq,      /* (I) Initialized Q-distribution parameters    */
    long           optType,  /* (I) Type of option (CRXQ_CALL or CRXQ_PUT only)*/
    double         strike,   /* (I) Option strike price                      */
    double         leverage, /* (I) Leverage for option                      */
    double         *price    /* (O) Option price, if returned SUCCESS        */
    );

/******************************************************************************
 *  CRXQCalibrateCRXQDATA
 *  Calibrate a CRXQDATA multi-q distribution structure from a given set of
 *  inputs. muQ is set to zero, and sigMQ is set to volATM*sqrt(optionExpiry).
 *  The input deltas are mapped to strike boundaries 
 *  k = exp(NormInv(delta) * sigMQ), so that delta=0.5 is ATM.
 *  There must be one fewer delta than there are qs: the middle value is used
 *  twice to ensure that the L and R sides match up.
 *  Returns SUCCESS or FAILURE
 *****************************************************************************/
int CRXQCalibrateCRXQDATA(
    double        optionExpiry,  /* (I) Option time to expiry in years        */
    double        forwardRate,   /* (I) Forward rate                          */
    double        volATM,        /* (I) ATM option annualized BS volatility   */
    int           numberOfQs,    /* (I) Number of Qs: must be an even number  */
    const double* inputQ,        /* (I) Vector of q-values                    */
    const double* inputDelta,    /* (I) Vector of strike deltas               */
    CRXQDATA*      mq            /* (O) CRXQDATA structure to calibrate       */
    );	  

/******************************************************************************
 * CRXQMap
 * Maps the Gaussian driver value, xVal to the MQ-mapped yield value
 * yield = fwdRate (volC (H(xVal) - 1) + fwdC), where H(.) is the q-mapping
 * function defined as H(x) = ki + di(exp(q(x-xi)) - 1)/qi
 *****************************************************************************/
int CRXQMap (
    const CRXQDATA* mq,     /* (I) Q-mapping definiing structure              */
    double         xVal,    /* (I) Gaussian-space x-value                     */
    double*        yield    /* (O) Resulting yield, if returned SUCCESS       */
    );


/******************************************************************************
 * CRXQGridPricer
 * Price a payoff by integrating over a simple multi-q measure
 *****************************************************************************/
int CRXQGridPricer(
    CRXQPAYOFF*  optPayoff, /* (I) Option payoff                            */
    CRXQFPAYOFF* payFunc,   /* (I) Option payoff function                   */
    CRXQDATA*    mq,        /* (I) Multi-q measure                          */
    double*      premium    /* (O) Option price                             */
    );

/******************************************************************************
 * CRXQDens
 * Numerically create a distribution density for a multi-q distribution
 *****************************************************************************/
int CRXQDens(
    const CRXQDATA* mq,         /* (I) Multi-q distribution                 */
    double          start,      /* (I) Start in std. normal space.          */
    double          end,        /* (I) End in std. normal space.            */
    long            numGridPts, /* (I) No. points in grid, yields, dens     */
    double*         grid,       /* (O) Grid of x pts ~ N(muMQ, sigMQ)       */
    double*         yields,     /* (O) Grid of Q-mapped yield points        */
    double*         dens,       /* (O) Probability density function         */
    double*         normC       /* (O) Normalization constant               */
    );    

/******************************************************************************
 * CRXQApproximateFAWithMQ
 * Calibrate a CRXQDATA multi-q measure to approximate a forward-adjusted 
 * measure (used for faster bivariate pricing).
 * Uses the number of qs on each side of the forward specified up to a 
 * maximum CRXQ_NBQ to fit the distribution numerically. The strike points (k)
 * are evenly spaced on each side from the minimum yield to the forward, 
 * and from the forward to the maximum yield in the FA distribution. 
 * The distribution is set-up to match the CDF at the strike points exactly,
 * then the affine adjustment parameters (K,C) are set to ensure that the
 * forward and the ATM call price match. In principle, this messes-up the
 * fit of the CDF but, in practice, the difference is very small if numQ is
 * large enough.
 *****************************************************************************/ 
int CRXQApproximateFAWithMQ(
    const CRXQFADATA* fa, /* (I) Given FA measure that is to be approximated*/   
    int numQ,             /* (I) Number of qs on each side of fwd to fit    */
    CRXQDATA*         pa  /* (O) Resulting multi-q measure                  */    
    );         

/*-----------------------------------------------------------------------------
 * SOURCE FILE: CRXQfapricer.c
 *---------------------------------------------------------------------------*/

/******************************************************************************
 * CRXQFAPricer
 * Price an option by numerically integrating over a forward-adjusted measure
 *****************************************************************************/
int CRXQFAPricer(
    CRXQPAYOFF* optPayoff,
    CRXQFADATA* fa,
    double*     premium
    );

/******************************************************************************
 * CRXQFAGridPricer
 * Computes the change of measure to FA, and then does 1D numerical
 * integration to find the option price. Bivariate options are priced
 * by the inner integral being explicitly evaluated inside the payoff
 * function. Assumes that the multi-q measure is already calibrated. 
 *****************************************************************************/
int CRXQFAGridPricer (
    CRXQPAYOFF*  optPayoff, /* (I) Option payff structure                    */
    CRXQFPAYOFF* payFunc,   /* (I) Option payoff function                    */
    CRXQFADATA*  fa,        /* (I) Measure data                              */
    double*      premium    /* (O) Option price                              */
    );

/******************************************************************************
 * CRXQFACreditInit
 * If the rate is a credit spread rather than an interest-rate then, in order
 * to calculate the annuity for the measure change, you must know the IR
 * zero rate for the spread period, and the recovery rate.
 *****************************************************************************/
int CRXQFACreditInit(
    double recoveryRate,
    double accrualFactor,
    double interestRate,
    double interestRateFreq,
    double riskyAnnuity,  
    CRXQFADATA* fa
    );

/******************************************************************************
 * The annuity-zero and delay-zero rates/spreads are both calculated in the
 * same form for a CDS spread and for an interest rate (both called R here).
 * ZeroRate = alpha * R^p
 * where alpha = Zfwd/(Rfwd)^p * exp(0.5 * sigmaR^T * p * (1-p))
 * where p = correlation*sigmaZ/sigmaR
 * The correlation and the sigmaZ are calculated from a VNFM approximation
 * This function initializes the following variables, which are the same for
 * credit and rates
 *     - freqDel, matDel, alphaDel, powerDel (Delay zero rate)
 *     - freqAnn, matAnn, alphaAnn, powerAnn (Annuity zero)
 *     .
 * Note that you need to call CRXQFACreditInit afterwards if your rate is
 * a credit spread rather than an interest rate.
 *****************************************************************************/
int CRXQFASmileInit(
    long          rateType,      /* (I) CRXQ_CREDIT_SPREAD or CRXQ_INTEREST RATE */
    double        expiry,        /* (I) observation date in yrs             */
    double        sigATM,        /* (I) ATM vol of rate/spread              */
    double        start,         /* (I) swap/zero rate start in yrs         */
    long          freq,          /* (I) comp. freq. of all rates/spreads    */     
    double        swapMat,       /* (I) rate/spread tenor in yrs            */
    double        cdsParSpread,  /* (I) fwd par spread                      */
    double        swapRate,      /* (I) fwd par rate                        */
    double        fwdAnnuity,    /* (I) forward (swap/spread) annuity       */
    double        zeroRateSwap,  /* (I) zero rate for same interval         */
    double        payDelay,      /* (I) payment delay                       */   
    double        zeroRatePay,   /* (I) zero rate for delay interval        */    
    long          numVnfmParams, /* (I) number of VNFM params               */
    double* vnfmParams,          /* (I) VNFM model parameters               */
    long          convexAdjSetl, /* (I) convexity settlement type           */
    long          delayAdjSetl,  /* (I) delay settlement type               */
    CRXQFADATA*    fa            /* (I/O) FA data structure to calibrate    */
    );

/******************************************************************************
 * CRXQVNFMZero2Swap
 * 
 * 
 *****************************************************************************/
int CRXQVNFMZero2Swap(
    double       expiry,
    double       volStart,
    double       rateStart,
    long         freqSwap,
    double       swapMat,
    double       swapRate,
    double       fwdAnnuity,
    double       zeroMat,
    double       zeroRate,
    long         numVnfmParams,
    double       *vnfmParams,
    double       *swapVol,
    double       *zeroVol,
    double       *zeroSwapCorr
    );

/******************************************************************************
 * CRXQFADens
 * Computes the adjusted-forward (or payment) measure probability density
 * explicitly for a given rate (or credit spread) with a given multi-Q
 * distribution. The form of the annuity will depend on whether the rate is
 * an interest rate or a credit spread
 *****************************************************************************/
int CRXQFADens (
    const CRXQFADATA *fa,      /* (I) Rate/spd type, delay and Q-dist info  */
    double         start,      /* (I) Start in std. normal space.           */
    double         end,        /* (I) End in std. normal space.             */
    long           numGridPts, /* (I) No. of points to generate.            */
    double         *grid,      /* (O) Grid of x points ~ N(muMQ, sigMQ)     */
    double         *yields,    /* (O) Grid of Q-mapped yield points         */
    double         *faDens,    /* (O) Probability density function          */
    double         *normC      /* (O) Normalization factor for faDens       */
    );

/*-----------------------------------------------------------------------------
 * SOURCE FILE: CRXQgrid1d.c
 *---------------------------------------------------------------------------*/

/******************************************************************************
 * CRXQSimpsonPricer1D
 * Calculate option price by using Simpson's rule to integrate the payoff
 * over the calibrated adjusted-forward measure.
 *****************************************************************************/
int CRXQSimpsonPricer1D (
    CRXQPAYOFF  *optPayoff, /* (I) Option type and parameters.               */
    CRXQFPAYOFF *payFunc,   /* (I) Option payoff function                    */
    double     *grid,       /* (I) Gaussian grid points                      */
    double     *yield,      /* (I) Q-mapped rate/spread values               */
    double     *dens,       /* (I) FA measure probability densities          */
    double     densNorm,    /* (I) Normalization factor so probs sum to 1    */
    double     *premium     /* (O) Option price                              */
    );    

/* 1d payoffs */                            
CRXQFPAYOFF CRXQPay1D_Yield;
CRXQFPAYOFF CRXQPay1D_Vnl;
CRXQFPAYOFF CRXQPay1D_Bin;
CRXQFPAYOFF CRXQPay1D_Ann;
CRXQFPAYOFF CRXQPay1D_Tec;
CRXQFPAYOFF CRXQPay1D_TecFwd;
CRXQFPAYOFF CRXQPay1D_FixedRibIn;
CRXQFPAYOFF CRXQPay1D_FixedRibOut;
CRXQFPAYOFF CRXQPay1D_Pow;
CRXQFPAYOFF CRXQPay1D_Quad;

/* quasi-1d payoffs */
CRXQFPAYOFF CRXQPay1D_Sum;
CRXQFPAYOFF CRXQPay1D_Prod;
CRXQFPAYOFF CRXQPay1D_Perc;
CRXQFPAYOFF CRXQPay1D_PercWgt;
CRXQFPAYOFF CRXQPay1D_YldYld;
CRXQFPAYOFF CRXQPay1D_YldNull;
CRXQFPAYOFF CRXQPay1D_VnlNull;
CRXQFPAYOFF CRXQPay1D_BinNull;
CRXQFPAYOFF CRXQPay1D_InRIB; 
CRXQFPAYOFF CRXQPay1D_InRIB_EPS; 
CRXQFPAYOFF CRXQPay1D_OutRIB;
CRXQFPAYOFF CRXQPay1D_OutRIB_EPS;
CRXQFPAYOFF CRXQPay1D_InSpdRIB;
CRXQFPAYOFF CRXQPay1D_InSpdRIB_EPS;
CRXQFPAYOFF CRXQPay1D_OutSpdRIB;
CRXQFPAYOFF CRXQPay1D_OutSpdRIB_EPS;
CRXQFPAYOFF CRXQPay1D_FlrFlrOrCapCap;
CRXQFPAYOFF CRXQPay1D_FlrFlrOrCapCapEmbedFlt;
CRXQFPAYOFF CRXQPay1D_FlrCapOrCapFlr;
CRXQFPAYOFF CRXQPay1D_FlrCapOrCapFlrEmbedFlt;
CRXQFPAYOFF CRXQPay1D_InAndIn;
CRXQFPAYOFF CRXQPay1D_InOrIn;
CRXQFPAYOFF CRXQPay1D_InAndOut;
CRXQFPAYOFF CRXQPay1D_InOrOut;
CRXQFPAYOFF CRXQPay1D_OutAndIn;
CRXQFPAYOFF CRXQPay1D_OutOrIn;
CRXQFPAYOFF CRXQPay1D_OutAndOut;
CRXQFPAYOFF CRXQPay1D_OutOrOut;
CRXQFPAYOFF CRXQPay1D_ComplexSPD;


/*----------------------------------------------------------------------------- 
 * Monte Carlo related definitions
 */

#define CRXQ_DEV_DIM_MAX 2

/* Defines for sobseq/sobseq2. */
#define MAXBIT 30
#define MAXDIM 6

/* Defines for ran2. */
#define IM1  2147483563
#define IM2  2147483399
#define AM   (1.0/IM1)
#define IMM1 (IM1-1)
#define IA1  40014
#define IA2  40692
#define IQ1  53668
#define IQ2  52774
#define IR1  12211
#define IR2  3791
#define NTAB 32
#define NDIV (1+IMM1/NTAB)
#define EPS  1.2e-7
#define RNMX (1.0-EPS)

#define CRXQ_BV_RAN2_SEED -9971717

/* Integration types. */
#define CRXQ_BV_RAN2  1
#define CRXQ_BV_SOBOL 2
#define CRXQ_MC_NUM_PTS 100000


/*----------------------------------------------------------------------------- 
 * crxqgrid2d.c
 */
int     Q3MCPricer            (long      dimX,   
                               CRXQDATA  **smX,    
                               double   *rhoX,   
                               CRXQPAYOFF   *pfX,    
                               CRXQFPAYOFF  *pFcn,   
                               long      genTyp, 
                               long      npts,   
                               double   *result);

/* 2d Monte Carlo payoffs */
CRXQFPAYOFF CRXQPay2D_Sum;
CRXQFPAYOFF CRXQPay2D_Prod;
CRXQFPAYOFF CRXQPay2D_Perc;
CRXQFPAYOFF CRXQPay2D_PercWgt;
CRXQFPAYOFF CRXQPay2D_YldYld;
CRXQFPAYOFF CRXQPay2D_YldNull;
CRXQFPAYOFF CRXQPay2D_NullYld;
CRXQFPAYOFF CRXQPay2D_VnlNull;
CRXQFPAYOFF CRXQPay2D_NullVnl;
CRXQFPAYOFF CRXQPay2D_MinMaxIn;
CRXQFPAYOFF CRXQPay2D_MinMaxOut;
CRXQFPAYOFF CRXQPay2D_MinMaxSpdIn;
CRXQFPAYOFF CRXQPay2D_MinMaxSpdOut;
CRXQFPAYOFF CRXQPay2D_FlrFlrOrCapCap;
CRXQFPAYOFF CRXQPay2D_FlrFlrOrCapCapEmbedFlt;
CRXQFPAYOFF CRXQPay2D_FlrCapOrCapFlr;
CRXQFPAYOFF CRXQPay2D_FlrCapOrCapFlrEmbedFlt;
CRXQFPAYOFF CRXQPay2D_InAndIn;
CRXQFPAYOFF CRXQPay2D_InOrIn;
CRXQFPAYOFF CRXQPay2D_InAndOut;
CRXQFPAYOFF CRXQPay2D_InOrOut;
CRXQFPAYOFF CRXQPay2D_OutAndIn;
CRXQFPAYOFF CRXQPay2D_OutOrIn;
CRXQFPAYOFF CRXQPay2D_OutAndOut;
CRXQFPAYOFF CRXQPay2D_OutOrOut;
CRXQFPAYOFF CRXQPay2D_ComplexSPD;

void   *DR_Array              (int          type,
                               int          nl,
                               int          nh);

int     Free_DR_Array         (void         *Array,
                               int          type,
                               int          nl,
                               int          nh);



/*-----------------------------------------------------------------------------
 * mcutil.c
 */

double ran2(long *seed);

void ran2x2(
    long   *seed, /* (I) */ 
    double *X     /* (O)  Pair of uniform deviates.*/
    );

void sobseq(
    int    n,   /* (I) n<0 Initialize routine, n>0 generate n dim point */
    double X[]  /* (O) */
    );

int BoxMuller(
    long    genTyp, /* (I) Specifies generator for uniform deviates. */
    long    dimX,   /* (I) Dimensionality of deviates.               */
    long    numX,   /* (I) Number of deviates of dimension dimX.     */
    double *X       /* (O) Deviates.                                 */
    );

int Gauss(
    long    dimX, /* (I)   Requested dimensionality of deviates.            */
    long    numX, /* (I)   Requested number of deviates of dimension dimX.  */
    double *sigX, /* (I)   Requested std devs for corr. deviates.           */
    double *muX,  /* (I)   Requested means for corr. deviates.              */
    double *rhoX, /* (I)   Requested linear correlation structure.          */
    double *X     /* (I/O) dimX * numX n(0,1) deviates to transform on input*/
                  /*       numX corr. n(mu,sig) deviates of dimension dimX. */
    );


#endif
