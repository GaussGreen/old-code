/******************************************************************************
 * Module:      Q3
 * Submodule:
 * File: q3.h   
 * Function:    
 * Author:      Interest Rates DR
 * Revision:    $Header$
 *****************************************************************************/
#ifndef  Q3_H
#define  Q3_H


/* general constants */
#define  TRUE                   1
#define  FALSE                  0
#define  SUCCESS                0
#define  FAILURE                -1
#define  TINY                   1E-10
#define  INVSQRT2PI             0.398942280401433   /* 1/sqrt(2*pi) */

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

/* SV constants */
#define  Q3_SV_NUM_STDEV        5.
#define  Q3_SV_INT_PTS          25
#define  Q3_NUM_SOLVER_ITER     4

/* MQ constants */
#define  Q3_NBQ                 34
#define  Q3_Q_SHIFT             1E-6
#define  Q3_MQ_NORM_TAIL        3.
#define  Q3_CUM_NUM_STDEV       10.        /* max allowed arg for NormCum */
#define  Q3_MQ_RESN             1E-14      /* min df/dq resolution in NR  */

/* FA constants */
#define  Q3_1D_INT_PTS          500
#define  Q3_1D_NUM_STDEV        5
#define  Q3_1D_STEP             (2 * (double)Q3_1D_NUM_STDEV / ((double)Q3_1D_INT_PTS - 1))
#define  Q3_BETA_SHIFT          1E-4
#define  Q3_FA_BDRY             1E-2   

/* SV parameter ranges (rate, vol, vol of vol, skew) */
#define  Q3_MIN_FWD             1E-6
#define  Q3_MIN_VOL             1E-5
#define  Q3_MAX_VOL             5.00
#define  Q3_MAX_VOL_VOL         2.00
#define  Q3_MIN_SKEW            -5.
#define  Q3_MAX_SKEW            7.

/* MQ parameter ranges (deltas and taus) */
#define  Q3_MIN_DELTA_L         0.01
#define  Q3_MAX_DELTA_L         0.45
#define  Q3_MIN_DELTA_R         0.55
#define  Q3_MAX_DELTA_R         0.99
#define  Q3_MIN_TAU_L           0.
#define  Q3_MAX_TAU_R           10.

/* MQ Calibration constants for Newton-Raphson */
#define  Q3_MQ_FIX              1
#define  Q3_MQ_TOL              2
#define  Q3_MQ_FIX_SAFE         3
#define  Q3_MQ_TOL_SAFE         4

#define  Q3_Q_STEPS             20
#define  Q3_OTM_TOL             1E-7

#define  Q3_M_STEPS             20
#define  Q3_M_DELTA             1E-5
#define  Q3_FWD_TOL             1E-5

#define  Q3_S_STEPS             100
#define  Q3_S_DELTA             1E-5
#define  Q3_ATM_TOL             1E-4
#define  Q3_S_MAX_STEP          2E-1   

#define  Q3_MIN_VOL_CALIB       5E-3

/* PA constants */
#define  Q3_NUM_BIN_L1          7
#define  Q3_NUM_BIN_L2          7
#define  Q3_NUM_BIN_R1          7
#define  Q3_NUM_BIN_R2          7
#define  Q3_NUM_BIN_L           Q3_NUM_BIN_L1 + Q3_NUM_BIN_L2
#define  Q3_NUM_BIN_R           Q3_NUM_BIN_R1 + Q3_NUM_BIN_R2
#define  Q3_EQUIK_ITVL          0
#define  Q3_EQUIX_ITVL          1
#define  Q3_ITVL_L1             Q3_EQUIX_ITVL
#define  Q3_ITVL_L2             Q3_EQUIK_ITVL
#define  Q3_ITVL_R1             Q3_EQUIX_ITVL
#define  Q3_ITVL_R2             Q3_EQUIX_ITVL
#define  Q3_TAIL_PROB           0.0005
#define  Q3_PROB_STEP           .1
#define  Q3_TAIL_ITER           10
#define  Q3_TAIL_Q_L            0
#define  Q3_TAIL_Q_R            0
#define  Q3_QGUESS_LIMIT        10 
#define  Q3_QSOLV_TOL           1E6 * DBL_EPSILON
#define  Q3_XSOLV_TOL           1E3 * DBL_EPSILON
#define  Q3_QSOLV_ITER          120
#define  Q3_X_DELTA              .0001
#define  Q3_Q_XOVER              1e-7 
#define  Q3_VOLMXITER            10 
#define  Q3_VOLMXSTEP            0.1
#define  Q3_VOLDELTA             0.0001
#define  Q3_PREQERR              0.00000001
#define  Q3_MAX_ZERO_ADJ         1.0
#define  Q3_LGST_MU              1.0
#define  Q3_LGST_SIG             0.05
#define  Q3_LGST_MINX            0.0
#define  Q3_LGST_MAXX            10.0
#if (Q3_NUM_BIN_R + Q3_NUM_BIN_L + 6) > Q3_NBQ  
#   error Too many additional gaussian nodes.
#endif

/* Bivariate corellation limits */
#define Q3_MAX_CORR  0.999999
#define Q3_MIN_CORR -0.999999

/* Monte Carlo Stuff */
#define Q3_DEV_DIM_MAX 2

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

#define Q3_BV_RAN2_SEED -9971717

/* Integration types. */
#define Q3_BV_RAN2  1
#define Q3_BV_SOBOL 2
#define Q3_MC_NUM_PTS 100000

/* macros */


#ifndef  MAX
#define  MAX(a,b) ((a) > (b) ? (a) : (b))
#endif

#ifndef  MIN
#define  MIN(a,b) ((a) < (b) ? (a) : (b))
#endif

#ifndef  COLLAR
#define  COLLAR(a,b,x) (MIN((b), MAX((a), (x))))
#endif   COLLAR

#ifndef  SHIFT_ZERO
#define  SHIFT_ZERO(a) ((fabs((a)) < TINY) ? TINY : (a))
#endif   SHIFT_ZERO

#ifndef  SHIFT_EPS
#define  SHIFT_EPS(a) ((fabs((a)) < 1e-10) ? 1e-10 : (a))
#endif   SHIFT_EPS

/*t----------------------------------------------------------------------------
 * Smile measure data structure  
 */

/* SV (stochastic vol) parameters */
typedef struct {
    double  expiry;                 /* (I) Expiration time                  */
    double  fwdRate;                /* (I) Forward rate                     */
    double  sigATM;                 /* (I) Mkt ATM volatility               */
    double  q;                      /* (I) Mkt skew                         */
    double  volVolMkt;              /* (I) Mkt scaled vol of vol            */
    double  bbR;                    /* (I) Mkt rate backbone                */
    double  bbV;                    /* (I) Mkt vol backbone                 */
    double  sigSV;                  /* (O) Mdl internal vol                 */
    double  volVolSV;               /* (O) Mdl internal vol of vol          */
} SVDATA;

/* MQ (multi q) parameters  */
typedef struct {
    /* Market inputs */
    double  expiry;                 /* (I) Expiration time                  */
    double  fwdRate;                /* (I) Forward rate                     */
    double  sigATM;                 /* (I) Mkt ATM volatility               */
    double  tauL;                   /* (I) Mkt tau left of fwd              */
    double  tauR;                   /* (I) Mkt tau right of fwd             */
    double  tauN;                   /* (I) Mkt tau for negative rates       */
    double  callL[Q3_NBQ];          /* (I) Target call prices left of fwd   */
    double  putR[Q3_NBQ];           /* (I) Target put prices right of fwd   */
    double  optATM;                 /* (I) Target ATM option price          */

    /* Internal MultiQ parameters */
    long    nbQL;                   /* (I) Mdl number of left intervals     */
    long    nbQR;                   /* (I) Mdl number of right interval     */
    double  kL[Q3_NBQ];             /* (I) Mdl strike/fwd left of fwd       */
    double  kR[Q3_NBQ];             /* (I) Mdl strike/fwd right of fwd      */
    double  qL[Q3_NBQ];             /* (O) Mdl q left of fwd                */
    double  qR[Q3_NBQ];             /* (O) Mdl q right of fwd               */
    double  muMQ;                   /* (O) Mdl mean of normal dist          */
    double  sigMQ;                  /* (O) Mdl vol of normal dist           */

    /* Affine adjustment */
    long    affineAdj;              /* (I) Use affine adjustment            */
    double  fwdC;                   /* (O) Mdl forward adjustment           */
    double  volC;                   /* (O) Mdl ATM vol adjustment           */

    /* Calibration setup (encoded in NCK) */
    long    calibType;              /* (I) Calibration type (FIX or TOL)    */
    int     calibSafe;              /* (I) Apply safe mu and sig steps or no*/
    long    mSteps;                 /* (I) Number of NR steps for muMQ      */
    long    sSteps;                 /* (I) Number of NR steps for sigMQ     */
    long    qSteps;                 /* (I) Number of NR steps for q's       */
    double  mDelta;                 /* (I) Step for muMQ num derivative     */
    double  sDelta;                 /* (I) Step fog sigMQ num derivative    */    
    double  fwdTol;                 /* (I) Calib tolerance for forward      */
    double  atmTol;                 /* (I) Calib tolerance for ATM vol      */  
    double  otmTol;                 /* (I) Calib tolerance for OTM vol      */

    /* bounds for positivity of call spreads */
    double  muOverSigMax, muOverSigMin;

    /* Other selectable behavior */
    long    calibQs;                /* (I) Calibrate Q values               */
    long    calcFwd;                /* (I) Calc. fwd, ignore struct value   */       
} MQDATA;

/* FA (forward adjusted) parameters */
typedef struct {
    double  freqAnn;                /* (I) Frequency of swap rate           */
    double  freqDel;                /* (I) Frequency of delay rate          */
    double  matAnn;                 /* (I) Maturity (tenor) of swap rate    */
    double  matDel;                 /* (I) Delay interval                   */
    double  alphaAnn;               /* (O) Annuity alpha for ann zero rate  */
    double  alphaDel;               /* (O) Annuity alpha for delay zero rate*/
    double  powerAnn;               /* (O) VNFM power for annuity calc      */
    double  powerDel;               /* (O) VNFM power for delay calc        */
    MQDATA  *mq;                    /* (O) Underlying annuity measure MQ    */ 
} FADATA;

/*e*/


/* instrument description: option types */
#define  Q3_ADJ_FWD             0
#define  Q3_CALL                1
#define  Q3_PUT                 2
#define  Q3_IN                  3
#define  Q3_OUT                 4

#define  Q3_VNL                 0
#define  Q3_BIN                 10
#define  Q3_ANN                 100   
#define  Q3_TEC                 200
#define  Q3_VEGA                1000
#define  Q3_POW                 1100

/* bivariate option types */
#define  Q3_CALL_SUM            301
#define  Q3_PUT_SUM             302

#define  Q3_CALL_PROD           331     
#define  Q3_PUT_PROD            332     
#define  Q3_CALL_PERC           341     
#define  Q3_PUT_PERC            342     

#define  Q3_FLR_W_FLR           305
#define  Q3_FLR_W_CAP           306
#define  Q3_CAP_W_FLR           307
#define  Q3_CAP_W_CAP           308

#define  Q3_FLR_W_FLR_EMBED     355
#define  Q3_FLR_W_CAP_EMBED     356
#define  Q3_CAP_W_FLR_EMBED     357
#define  Q3_CAP_W_CAP_EMBED     358

#define  Q3_IN_BIRIB            313
#define  Q3_OUT_BIRIB           314
#define  Q3_IN_SPDRIB           323
#define  Q3_OUT_SPDRIB          324

#define  Q3_IN_BIRIB_EPS        413
#define  Q3_OUT_BIRIB_EPS       414
#define  Q3_IN_SPDRIB_EPS       423
#define  Q3_OUT_SPDRIB_EPS      424

#define  Q3_JOINT_FWD           -999

#define  Q3_BS_SPRD_CALL        -1 
#define  Q3_BS_SPRD_PUT         -2
#define  Q3_BS_PERC_CALL        -3
#define  Q3_BS_PERC_PUT         -4
#define  Q3_BS_SPRD_RATE_CALL   -11
#define  Q3_BS_SPRD_RATE_PUT    -12
#define  Q3_BS_PERC_RATE_CALL   -13
#define  Q3_BS_PERC_RATE_PUT    -14
#define  Q3_BS_PERC_SPRD_CALL   -21
#define  Q3_BS_PERC_SPRD_PUT    -22

/* basis options calibration type */
#define Q3_BS_CAL_FWD           1
#define Q3_BS_CAL_VOL           2

/* arithmetic on option types */
#define  Q3_PAY_TYPE(x)         ((x)/10 * 10)
#define  Q3_COP_TYPE(x)         ((x) - Q3_PAY_TYPE(x))
#define  Q3_COP(cop)            (((cop) == 1) ? Q3_CALL : Q3_PUT)
#define  Q3_COP_TO_COP(q3_cop)  (((q3_cop) == Q3_CALL) ? 1. : -1.)

/* instrument description: settlement type */
#define  Q3_CASH_SETL           1
#define  Q3_PHYS_SETL           0      

/* instrument description: payoff parameters */
#define  Q3_MAX_PAY_PARAMS      10
#define  Q3_MAX_V_PARAMS         2

/*t----------------------------------------------------------------------------
 * Option type data structure 
 */
typedef struct {
    double params[Q3_MAX_PAY_PARAMS];  /* Scalar parameters           */
    double *vParam[Q3_MAX_V_PARAMS];   /* Vector parameters           */
    long   vParamLen[Q3_MAX_V_PARAMS]; /* Vector parameter lengths    */
    MQDATA *mq[2];                     /* Measure data (bivariate)    */
    double strike;                     /* Strike                      */
    double cop;                        /* (1,-1) = (Call,Put)         */
    double corr;                       /* Correlation  (bivariate)    */
    long   optType;                    /* See instruments above       */
} PAYOFF;

/*t----------------------------------------------------------------------------
 * Payoff function signature. 
 */
typedef int FPAYOFF(
    PAYOFF *prm, 
    double *pt, 
    double smooth,
    double *payoff
    );

/*e*/

/* instrument description: payoff smoothing */
#define  Q3_SMOOTH_FACTOR       2.    


/*-----------------------------------------------------------------------------
 * function prototypes
 */

/*-----------------------------------------------------------------------------
 * multiq.c
 */
int     Q3MQPricer            (MQDATA       *mq,
                               long         type,
                               double       strike,
                               double       *price);

int     Q3MQLevPricer         (MQDATA       *mq,
                               long         type,
                               double       strike,
                               double       leverage,
                               double       *price);

int     Q3MQBinPricer         (MQDATA       *mq,    
                               long         optType, 
                               double       strike,  
                               double       *price);

int     Q3MQMap               (MQDATA       *mq,
                               double       xval,
                               double       *yield);

int     Q3MQSmileInit         (double       fwdRate,
                               double       sigATM,
                               double       expiry,
                               double       *smile,
                               MQDATA       *mq);

int     Q3MQSmileInitFromMQ   (MQDATA       *mq1,
                               MQDATA       *mq2);                                        

void    Q3MQCopyMQ            (MQDATA       *mq1,
                               MQDATA       *mq2);                                        


int     Q3MQTargetSV          (SVDATA       *sv,
                               double       *mqGuess,
                               MQDATA       *mq);

int     Q3MQTargetFA          (FADATA       *fa,
                               MQDATA       *pa);  

int     Q3MQCalib             (MQDATA       *mq);

int     Q3MQCalibMean         (MQDATA       *mq);

int     Q3MQCalibMidQs        (MQDATA       *mq);

int     Q3MQCalibSpread       (double       fwd,
                               double       vol,
                               double       expiry,
                               double      *smile,
                               MQDATA      *mq);

int     Q3DecodeNCK           (double       nck,
                               MQDATA       *mq);

int     Q3SVToMQ              (double       fwdRate,
                               double       sigATM,
                               double       expiry,
                               double       *smile,
                               MQDATA       *mq);

int     Q3MQBootstrapQ        (MQDATA       *mq);
 
int     Q3MQSolveMap4Q        (double       A,   
                               double       B,
                               double       *q);

/*----------------------------------------------------------------------------- 
 * stochvol.c
 */
int     Q3SVPricer            (SVDATA       *sv,
                               long         type,
                               double       strike,
                               double       *price);

int     Q3SVSmileInit         (double       fwdRate,
                               double       sigATM,
                               double       expiry,
                               double       *smile,
                               SVDATA       *sv);

int     Q3SVCalib             (SVDATA       *sv);

int     Q3SVBinPricer         (SVDATA       *sv,    
                               double       strike, 
                               double       *price);

/*-----------------------------------------------------------------------------
 * fapricer.c
 */
int     Q3FAPricer            (PAYOFF       *optPayoff,
                               FADATA       *fa,
                               double       *premium);

int     Q3FAGridPricer        (PAYOFF       *optPayoff, 
                               FPAYOFF      *payFunc,   
                               FADATA       *fa,        
                               double       *premium);

int     Q3FATailStrikes       (FADATA       *fa,
                               double       *tailProb,   
                               double       *tailStrike);

int     Q3FASmileInit         (double       expiry,              
                               double       sigATM,
                               double       start,              
                               long         freq,                
                               double       swapMat,           
                               double       swapRate,           
                               double       fwdAnnuity,        
                               double       zeroRateSwap,       
                               double       payDelay,           
                               double       zeroRatePay,            
                               long         numVnfmParams,    
                               double       *vnfmParams,    
                               long         cashPhysSetl,     
                               FADATA       *fa);    

int     Q3FADens              (FADATA       *fa,
                               double       start,
                               double       end,
                               long         numGridPts,
                               double       *grid,
                               double       *yields,
                               double       *faDens,
                               double       *normC);

int     Q3VNFMZero2Swap       (double       expiry,
                               double       volStart,
                               double       rateStart,
                               long         rateFreq,
                               double       swapMat,
                               double       swapRate,
                               double       fwdAnnuity,
                               double       zeroMat,
                               double       zeroRate,
                               long         numVnfmParams,
                               double       *vnfmParams,
                               double       *swapVol,
                               double       *zeroVol,
                               double       *zeroSwapCorr);


/*-----------------------------------------------------------------------------
 * mqpricer.c
 */
int     Q3MQGridPricer        (PAYOFF       *optPayoff,
                               FPAYOFF      *payFunc,
                               MQDATA       *mq,
                               double       *premium);

int     Q3MQDens              (MQDATA       *mq,
                               double       start,
                               double       end,
                               long         numGridPts,
                               double       *grid,
                               double       *yields,
                               double       *dens,
                               double       *normC);    

/*-----------------------------------------------------------------------------
 * grid1d.c
 */
int     Q3SimpsonPricer1D     (PAYOFF       *optPayoff, 
                               FPAYOFF      *payFunc,
                               double       *grid,  
                               double       *yield, 
                               double       *dens,  
                               double       densNorm, 
                               double       *premium);    

int     Q3Payoff1D            (PAYOFF       *optPayoff,
                               double       yield,
                               double       smoothStep,
                               double       *payoff);

/* 1d payoffs */                            
FPAYOFF Q3Pay1D_Yield;
FPAYOFF Q3Pay1D_Vnl;
FPAYOFF Q3Pay1D_Bin;
FPAYOFF Q3Pay1D_Ann;
FPAYOFF Q3Pay1D_Tec;
FPAYOFF Q3Pay1D_TecFwd;
FPAYOFF Q3Pay1D_FixedRibIn;
FPAYOFF Q3Pay1D_FixedRibOut;
FPAYOFF Q3Pay1D_Pow;

/* quasi-1d payoffs */
FPAYOFF Q3Pay1D_Sum;
FPAYOFF Q3Pay1D_Prod;
FPAYOFF Q3Pay1D_Perc;
FPAYOFF Q3Pay1D_YldYld;
FPAYOFF Q3Pay1D_YldNull;
FPAYOFF Q3Pay1D_VnlNull;
FPAYOFF Q3Pay1D_MinMaxIn;
FPAYOFF Q3Pay1D_MinMaxIn_PayInd; 
FPAYOFF Q3Pay1D_MinMaxOut;
FPAYOFF Q3Pay1D_MinMaxOut_PayInd;
FPAYOFF Q3Pay1D_MinMaxSpdIn;
FPAYOFF Q3Pay1D_MinMaxSpdOut;
FPAYOFF Q3Pay1D_MinMaxSpdIn_EPS;
FPAYOFF Q3Pay1D_MinMaxSpdOut_EPS;
FPAYOFF Q3Pay1D_FlrFlrOrCapCap;
FPAYOFF Q3Pay1D_FlrFlrOrCapCapEmbedFlt;
FPAYOFF Q3Pay1D_FlrCapOrCapFlr;
FPAYOFF Q3Pay1D_FlrCapOrCapFlrEmbedFlt;

/*----------------------------------------------------------------------------- 
 * grid2d.c
 */
int     Q3MCPricer            (long      dimX,   
                               MQDATA  **smX,    
                               double   *rhoX,   
                               PAYOFF   *pfX,    
                               FPAYOFF  *pFcn,   
                               long      genTyp, 
                               long      npts,   
                               double   *result);

/* 2d Monte Carlo payoffs */
FPAYOFF Q3Pay2D_Sum;
FPAYOFF Q3Pay2D_Prod;
FPAYOFF Q3Pay2D_Perc;
FPAYOFF Q3Pay2D_YldYld;
FPAYOFF Q3Pay2D_YldNull;
FPAYOFF Q3Pay2D_NullYld;
FPAYOFF Q3Pay2D_VnlNull;
FPAYOFF Q3Pay2D_NullVnl;
FPAYOFF Q3Pay2D_MinMaxIn;
FPAYOFF Q3Pay2D_MinMaxOut;
FPAYOFF Q3Pay2D_MinMaxSpdIn;
FPAYOFF Q3Pay2D_MinMaxSpdOut;
FPAYOFF Q3Pay2D_FlrFlrOrCapCap;
FPAYOFF Q3Pay2D_FlrFlrOrCapCapEmbedFlt;
FPAYOFF Q3Pay2D_FlrCapOrCapFlr;
FPAYOFF Q3Pay2D_FlrCapOrCapFlrEmbedFlt;


/*----------------------------------------------------------------------------- 
 * alloc.c
 */
void   *DR_Array              (int          type,
                               int          nl,
                               int          nh);

int     Free_DR_Array         (void         *Array,
                               int          type,
                               int          nl,
                               int          nh);

/*----------------------------------------------------------------------------- 
 * bas.c
 */
double  NormDens              (double       x);

double  NormCum               (double       x);

double  NormCumInv            (double       prob); 

double  BSQInt                (double       a,
                               double       b,
                               double       x1,
                               double       x2);

int     BSQPricer             (double       yield,
                               double       strike,
                               double       expiry,
                               double       sigATM,
                               double       Q,
                               long         optType,
                               double       *price);

int     BSQImpVol             (double       yield,
                               double       strike,
                               double       expiry,
                               double       Q,
                               double       price,
                               long         optType,
                               double       volGuess,
                               double       *impVol);

double  ExpCErrFcn            (double       a,
                               double       b,
                               double       x);


/*----------------------------------------------------------------------------- 
 * util.c
 */

double  Q3SimpsIntegral       (double       step,
                               long         numPoints,
                               double       *fValues);

double  Q3SmoothStepFcn       (double       x,
                               double       step);

double  Q3SmoothMAX           (double       x,
                               double       step);

double  Q3ProbInBarriers      (double x,      /* (I) input value            */
                               double loBar,  /* (I) low barrier            */
                               double hiBar,  /* (I) high barrier           */
                               double loEps,  /* (I) low barrier smoothing  */
                               double hiEps,  /* (I) high barrier smoothing */
                               double smooth);/* (I) smoothing step size    */

double  Q3ProbOutBarriers     (double x,      /* (I) input value            */
                               double loBar,  /* (I) low barrier            */
                               double hiBar,  /* (I) high barrier           */
                               double loEps,  /* (I) low barrier smoothing  */
                               double hiEps,  /* (I) high barrier smoothing */
                               double smooth);/* (I) smoothing step size    */

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

/*----------------------------------------------------------------------------- 
 * error.c
 */

void    Q3ErrMsg             (char *fmt, ...);
int     Q3ErrCallbackSet     (int (*cbfunc)(const char *string));


#endif


