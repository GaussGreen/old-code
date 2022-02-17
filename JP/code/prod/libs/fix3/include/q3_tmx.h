/******************************************************************************
 * Module:      Q3TMX
 * Submodule:
 * File: q3.h   
 * Function:    
 * Author:      Interest Rates DR
 * Revision:    $Header$
 *****************************************************************************/
#ifndef  Q3TMX_H
#define  Q3TMX_H


/* general constants */
#define  TRUE                   1
#define  FALSE                  0
#define  SUCCESS                0
#define  FAILURE                -1
#define  Q3TMX_TINY                1E-14
#define  Q3TMX_INVSQRT2PI          0.398942280401433   /* 1/sqrt(2*pi) */

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
#define  Q3TMX_SV_NUM_STDEV        5.
#define  Q3TMX_SV_INT_PTS          25
#define  Q3TMX_NUM_SOLVER_ITER     4

/* MQ constants */
#define  Q3TMX_NBQ                 34
#define  Q3TMX_Q_SHIFT             1E-6
#define  Q3TMX_MQ_NORM_TAIL        3.
#define  Q3TMX_CUM_NUM_STDEV       10.        /* max allowed arg for NormCum */
#define  Q3TMX_MQ_RESN             1E-14      /* min df/dq resolution in NR  */
#define  Q3TMX_6Q_DELTAN_DEFAULT   0.
#define  Q3TMX_6Q_TAUN_DEFAULT     0.25

/* FA constants */
#define  Q3TMX_1D_INT_PTS          500
#define  Q3TMX_1D_NUM_STDEV        5
#define  Q3TMX_1D_STEP             (2 * (double)Q3TMX_1D_NUM_STDEV / ((double)Q3TMX_1D_INT_PTS - 1))
#define  Q3TMX_BETA_SHIFT          1E-4
#define  Q3TMX_FA_BDRY             1E-2   

/* SV parameter ranges (rate, vol, vol of vol, skew) */
#define  Q3TMX_MIN_FWD             1E-6
#define  Q3TMX_MIN_VOL             1E-5
#define  Q3TMX_MAX_VOL             3.00
#define  Q3TMX_MAX_VOL_VOL         2.00
#define  Q3TMX_MAX_SKEW            1.65

/* MQ parameter ranges (deltas and taus) */
#define  Q3TMX_MIN_DELTA_L         0.01
#define  Q3TMX_MAX_DELTA_L         0.45
#define  Q3TMX_MIN_DELTA_R         0.55
#define  Q3TMX_MAX_DELTA_R         0.99
#define  Q3TMX_MIN_TAU_L           0.

/* MQ Calibration constants for Newton-Raphson */
#define  Q3TMX_MQ_FIX              1
#define  Q3TMX_MQ_TOL              2
#define  Q3TMX_MQ_FIX_SAFE         3
#define  Q3TMX_MQ_TOL_SAFE         4

#define  Q3TMX_Q_STEPS             20
#define  Q3TMX_OTM_TOL             1E-7

#define  Q3TMX_M_STEPS             20
#define  Q3TMX_M_DELTA             1E-5
#define  Q3TMX_FWD_TOL             1E-5

#define  Q3TMX_S_STEPS             100
#define  Q3TMX_S_DELTA             1E-5
#define  Q3TMX_ATM_TOL             1E-4
#define  Q3TMX_S_MAX_STEP          2E-1   

#define  Q3TMX_MIN_VOL_CALIB       5E-3

/* PA constants */
#define  Q3TMX_NUM_BIN_L1          7
#define  Q3TMX_NUM_BIN_L2          7
#define  Q3TMX_NUM_BIN_R1          7
#define  Q3TMX_NUM_BIN_R2          7
#define  Q3TMX_NUM_BIN_L           Q3TMX_NUM_BIN_L1 + Q3TMX_NUM_BIN_L2
#define  Q3TMX_NUM_BIN_R           Q3TMX_NUM_BIN_R1 + Q3TMX_NUM_BIN_R2
#define  Q3TMX_EQUIK_ITVL          0
#define  Q3TMX_EQUIX_ITVL          1
#define  Q3TMX_ITVL_L1             Q3TMX_EQUIX_ITVL
#define  Q3TMX_ITVL_L2             Q3TMX_EQUIK_ITVL
#define  Q3TMX_ITVL_R1             Q3TMX_EQUIX_ITVL
#define  Q3TMX_ITVL_R2             Q3TMX_EQUIX_ITVL
#define  Q3TMX_TAIL_PROB           0.0005
#define  Q3TMX_PROB_STEP           .1
#define  Q3TMX_TAIL_ITER           10
#define  Q3TMX_TAIL_Q_L            0
#define  Q3TMX_TAIL_Q_R            0
#define  Q3TMX_QGUESS_LIMIT        10 
#define  Q3TMX_QSOLV_TOL           1E6 * DBL_EPSILON
#define  Q3TMX_XSOLV_TOL           1E3 * DBL_EPSILON
#define  Q3TMX_QSOLV_ITER          120
#define  Q3TMX_X_DELTA              .0001
#define  Q3TMX_Q_XOVER              1e-7 
#define  Q3TMX_VOLMXITER            10 
#define  Q3TMX_VOLMXSTEP            0.1
#define  Q3TMX_VOLDELTA             0.0001
#define  Q3TMX_PREQERR              0.00000001
#define  Q3TMX_MAX_ZERO_ADJ         1.0
#define  Q3TMX_LGST_MU              1.0
#define  Q3TMX_LGST_SIG             0.05
#define  Q3TMX_LGST_MINX            0.0
#define  Q3TMX_LGST_MAXX            10.0
#if (Q3TMX_NUM_BIN_R + Q3TMX_NUM_BIN_L + 6) > Q3TMX_NBQ  
#   error Too many additional gaussian nodes.
#endif

/* Bivariate corellation limits */
#define Q3TMX_MAX_CORR  0.999999
#define Q3TMX_MIN_CORR -0.999999

/* Monte Carlo Stuff */
#define Q3TMX_DEV_DIM_MAX 2

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

#define Q3TMX_BV_RAN2_SEED -9971717

/* Integration types. */
#define Q3TMX_BV_RAN2  1
#define Q3TMX_BV_SOBOL 2
#define Q3TMX_MC_NUM_PTS 100000

/* macros */


#ifndef  MAX
#define  MAX(a,b) ((a) > (b) ? (a) : (b))
#endif

#ifndef  MIN
#define  MIN(a,b) ((a) < (b) ? (a) : (b))
#endif

#ifndef  COLLAR
#define  COLLAR(a,b,x) (MIN((b), MAX((a), (x))))
#endif   

#ifndef  SHIFT_ZERO
#define  SHIFT_ZERO(a) ((fabs((a)) < TINY) ? TINY : (a))
#endif   

#ifndef  SHIFT_EPS
#define  SHIFT_EPS(a) ((fabs((a)) < 1e-10) ? 1e-10 : (a))
#endif   

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
    double  callL[Q3TMX_NBQ];          /* (I) Target call prices left of fwd   */
    double  putR[Q3TMX_NBQ];           /* (I) Target put prices right of fwd   */
    double  optATM;                 /* (I) Target ATM option price          */

    /* Internal MultiQ parameters */
    long    nbQL;                   /* (I) Mdl number of left intervals     */
    long    nbQR;                   /* (I) Mdl number of right interval     */
    double  kL[Q3TMX_NBQ];             /* (I) Mdl strike/fwd left of fwd       */
    double  kR[Q3TMX_NBQ];             /* (I) Mdl strike/fwd right of fwd      */
    double  qL[Q3TMX_NBQ];             /* (O) Mdl q left of fwd                */
    double  qR[Q3TMX_NBQ];             /* (O) Mdl q right of fwd               */
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
#define  Q3TMX_ADJ_FWD             0
#define  Q3TMX_CALL                1
#define  Q3TMX_PUT                 2
#define  Q3TMX_IN                  3
#define  Q3TMX_OUT                 4

#define  Q3TMX_VNL                 0
#define  Q3TMX_BIN                 10
#define  Q3TMX_ANN                 100   
#define  Q3TMX_TEC                 200
#define  Q3TMX_VEGA                1000
#define  Q3TMX_POW                 1100

/* bivariate option types */
#define  Q3TMX_CALL_SUM            301
#define  Q3TMX_PUT_SUM             302

#define  Q3TMX_CALL_PROD           331     
#define  Q3TMX_PUT_PROD            332     
#define  Q3TMX_CALL_PERC           341     
#define  Q3TMX_PUT_PERC            342     

#define  Q3TMX_FLR_W_FLR           305
#define  Q3TMX_FLR_W_CAP           306
#define  Q3TMX_CAP_W_FLR           307
#define  Q3TMX_CAP_W_CAP           308

#define  Q3TMX_FLR_W_FLR_EMBED     355
#define  Q3TMX_FLR_W_CAP_EMBED     356
#define  Q3TMX_CAP_W_FLR_EMBED     357
#define  Q3TMX_CAP_W_CAP_EMBED     358

#define  Q3TMX_IN_BIRIB            313
#define  Q3TMX_OUT_BIRIB           314
#define  Q3TMX_IN_SPDRIB           323
#define  Q3TMX_OUT_SPDRIB          324

#define  Q3TMX_IN_BIRIB_EPS        413
#define  Q3TMX_OUT_BIRIB_EPS       414
#define  Q3TMX_IN_SPDRIB_EPS       423
#define  Q3TMX_OUT_SPDRIB_EPS      424

#define  Q3TMX_JOINT_FWD           -999

#define  Q3TMX_BS_SPRD_CALL        -1 
#define  Q3TMX_BS_SPRD_PUT         -2
#define  Q3TMX_BS_PERC_CALL        -3
#define  Q3TMX_BS_PERC_PUT         -4
#define  Q3TMX_BS_SPRD_RATE_CALL   -11
#define  Q3TMX_BS_SPRD_RATE_PUT    -12
#define  Q3TMX_BS_PERC_RATE_CALL   -13
#define  Q3TMX_BS_PERC_RATE_PUT    -14
#define  Q3TMX_BS_PERC_SPRD_CALL   -21
#define  Q3TMX_BS_PERC_SPRD_PUT    -22

/* basis options calibration type */
#define Q3TMX_BS_CAL_FWD           1
#define Q3TMX_BS_CAL_VOL           2

/* arithmetic on option types */
#define  Q3TMX_PAY_TYPE(x)         ((x)/10 * 10)
#define  Q3TMX_COP_TYPE(x)         ((x) - Q3TMX_PAY_TYPE(x))
#define  Q3TMX_COP(cop)            (((cop) == 1) ? Q3TMX_CALL : Q3TMX_PUT)
#define  Q3TMX_COP_TO_COP(q3_cop)  (((q3_cop) == Q3TMX_CALL) ? 1. : -1.)

/* instrument description: settlement type */
#define  Q3TMX_CASH_SETL           1
#define  Q3TMX_PHYS_SETL           0      

/* instrument description: payoff parameters */
#define  Q3TMX_MAX_PAY_PARAMS      10
#define  Q3TMX_MAX_V_PARAMS         2

/*t----------------------------------------------------------------------------
 * Option type data structure 
 */
typedef struct {
    double params[Q3TMX_MAX_PAY_PARAMS];  /* Scalar parameters           */
    double *vParam[Q3TMX_MAX_V_PARAMS];   /* Vector parameters           */
    long   vParamLen[Q3TMX_MAX_V_PARAMS]; /* Vector parameter lengths    */
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
#define  Q3TMX_SMOOTH_FACTOR       2.    


/*-----------------------------------------------------------------------------
 * function prototypes
 */

/*-----------------------------------------------------------------------------
 * multiq.c
 */
int     Q3TMXMQPricer            (MQDATA       *mq,
                               long         type,
                               double       strike,
                               double       *price);

int     Q3TMXMQLevPricer         (MQDATA       *mq,
                               long         type,
                               double       strike,
                               double       leverage,
                               double       *price);

int     Q3TMXMQBinPricer         (MQDATA       *mq,    
                               long         optType, 
                               double       strike,  
                               double       *price);

int     Q3TMXMQMap               (MQDATA       *mq,
                               double       xval,
                               double       *yield);

int     Q3TMXMQSmileInit         (double       fwdRate,
                               double       sigATM,
                               double       expiry,
                               double       *smile,
                               MQDATA       *mq);

int     Q3TMXMQSmileInitFromMQ   (MQDATA       *mq1,
                               MQDATA       *mq2);                                        

void    Q3TMXMQCopyMQ            (MQDATA       *mq1,
                               MQDATA       *mq2);                                        


int     Q3TMXMQTargetSV          (SVDATA       *sv,
                               double       *mqGuess,
                               MQDATA       *mq);

int     Q3TMXMQTargetFA          (FADATA       *fa,
                               MQDATA       *pa);  

int     Q3TMXMQCalib             (MQDATA       *mq);

int     Q3TMXMQCalibMean         (MQDATA       *mq);

int     Q3TMXMQCalibMidQs        (MQDATA       *mq);

int     Q3TMXMQCalibSpread       (double       fwd,
                               double       vol,
                               double       expiry,
                               double      *smile,
                               MQDATA      *mq);

int     Q3TMXDecodeNCK           (double       nck,
                               MQDATA       *mq);

int     Q3TMXSVToMQ              (double       fwdRate,
                               double       sigATM,
                               double       expiry,
                               double       *smile,
                               MQDATA       *mq);

int     Q3TMXMQBootstrapQ        (MQDATA       *mq);
 
int     Q3TMXMQSolveMap4Q        (double       A,   
                               double       B,
                               double       *q);

/*----------------------------------------------------------------------------- 
 * stochvol.c
 */
int     Q3TMXSVPricer            (SVDATA       *sv,
                               long         type,
                               double       strike,
                               double       *price);

int     Q3TMXSVSmileInit         (double       fwdRate,
                               double       sigATM,
                               double       expiry,
                               double       *smile,
                               SVDATA       *sv);

int     Q3TMXSVCalib             (SVDATA       *sv);

int     Q3TMXSVBinPricer         (SVDATA       *sv,    
                               double       strike, 
                               double       *price);

/*-----------------------------------------------------------------------------
 * fapricer.c
 */
int     Q3TMXFAPricer            (PAYOFF       *optPayoff,
                               FADATA       *fa,
                               double       *premium);

int     Q3TMXFAGridPricer        (PAYOFF       *optPayoff, 
                               FPAYOFF      *payFunc,   
                               FADATA       *fa,        
                               double       *premium);

int     Q3TMXFATailStrikes       (FADATA       *fa,
                               double       *tailProb,   
                               double       *tailStrike);

int     Q3TMXFASmileInit         (double       expiry,              
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

int     Q3TMXFADens              (FADATA       *fa,
                               double       start,
                               double       end,
                               long         numGridPts,
                               double       *grid,
                               double       *yields,
                               double       *faDens,
                               double       *normC);

int     Q3TMXVNFMZero2Swap       (double       expiry,
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
int     Q3TMXMQGridPricer        (PAYOFF       *optPayoff,
                               FPAYOFF      *payFunc,
                               MQDATA       *mq,
                               double       *premium);

int     Q3TMXMQDens              (MQDATA       *mq,
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
int     Q3TMXSimpsonPricer1D     (PAYOFF       *optPayoff, 
                               FPAYOFF      *payFunc,
                               double       *grid,  
                               double       *yield, 
                               double       *dens,  
                               double       densNorm, 
                               double       *premium);    

int     Q3TMXPayoff1D            (PAYOFF       *optPayoff,
                               double       yield,
                               double       smoothStep,
                               double       *payoff);

/* 1d payoffs */                            
FPAYOFF Q3TMXPay1D_Yield;
FPAYOFF Q3TMXPay1D_Vnl;
FPAYOFF Q3TMXPay1D_Bin;
FPAYOFF Q3TMXPay1D_Ann;
FPAYOFF Q3TMXPay1D_Tec;
FPAYOFF Q3TMXPay1D_TecFwd;
FPAYOFF Q3TMXPay1D_FixedRibIn;
FPAYOFF Q3TMXPay1D_FixedRibOut;
FPAYOFF Q3TMXPay1D_Pow;

/* quasi-1d payoffs */
FPAYOFF Q3TMXPay1D_Sum;
FPAYOFF Q3TMXPay1D_Prod;
FPAYOFF Q3TMXPay1D_Perc;
FPAYOFF Q3TMXPay1D_YldYld;
FPAYOFF Q3TMXPay1D_YldNull;
FPAYOFF Q3TMXPay1D_VnlNull;
FPAYOFF Q3TMXPay1D_MinMaxIn;
FPAYOFF Q3TMXPay1D_MinMaxIn_PayInd; 
FPAYOFF Q3TMXPay1D_MinMaxOut;
FPAYOFF Q3TMXPay1D_MinMaxOut_PayInd;
FPAYOFF Q3TMXPay1D_MinMaxSpdIn;
FPAYOFF Q3TMXPay1D_MinMaxSpdOut;
FPAYOFF Q3TMXPay1D_MinMaxSpdIn_EPS;
FPAYOFF Q3TMXPay1D_MinMaxSpdOut_EPS;
FPAYOFF Q3TMXPay1D_FlrFlrOrCapCap;
FPAYOFF Q3TMXPay1D_FlrFlrOrCapCapEmbedFlt;
FPAYOFF Q3TMXPay1D_FlrCapOrCapFlr;
FPAYOFF Q3TMXPay1D_FlrCapOrCapFlrEmbedFlt;

/*----------------------------------------------------------------------------- 
 * grid2d.c
 */
int     Q3TMXMCPricer            (long      dimX,   
                               MQDATA  **smX,    
                               double   *rhoX,   
                               PAYOFF   *pfX,    
                               FPAYOFF  *pFcn,   
                               long      genTyp, 
                               long      npts,   
                               double   *result);

/* 2d Monte Carlo payoffs */
FPAYOFF Q3TMXPay2D_Sum;
FPAYOFF Q3TMXPay2D_Prod;
FPAYOFF Q3TMXPay2D_Perc;
FPAYOFF Q3TMXPay2D_YldYld;
FPAYOFF Q3TMXPay2D_YldNull;
FPAYOFF Q3TMXPay2D_NullYld;
FPAYOFF Q3TMXPay2D_VnlNull;
FPAYOFF Q3TMXPay2D_NullVnl;
FPAYOFF Q3TMXPay2D_MinMaxIn;
FPAYOFF Q3TMXPay2D_MinMaxOut;
FPAYOFF Q3TMXPay2D_MinMaxSpdIn;
FPAYOFF Q3TMXPay2D_MinMaxSpdOut;
FPAYOFF Q3TMXPay2D_FlrFlrOrCapCap;
FPAYOFF Q3TMXPay2D_FlrFlrOrCapCapEmbedFlt;
FPAYOFF Q3TMXPay2D_FlrCapOrCapFlr;
FPAYOFF Q3TMXPay2D_FlrCapOrCapFlrEmbedFlt;


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
double  Q3TMXNormDens              (double       x);

double  Q3TMXNormCum               (double       x);

double  Q3TMXNormCumInv            (double       prob); 

double  Q3TMXBSQInt                (double       a,
                               double       b,
                               double       x1,
                               double       x2);

int     Q3TMXBSQPricer             (double       yield,
                               double       strike,
                               double       expiry,
                               double       sigATM,
                               double       Q,
                               long         optType,
                               double       *price);

int     Q3TMXBSQImpVol             (double       yield,
                               double       strike,
                               double       expiry,
                               double       Q,
                               double       price,
                               long         optType,
                               double       volGuess,
                               double       *impVol);

double  Q3TMXExpCErrFcn            (double       a,
                               double       b,
                               double       x);


/*----------------------------------------------------------------------------- 
 * util.c
 */

double  Q3TMXSimpsIntegral       (double       step,
                               long         numPoints,
                               double       *fValues);

double  Q3TMXSmoothStepFcn       (double       x,
                               double       step);

double  Q3TMXSmoothMAX           (double       x,
                               double       step);

double  Q3TMXProbInBarriers      (double x,      /* (I) input value            */
                               double loBar,  /* (I) low barrier            */
                               double hiBar,  /* (I) high barrier           */
                               double loEps,  /* (I) low barrier smoothing  */
                               double hiEps,  /* (I) high barrier smoothing */
                               double smooth);/* (I) smoothing step size    */

double  Q3TMXProbOutBarriers     (double x,      /* (I) input value            */
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

void    Q3TMXErrMsg             (char *fmt, ...);
int     Q3TMXErrCallbackSet     (int (*cbfunc)(const char *string));


#endif


