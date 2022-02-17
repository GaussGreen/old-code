/******************************************************************************
 * Module:      Q3
 * Submodule:
 * File: q3.h   
 * Function:    
 * Author:      Interest Rates DR
 * Revision:    $Header$
 *****************************************************************************/
#ifndef  CRMULTIQ_H
#define  CRMULTIQ_H

#include <common/include/drmacros.h>
#include <crxflow/include/crxdata.h>
#include <crxflow/include/crxutil.h>
#include <crxflow/include/optbsq.h>

/* types for memory allocation */

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
#define  Q3_MIN_VOL             1E-8
#define  Q3_MAX_VOL             2.50
#define  Q3_MAX_VOL_VOL         2.00
#define  Q3_MIN_SKEW            0.
#define  Q3_MAX_SKEW            2.

/* MQ parameter ranges (deltas and taus) */
#define  Q3_MIN_DELTA_L         0.01
#define  Q3_MAX_DELTA_L         0.45
#define  Q3_MIN_DELTA_R         0.55
#define  Q3_MAX_DELTA_R         0.99
#define  Q3_MIN_TAU_L           0.
#define  Q3_MAX_TAU_R           1.

/* MQ Calibration constants for Newton-Raphson */
#define  Q3_MQ_FIX              1
#define  Q3_MQ_TOL              2

#define  Q3_S_STEPS             20
#define  Q3_FWD_TOL             1E-5

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

/* Bivariate epsilon limits */
#define Q3_MIN_EPS   1e-6

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

#ifndef  COLLAR
#define  COLLAR(a,b,x) (MIN((b), MAX((a), (x))))
#endif   

#ifndef  SHIFT_ZERO
#define  SHIFT_ZERO(a) ((fabs((a)) < TINY) ? TINY : (a))
#endif   

#ifndef  SHIFT_EPS
#define  SHIFT_EPS(a) ((fabs((a)) < 1e-6) ? 1e-6 : (a))
#endif   

/*t----------------------------------------------------------------------------
 * Smile measure data structure  
 */

/* MQ (multi q) parameters  */
typedef struct {
    /* Market inputs */
    double  optExpy; 					  /* (I) Expiration time              */
    double  fwdRate;					  /* (I) Forward rate                 */
    double  sigATM;						  /* (I) Mkt ATM volatility           */
    double  optATM;						  /* (I) Target ATM option price      */
	long    optType;                      /* (I) Target ATM option type       */

    /* Internal MultiQ parameters */
    long    nbQL;						  /* (I) Mdl number of left intervals */
    long    nbQR;						  /* (I) Mdl number of right interval */
    double  kL[Q3_NBQ];					  /* (I) Mdl strike/fwd left of fwd   */
    double  dL[Q3_NBQ];					  /* (I) Mdl strike/fwd left of fwd   */
    double  xL[Q3_NBQ];					  /* (I) Mdl strike/fwd left of fwd   */
    double  bL[Q3_NBQ];					  /* (I) Mdl strike/fwd left of fwd   */

    double  kR[Q3_NBQ];					  /* (I) Mdl strike/fwd right of fwd  */
    double  dR[Q3_NBQ];					  /* (I) Mdl strike/fwd left of fwd   */
    double  xR[Q3_NBQ];					  /* (I) Mdl strike/fwd right of fwd  */
    double  bR[Q3_NBQ];					  /* (I) Mdl strike/fwd right of fwd  */
    double  qL[Q3_NBQ];					  /* (O) Mdl q left of fwd            */
    double  qR[Q3_NBQ];					  /* (O) Mdl q right of fwd           */

    double  sigMQ;						  /* (O) Mdl vol of normal dist       */
	double  optMQ;						  /* (O) calibrated atm opt price     */
	double  C;							  /* (O) 1/E[H(x)]                    */
	double  K;							  /* (O) F(x,C,K)=K(CH(x)-1)+1        */

    /* Calibration setup (encoded in NCK) */
    long    calibType;					  /* (I) Calibration type (FIX or TOL)*/
    long    sSteps;						  /* (I) Number of NR steps for sigMQ */

    double  sDelta;						  /* (I) Step fog sigMQ num derivative*/    
    double  fwdTol;						  /* (I) Calib tolerance for forward  */
    double  atmTol;						  /* (I) Calib tolerance for ATM vol  */  

    /* Other selectable behavior */
    long    calcFwd;					  /* (I) Calc. fwd, ignore struct value*/       

	double muMQ;                          /* to be deleted                    */
} MQDATA;

/* FA (forward adjusted) parameters */
typedef struct {
    double  freqSwap;                     /* (I) Frequency of swap rate       */
	double  freqAnn;                      /* (I) Frequency of cds spread      */
    double  matAnn;                       /* (I) Maturity (tenor) of CDS spread*/
    double  matDel;                       /* (I) Delay interval               */
    double  alphaAnn;                     /* (O) Annuity alpha for ann zero 
												 rate                         */
    double  alphaDel;                     /* (O) Annuity alpha for delay zero 
												 rate                         */
    double  powerAnn;                     /* (O) VNFM power for annuity calc  */
    double  powerDel;                     /* (O) VNFM power for delay calc    */
	double  recovery;
	double  ryAnnuity;
	double  fwdRate;                      /* (O) fwd par spread               */
	double  IRate;					      /* (O) equivalent IR (for protection 
												 leg)				          */

	double  zeroSwapVolRatio;             /* (O) ratio between zero rate to delay 
												 vol and vol of cds spread    */ 
	double  accrFactor;                   /* (I) ratio between accr and no accr
											     fwd par spreads              */
    MQDATA  *mq;                          /* (O) Underlying annuity measure MQ*/ 
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
 * crmultiq.c
 */

int     Q3MQPricerCR          (MQDATA       *mq,
                               long         type,
                               double       strike,
                               double       *price);

int Q3MQCAPricerCR(
	MQDATA            *mq,           /* (I) MQ data  */
	FADATA			  *fa,
	long              optType,       /* (I) option type */
	double            strike,        /* (I) option strike */
	double	          zeroRatio,
	double		      volRatio,
	double            *price);        /* (O) option price & vega */

MQDATA* Q3MQMakeCR(
/** Forward value */
    double     fwd,
/** At the money volatility */
    double     vol,
/** Option type, i.e. call or put. Use CRX_OPTION_TYPE_... */
    long       optType,
/** Time to expiry (in years) */
    double     time,
/** Multi Q-distribution user parameters */
    CrxTQDist *qdist);

int Q3MQInitCR(MQDATA     *mq,			  /* (O) MQDATA structure             */
			   double     fwdRate,		  /* (I) fwd rate                     */
			   double     sigATM,		  /* (I) atm vol                      */
			   long       optType,		  /* (I) atm opt type                 */
			   double     optMat,         /* (I) atm opt mat                  */
			   int        numQ,			  /* (I) number of q's                */
			   double     *q,			  /* (I) a list of q's                */
			   int        numD,			  /* (I) number of d's                */
			   double     *d);			  /* (I) a list of d's                */

int Q3MQUpdateCR(MQDATA     *mq);			  /* (O) MQDATA structure             */
int Q3MQSetCandKCR(MQDATA   *mq);		  /* (I) MQ structure                  */
int     Q3MQSetCCR              (MQDATA       *mq);

/** Computes mapping function $y=H(x)$ where x is drawn from N(0,1) */
int     Q3MQMapCR             (MQDATA       *mq,
                               double        xstd,
                               double       *yield);

/** Computes inverse function $x=H^{-1}(y)$ where x is returned as a number
    of standard deviations from 0 (i.e. as if drawn from N(0,1)) */
int Q3MQMapInverseCR(
    MQDATA            *mq,
    double             yval,
    double            *xstd);

/**
 * Computes the probability density in y-space.
 */
int Q3MQDensityCR(
    MQDATA            *mq,
    double             yval,
    double            *probDensity);
    
/**
 * Computes the cumulative probability in y-space.
 */
int Q3MQCumCR(
    MQDATA            *mq,
    double             yval,
    double            *cumProbability);

/*-----------------------------------------------------------------------------
 * fapricer.c
 */

double Zero(double parSpread, FADATA *fa);

double Protection(double parSpread, double zeroRatio, FADATA *fa);


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
							   long         freqSwap,
                               long         freqCDS,                
                               double       swapMat,           
                               double       swapRate,
							   double       cdsParSpread,
							   double		recovery,
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

int     Q3VNFMZero2SwapCR     (double       expiry,
                               double       volStart,
                               double       rateStart,
                               long         freqSwap,
                               double       swapMat,
                               double       swapRate,
							   double		recovery,
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

/* quasi-1d payoffs */

FPAYOFF Q3Pay1D_VnlNull;


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

double  Q3ProbInBarriersSmooth(double x,      /* (I) input value            */
                               double loBar,  /* (I) low barrier            */
                               double hiBar,  /* (I) high barrier           */
                               double loEps,  /* (I) low barrier smoothing  */
                               double hiEps,  /* (I) high barrier smoothing */
                               double smooth);/* (I) smoothing step size    */

double  Q3ProbOutBarriersSmooth(double x,     /* (I) input value            */
                               double loBar,  /* (I) low barrier            */
                               double hiBar,  /* (I) high barrier           */
                               double loEps,  /* (I) low barrier smoothing  */
                               double hiEps,  /* (I) high barrier smoothing */
                               double smooth);/* (I) smoothing step size    */

int     Q3SwapDouble          (double *,      /*(I/O) Swop two elements     */
                               double *);     

int     Q3SwapInt             (int *,         /*(I/O) Swop two elements     */
                               int *);     


#endif
