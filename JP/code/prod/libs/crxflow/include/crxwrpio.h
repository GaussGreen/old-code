/*********************************************************************************
 * CRXWRPIO.H 
 * crx io wrappers
 *
 ********************************************************************************/

#ifndef __CRXWRPIO_H__
#define __CRXWRPIO_H__

#ifdef __cplusplus
extern "C" {
#endif

#include "l_date.h"
#include "t_curve.h"

#include "crxdata.h"
#include "crxutil.h"    
#include "crxutilio.h"

#include <alib/eqsettle.h>     /* ALib TEqStatData */
#include <alib/eqdiv.h>
#include <alib/datelist.h>     /* GtoNewDateList */
#include <alib/gtonpi.h>


/**
 * IR input structure
 */
typedef struct
{
    /* DR T_CURVE and Date -----------------------------------------------------*/
    T_CURVE  ZeroCurve[3];                /* deterministic yield curve          */

    long     BaseDate;                    /* Volatility curve base date         */
    long     SwapSt[MAX_VOL];             /* Underlying swap start              */
    long     SwapMat[MAX_VOL];            /* Underlying swap maturity           */
    long     VolDate[MAX_VOL];            /* Volatility dates                   */
    long     CorrSwapSt;
    long     CorrSwapMat;

    /* Alib TCurve and TDate ---------------------------------------------------*/
    TCurve   *AlibZeroCurve[3];            /* deterministic yield curve          */

    TDate    AlibBaseDate;                /* Volatility curve base date         */
    TDate    AlibSwapSt[MAX_VOL];         /* Underlying swap start              */
    TDate    AlibSwapMat[MAX_VOL];        /* Underlying swap maturity           */
    TDate    AlibVolDate[MAX_VOL];        /* Volatility dates                   */
    TDate    AlibCorrSwapSt;
    TDate    AlibCorrSwapMat;

    /* Misc Parameter ----------------------------------------------------------*/
    char     SwapFreq;                    /* Benchmark swap frequency           */
    char     SwapDCC;                     /* Benchmark swap day count convention*/
    char     MMB;                         /* Money market basis (360 or 365)    */

    int      NbVol;                       /* Nb of points in vol curve          */
    double   Vol[MAX_VOL];                /* Vol curve                          */
    char     VolType[MAX_VOL];            /* Caplet or swaption vol             */
    char     Freq;                        /* Frequency of underlying rate       */
    char     DCC;                         /* Day count convention               */
    int      VolUnit;                     /* Volatility unit                    */
    
    /* Model parameters --------------------------------------------------------*/
    int      NbFactor;
  
    double   Alpha[3];                    /* Relative size factors              */
    double   Beta[3];                     /* Mean reversions                    */
    double   Rho[3];                      /* Correlation                        */
    double   QLeft;                       /* Left Q mapping coeff               */
    double   QRight;                      /* Right Q mapping coeff              */
    double   FwdShift;                    /* Shift mapping coeff                */
    double   BackboneCoeff;               /* Back Bone Coefficient              */

    /* swap defn for Correlation with other assets or IR's ---------------------*/
    char     CorrSwapDCC;
    char     CorrSwapFreq;

    /* Others ------------------------------------------------------------------*/
    int      CrvToDiff;                   /* 0,1,2 for sigmaR calculation       */
    int      CrvToDisc;                   /* 0,1,2 for discount curve           */
    int      CalibFlag;                   /* FALSE => set spot vol to 1.0       */
    int      SkipFlag;                    /* TRUE => skip vols that fail        */

} IR_INPUT;


/*********************************************************************************
 * CR input structure
 ********************************************************************************/
  
typedef struct
{
  
    /* DR T_CURVE and Date -----------------------------------------------------*/
    T_CURVE  ZeroCurve;                   /* deterministic yield curve          */
    
    long     BaseDate;                    /* Volatility curve base date         */
    long     SwapSt[MAX_VOL];             /* Underlying swap start              */
    long     SwapMat[MAX_VOL];            /* Underlying swap maturity           */
    long     VolDate[MAX_VOL];            /* Volatility dates                   */

    /* credit infor */
    double     RecoveryRate;              /* Recovery Rate                      */
    char       RecoveryRateConv;          /* Recovery Rate Convention           */
  
    long       SPN;                       /* SPN                                */
    char       CtpyName[MAXBUFF];         /* counterparty name                  */
  
    char       Currency[MAXBUFF];         /* currency code                      */
    char       Seniority[MAXBUFF];        /* seniority                          */
    long       CurveId;                   /* spread curve id                    */
    long       DefDate;                   /* default date, if not defaulted, 0  */
    double     CorrTwk  ;                 /* Name Default Correlation Twk       */

    /* Alib TCurve and TDate ---------------------------------------------------*/
    TCurve   *AlibZeroCurve;              /* deterministic yield curve          */

    TDate    AlibBaseDate;                /* Volatility curve base date         */
    TDate    AlibSwapSt[MAX_VOL];         /* Underlying swap start              */
    TDate    AlibSwapMat[MAX_VOL];        /* Underlying swap maturity           */
    TDate    AlibVolDate[MAX_VOL];        /* Volatility dates                   */

    /* Misc Parameter ----------------------------------------------------------*/
    char     SwapFreq;                    /* Benchmark swap frequency           */
    char     SwapDCC;                     /* Benchmark swap day count convention*/
    char     MMB;                         /* Money market basis (360 or 365)    */

    int      NbVol;                       /* Nb of points in vol curve          */
    double   Vol[MAX_VOL];                /* Vol curve                          */
    char     VolType[MAX_VOL];            /* Caplet or swaption vol             */
    char     Freq;                        /* Frequency of underlying rate       */
    char     DCC;                         /* Day count convention               */
    int      VolUnit;                     /* Volatility unit                    */

    /* Model parameters --------------------------------------------------------*/
    int      NbFactor;

    double   Alpha;                       /* Relative size factors              */
    double   Beta;                        /* Mean reversions                    */
    double   Rho;                         /* Correlation                        */
    double   QLeft;                       /* Left Q mapping coeff               */
    double   QRight;                      /* Right Q mapping coeff              */
    double   FwdShift;                    /* Shift mapping coeff                */
    double   BackboneCoeff;               /* Back Bone Coefficient              */

    long     BaseIR_id;                   /* Id number of IR environment        */
    int      CorrAdjFlag;                 /* clean spread adj from IR/CR corr, 0-no, 1-reuse, 2-overwrite */

    /* Others ------------------------------------------------------------------*/
    int      CalibFlag;                   /* FALSE => set spot vol to 1.0       */
    int      SkipFlag;                    /* TRUE => skip vols that fail        */

    double     RecoveryDispersion;        /* 0 for no rec vol, 1 for max vol    */
    double     RecoveryBeta;              /* between -1 and +1                  */
  

} CR_INPUT;


typedef struct
{
    /* model and vol ---------------------------------------------------------*/
    T_CURVE   BasisZCurve;


    long      BaseDate;                   /* Volatility curve base date       */
    int       NbVol;                      /* Nb of points in vol curve        */
    long      SwapSt[MAX_VOL];            /* Underlying swap start            */
    long      SwapMat[MAX_VOL];           /* Underlying swap maturity         */
    long      VolDate[MAX_VOL];           /* Volatility dates                 */
    double    Vol[MAX_VOL];               /* Vol curve                        */
    char      VolType[MAX_VOL];           /* Caplet or swaption vol           */
    int       VolUnit;                    /* Volatility unit                  */
                                           
    /* base ir information ---------------------------------------------------*/
    long      BaseIR_id;
    long      LiborCrv;                   /* curve id for libor               */
    long      DiscCrv;                    /* curve id for discount            */
    char      LiborDCC;                   /* underlying libor DCC             */
    char      LiborFreq;                  /* underlying libor Freq            */
    char      BasisDCC;                   /* basis DCC                        */
    char      BasisFreq;                  /* basis Freq                       */

    char      BasisType;                  /* specified as spread or percentage*/

    /* Alib TCurve and TDate -------------------------------------------------*/
    TCurve   *AlibZeroCurve;              /* deterministic yield curve        */

    TDate    AlibBaseDate;                /* Volatility curve base date       */
    TDate    AlibSwapSt[MAX_VOL];         /* Underlying swap start            */
    TDate    AlibSwapMat[MAX_VOL];        /* Underlying swap maturity         */
    TDate    AlibVolDate[MAX_VOL];        /* Volatility dates                 */

    /* Model parameters ------------------------------------------------------*/
    int      NbFactor;

    double   Alpha;                       /* Relative size factors            */
    double   Beta;                        /* Mean reversions                  */
    double   Rho;                         /* Correlation                      */
    double   QLeft;                       /* Left Q mapping coeff             */
    double   QRight;                      /* Right Q mapping coeff            */
    double   FwdShift;                    /* Shift mapping coeff              */
    double   BackboneCoeff;               /* Back Bone Coefficient            */

} SP_INPUT;


typedef struct
{
    /* model parameters  -------------------------------------------------------*/
    double     SpotFX;                    /* in units of dom ccy / for ccy      */
    long       FXValueDate;               /* ValueDate for todays' Spot FX      */

    long       NbCompFXVols;
    long       CompVolDate[MAX_VOL];
    long       CompVolMatDate[MAX_VOL];
    double     CompVol[MAX_VOL];

    long       NbSpotFXVols;              /* not used at the moment             */
    long       SpotVolDate[MAX_VOL];
    double     SpotVol[MAX_VOL];

    /* smile parameters --------------------------------------------------------*/
    long       NbFXSmileParams;
    long       FXSmileDate[MAX_VOL];
    double     FXSmile_a1[MAX_VOL];
    double     FXSmile_a2[MAX_VOL];
    double     FXSmile_a3[MAX_VOL];

    
    long       ForeignIR_id;              /* foreign ir information             */

    /* Others ------------------------------------------------------------------*/
    int        VolBootstrapMode;
    double     FxCutOffLevel;

} FX_INPUT;


typedef struct
{
    double    SpotEq;                     /* spot equity price at ValueDate     */
    long      EqValueDate;                /* ValueDate for today's Spot EQ      */

    /* discrete dividend info --------------------------------------------------*/
    long      NbDiscDivs;                 /* number of discrete div dates       */
    long      DivKnownDate[MAX_DISC_DIV]; /* dates discrete divs are evaluated  */
    long      ExDivDate[MAX_DISC_DIV];    /* start of ex-div period             */
    long      DivPmtDate[MAX_DISC_DIV];   /* discrete div payment dates         */
    double    DollarDiv[MAX_DISC_DIV];    /* dollar div amounts                 */
    double    YieldDiv[MAX_DISC_DIV];     /* discrete yield divs (eg 0.03)      */
    char      PseudoDDivFlag;             /* Pseudo Discrete Dividend flag (Y/N)*/
    long      PseudoDDivPeriod;           /* Pseudo DDiv Period (days)          */

    /* borrow curve info -------------------------------------------------------*/
    long      NbBorrows;                  /* nb of borrow rates                 */
    long      BorrowDate[MAX_CONT_DIV];   /* borrow zero rate end dates         */
    double    BorrowZeroRate[MAX_CONT_DIV]; /* borrow zero rates                */

    /* continuous dividend info ------------------------------------------------*/
    long      NbContDivs;                 /* number of continuous div dates     */
    long      ContDivDate[MAX_CONT_DIV];  /* continuous yield div start dates   */
    double    ContDiv[MAX_CONT_DIV];      /* continuous yeild divs (eg 0.03)    */

    /* settlement info, all arrays of size NbSettDates -------------------------*/
    char      SettleType;                 /* F=Fixed, R=Rolling                 */
    long      EqDaysToSpot;               /* nb spot days to settlement         */
    long      NbSettDates;                /* nb settlem't dates(=<0 if rolling) */
    long      SettleDate[MAX_SETTLE];     /* settlement dates                   */
    long      LastTradeDate[MAX_SETTLE];  /* last trading date of settlem't prd */

    /* model and vol -----------------------------------------------------------*/
    long      NbCompEQVols;               /* number of composite eq vol points  */
    long      CompVolDate[MAX_VOL];       /* volatility dates (option expiries) */
    long      CompVolReset[MAX_VOL];      /* underlying forward maturity        */
    double    CompVol[MAX_VOL];           /* level of comp eq vol (eg 0.10)     */

    long      BaseIR_id;                  /* ID for IR of equity denomination   */

    /* system event ID's (populated by engine) ---------------------------------*/
    long      DivKnownEvID;               /* event ID for div evaluation dates  */
    long      DivPmtEvID;                 /* event ID for div payment dates     */

    /* Others ------------------------------------------------------------------*/
    int        VolBootstrapMode;          /* O=nil, 3=const                     */
    double     EqCutOffLevel;             /* eg 0.03                            */

} EQ_INPUT;

/*********************************************************************************
 * INPUT DATA
 ********************************************************************************/
typedef struct{
    int             nbIRInput;            /* num of irInputs                    */
    IR_INPUT        *irInput;

    int             nbCRInput;            /* num of crInputs                    */
    CR_INPUT        *crInput;

    double          **corrInput;          /* correlation matrix */

} CRX_INPUT;

/*******************************************************************************
 * INPUT DATA
 ******************************************************************************/
typedef struct{
    int             nbIRInput;            /* num of irInputs                 */
    IR_INPUT        *irInput;

    int             nbSPInput;            /* num of spInputs                 */
    SP_INPUT        *spInput;

    double          **corrInput;          /* correlation matrix              */

    int             PPY;
    int             nbCetIter;
    double          nbCutoff;

} BS_INPUT;

/*********************************************************************************
 * Vol term structure
 ********************************************************************************/
typedef struct
{
    long     BaseDate;                    /* Vol Base Date                      */

    /* input vol parameters ----------------------------------------------------*/
    int      NbVol;                       /* Number of vol points               */
    char     Freq;                        /* Payment frequency                  */
    char     DCC;                         /* Payment day count convention       */
    long     VolDate[MAX_VOL];            /* Vol Date                           */
    double   Vol[MAX_VOL];                /* BS Vol                             */
    long     SwapSt[MAX_VOL];             /* Underlying CDS Start Date          */
    long     SwapMat[MAX_VOL];            /* Underlying CDS Mat Date            */
    double   QLeft;                       /* Left Q mapping coef.               */
    double   QRight;                      /* Right Q mapping coef.              */
    double   Alpha;                       /* Relative size factor               */
    double   Beta;                        /* Mean reversion                     */
    double   FwdSh;                       /* Fwd shift mapping coef.            */
  
    /* numerical parameters ----------------------------------------------------*/
    int      VolUsed[MAX_VOL];            /* TRUE if used in calib.             */
    char     MeshFreq;                    /* Freq. used to calculate protection */
  
    /* calculated values -------------------------------------------------------*/
    double   Aweight[MAX_VOL];            /* Bootstrapped spot vol              */
    double   ParYield[MAX_VOL];           /* Fwd Par Spread                     */
    double   Annuity[MAX_VOL];            /* Fwd Annuity Array                  */
} VOL_DATA;

/*********************************************************************************
 *    irInput
 *    Read the following files
 *      1) the 3 zero curves
 *      2) IR info file
 *      3) IR volatility file
 *      4) IR parameter file
 *
 ********************************************************************************/
int  irInputEnvData (
                IR_INPUT  *irInput,       /* (O) ir Input Structure             */
                char      *FNameZC0,      /* (I) file name for zero curve 0     */
                char      *FNameZC1,      /* (I) file name for zero curve 1     */
                char      *FNameZC2,      /* (I) file name for zero curve 2     */
                char      *FNameInfo,     /* (I) IR info file name              */
                char      *FNameVoldiag,  /* (I) vol diagonal file name         */
                char      *FNameModlPar,  /* (I) model para file name           */
                char      *CalibFlag,     /* (I) Calib Flag                     */
                char      *NbFactorOWS,   /* (I) Nb of factor string            */
                char      *AlphaOWS,      /* (I) Alpha string                   */
                char      *BetaOWS,       /* (I) Beta string                    */
                char      *RhoOWS,        /* (I) Rho string                     */
                char      *BackboneOWS,   /* (I) Backbone string                */
                char      *SmileOWS);     /* (I) 2q string                      */

/*********************************************************************************
 *    crInput
 *    Read the following files
 *      1) the credit zero curves
 *      2) CR info file
 *      3) CR volatility file
 *      4) CR parameter file
 * 
 *    The following inputs should be hardcoded when calling this function
 *      CalibFlag       = "N"
 *      NbFactorOWS     = "1"
 *      RhoOWS          = "nil"
 *      BackboneOWS     = "0" (to avoid reading file if all other OWS are given)
 *
 ********************************************************************************/
int  crInputEnvData (
                CR_INPUT  *crInput,       /* (O) cr Input Structure             */
                char      *FNameZC,       /* (I) file name for zero curve 0     */
                char      *FNameInfo,     /* (I) CR info file name              */
                char      *FNameVoldiag,  /* (I) vol diagonal file name         */
                char      *FNameModlPar,  /* (I) model para file name           */
                char      *CalibFlag,     /* (I) Calib Flag                     */
                char      *NbFactorOWS,   /* (I) Nb of factor string            */
                char      *AlphaOWS,      /* (I) Alpha string                   */
                char      *BetaOWS,       /* (I) Beta string                    */
                char      *RhoOWS,        /* (I) Rho string                     */
                char      *BackboneOWS,   /* (I) Backbone string                */
                char      *SmileOWS);     /* (I) 2q string                      */

/*********************************************************************************
 *    spInput
 *    Read the following files
 *      1) the credit zero curves
 *      2) SP info file
 *      3) SP volatility file
 *      4) SP parameter file
 * 
 *    The following inputs should be hardcoded when calling this function
 *      CalibFlag       = "N"
 *      NbFactorOWS     = "1"
 *      RhoOWS          = "nil"
 *      BackboneOWS     = "0" (to avoid reading file if all other OWS are given)
 *
 ********************************************************************************/
int  spInputEnvData (
                SP_INPUT  *spInput,       /* (O) sp Input Structure             */
                char      *FNameZC,       /* (I) file name for zero curve 0     */
                char      *FNameInfo,     /* (I) CR info file name              */
                char      *FNameVoldiag,  /* (I) vol diagonal file name         */
                char      *FNameModlPar,  /* (I) model para file name           */
                char      *CalibFlag,     /* (I) Calib Flag                     */
                char      *NbFactorOWS,   /* (I) Nb of factor string            */
                char      *AlphaOWS,      /* (I) Alpha string                   */
                char      *BetaOWS,       /* (I) Beta string                    */
                char      *RhoOWS,        /* (I) Rho string                     */
                char      *BackboneOWS,   /* (I) Backbone string                */
                char      *SmileOWS);     /* (I) 2q string                      */

/*********************************************************************************
 *       This function does:
 *       1) Mode 0: Find and skip mode
 *           Skip everything until it finds line start with '#' (ignore spaces)
 *           Flag an error message if it can't find
 *       2) Mode 1: Skip mode
 *           Skip all the empty lines and spaces
 *           Read the first non-empty character and check it is '#' 
 *           Flag an error message if it isn't.
 *
 *       The routine returns with the file pointer at the start of the next line
 *
 ********************************************************************************/
int     FindAndSkipComLine (int      Mode,
                            FILE    *stream,
                            char     *SectionLabel,
                            char     *Routine,
                            char     *FileName);


/*********************************************************************************
 *      Functioning the same as FindAndSkipComLine() with added feature:
 *      3) Mode 2: same as Mode 0 but suppressing error message
 *      4) check optionalSeq number '0','1', ... '9',
 *      if '1' is passed in, only line begin with #1 will return SUCCESS
 ********************************************************************************/
int     FindAndSkipComLineOptional (int      Mode,
                                    FILE     *stream,
                                    char     *SectionLabel,
                                    char     *Routine,
                                    char     *FileName,
                                    char     optionalSeq);

/*********************************************************************************
 *       Read term structure input for DR Wrapper.
 *
 ********************************************************************************/
int  Term_Input_W (T_CURVE   *t_curve,    /* (O) Structure of zero curve data   */
                   char      *FileName);  /* (I) File name including extension  */

/*********************************************************************************
 *       Read volatility input
 *
 ********************************************************************************/
int     VolDiag_Input_W (
            long    *BaseDate,            /* (O) Volatility data                */
            int     *VolUnit,
            int     *NbVol,
            long    *VolDate,
            long    *SwapSt,
            long    *SwapMat,
            double  *Vol,
            char    *VolType,
            char    *FileName);           /* (I) File name including extension  */

/*********************************************************************************
 *  	Read model parameters
 *
 *       Alpha, Beta and Rho must have enough memory for MaxNbFactor before 
 *       calling this function
 *
 ********************************************************************************/
int     Param_Input (
    int        *NbFactor,                 /* (O) Model Parameters               */
    double     *Alpha,                    /* (O) NULL if MaxNbFactor = 0        */
    double     *Beta,                     /* (O) NULL if MaxNbFactor = 0        */
    double     *Rho,                      /* (O) NULL if MaxNbFactor = 0 or 1   */
    double     *BackboneCoeff,
    double     *QLeft,         
    double     *QRight,        
    double     *FwdShift,
    int        MaxNbFactor,               /* (I) Max nb of factor allowed       */
    char       *FileName,                 /* (I) File name including extension  */
    int        MAWFlg);                   /* (I) Multi-asset wrapper flag       */

int     Param_Input_Classic (
    int        *NbFactor,                 /* (O) Model Parameters               */
    double     *Alpha,                    /* (O) NULL if MaxNbFactor = 0        */
    double     *Beta,                     /* (O) NULL if MaxNbFactor = 0        */
    double     *Rho,                      /* (O) NULL if MaxNbFactor = 0 or 1   */
    double     *BackboneCoeff,
    double     *QLeft,         
    double     *QRight,        
    double     *FwdShift,
    int        MaxNbFactor,               /* (I) Max nb of factor allowed       */
    char       *FileName);                /* (I) File name including extension  */

int     Param_Input_MAW (
    int        *NbFactor,                 /* (O) Model Parameters               */
    double     *Alpha,                    /* (O) NULL if MaxNbFactor = 0        */
    double     *Beta,                     /* (O) NULL if MaxNbFactor = 0        */
    double     *Rho,                      /* (O) NULL if MaxNbFactor = 0 or 1   */
    double     *BackboneCoeff,
    double     *QLeft,         
    double     *QRight,        
    double     *FwdShift,
    int        MaxNbFactor,               /* (I) Max nb of factor allowed       */
    char       *FileName);                /* (I) File name including extension  */

/*********************************************************************************
 *       Check validity of DR Wrapper term structure inputs.
 *
 ********************************************************************************/
int     Term_Check_W (T_CURVE  *t_curve); /* (I) Structure of zero curve data   */

/*********************************************************************************
 *       Check validity of volatility inputs.
 *
 ********************************************************************************/
int     VolDiag_Check_W (long    BaseDate, /* Volatility data                   */
                         int     VolUnit,
                         int     NbVol,
                         long    *VolDate,
                         long    *SwapSt,
                         long    *SwapMat,
                         double  *Vol,
                         char    *VolType);

/*********************************************************************************
 *  	Check model parameters
 *
 ********************************************************************************/
int     Param_Check (long       NbFactor,   
                     double     *Alpha,      
                     double     *Beta,      
                     double     *Rho,       
                     double     BackboneCoeff,
                     double     FwdShift);

/*********************************************************************************
 *       Read volatility input
 *
 ********************************************************************************/
int     VolDiagCR_Input_W (
    long    *BaseDate,                    /* (O) Volatility data                */
    int     *VolUnit,
    int     *NbVol,
    long    *VolDate,
    long    *SwapSt,
    long    *SwapMat,
    double  *Vol,
    char    *FileName);                   /* (I) File name including extension  */

/*********************************************************************************
 *       Check validity of credit volatility inputs.
 *
 ********************************************************************************/
int     VolDiagCR_Check_W (
    long    BaseDate,                     /* Volatility data                    */
    int     VolUnit,
    int     NbVol,
    long    *VolDate,
    long    *SwapSt,
    long    *SwapMat,
    double  *Vol);

/*********************************************************************************
 *       This function does:
 *       1) Mode 0: Find and skip mode
 *           Skip everything until it finds line start with '###' (ignore spaces)
 *           Flag an error message if it can't find
 *       2) Mode 1: Skip mode
 *           Skip all the empty lines and spaces
 *           Read the first non-empty character and check it is '#' followed by '##' 
 *           Flag an error message if it isn't.
 *       3) Mode 2: Find and Skip mode with no error message printed
 *           Same as Mode 0, but no error message is printed if failed to find
 *           next section line.
 *           Note the status is still FAILURE in this case
 *
 *       The routine returns with the file pointer at the start of the next line
 *
 ********************************************************************************/
int     FindAndSkipSectionLine (int      Mode,
                                FILE    *stream,
                                char     *SectionLabel,
                                char     *Routine,
                                char     *FileName);

/*********************************************************************************
 *      Read summary.dat
 *
 *      Check all inputs  
 *
 ********************************************************************************/
int    summaryInputEnvData(
    char     *FNameSummary,               /* (I) Summary filename               */
    long     *Today,                      /* (O)                                */
    int      *NbIrInp,                    /* (O)                                */
    IR_INPUT **IrInp,                     /* (O)                                */
    long     *NbSpInp,                    /* (O)                                */
    SP_INPUT **SpInp,                     /* (O)                                */
    long     *NbFxInp,                    /* (O)                                */
    FX_INPUT **FxInp,                     /* (O)                                */
    long     *NbEqInp,                    /* (O)                                */
    EQ_INPUT **EqInp,                     /* (O)                                */
    int      *NbCrInp,                    /* (O)                                */
    CR_INPUT **CrInp);                    /* (O)                                */



/**
 *    Read corr matrix from correlation.dat.
 *    Memory should have been allocated already.
 *    Replaced with OverWriteString when given
 *    Check inputs
 *
 */
int corrInputEnvData(double  **corr,            /* (O) corr matrix */
                     long     NbAssets,         /* (I) nb of assets    */
                     char     *FName,           /* (I) Corr file name  */
                     char   ***CorrOWS);        /* (I) Corr overwrite string */


/**
 *    allocate memory for correlation overwrite matrix
 */

int Alloc_CorrOWS   (long       NbAssets,   /* (I) */
                     char   ****CorrOWS);   /* (O) */


/**
 *    free memory for correlation overwrite matrix
 */

int Free_CorrOWS   (long   NbAssets,       /* (I) */
                    char   ***CorrOWS);



/* free MAW Input */
void CrxMAWFreeInput(
    CRX_INPUT     *crxInput);             /* (I)  CRX Input                     */

/* free MAW Input */
void BSMAWFreeInput(
    BS_INPUT     *bsInput);             /* (I)  CRX Input                     */

/**
 * Read input from wrapper files
 * - Input
 *   -# pathDir:                                  path to input files
 * - Output
 *   -# Today:                                    today
 *   -# crxInput                                  CRX Input
 */
int CrxMAWReadInput(
    long          *Today,                 /* (O)  today                         */
    CRX_INPUT     *crxInput,              /* (O)  CRX Input                     */
    char          *pathDir);              /* (I)  path to input files           */


/**
 * Read input from wrapper files
 * - Input
 *   -# pathDir:                                  path to input files
 * - Output
 *   -# Today:                                    today
 *   -# crxInput                                  CRX Input
 */
int BSMAWReadInput(
    long          *Today,                 /* (O)  today                         */
    BS_INPUT      *bsInput,               /* (O)  basis Input                     */
    char          *pathDir);              /* (I)  path to input files           */


/*t-@CDOC(idxn=CrxTDrWrapperData,catn=structdef)---------------------
 * A data structure that holds the market information of
 * type II DR wrapper.
 */

typedef struct  {
        TDate           fToday;                 /* today's date */
        TCurve          *fDiscZcCurve;          /* Discount zero curve */
        TCurve          *fZcCurve;              /* Index zero curve */
        TCurve          *fRiskZcCurve;          /* Risky zero curve */
        TCurve          *fBasisZcCurve;         /* Basis zero curve */
        TCurve          *fBvCurve;              /* Base volatility curve */
        TCurve          *fBSVolCurve;           /* Base volatility curve */
        TSwaptionMatrix2D *fCmsSwMat;           /* Swaption matrix */
        long            fMMDenom;               /* 360 or 365 */
        TDayCount       fSwDcc;                 /* Swap day count conv */

        double          f1Beta;                 /* 1F parameters */
        double          f1Weight;
        int             f1Ppy;

        double          f2Beta1;                /* 2F parameters */
        double          f2Beta2;
        double          f2Weight1;
        double          f2Weight2;
        double          f2Corr12;
        int             f2Ppy;

        double          f3Beta1;                /* 3F parameters */
        double          f3Beta2;
        double          f3Beta3;
        double          f3Weight1;
        double          f3Weight2;
        double          f3Weight3;
        double          f3Corr12;
        double          f3Corr13;
        double          f3Corr23;
        int             f3Ppy;

} CrxTDrWrapperData;

#define CRX_DRW_TYPE2_2CURVES           (0x0001L)
#define CRX_DRW_TYPE2_3CURVES           (0x0002L)



/**------------------------------------------------
 * Read in Kapital IR wrapper environment data.
 * Support Type-2 and Type-3 wrapper, as specified by options.
 * IR volatility curve is NOT read, but can be added easily.
 */
int CrxTDrWrapperDataGetFull(
        char *pathdir,
        long options,
        CrxTDrWrapperData **that);

int CrxTDrWrapperDataFree(CrxTDrWrapperData *that);


/**--------------------------------------------------------------
 * Creates the output price file for a DRW type II.
 */
int CrxTDrWrapperDataPutPrice(double value);



/**-------------------------------------------------
 * Equaty static file data structure.
 *
 */
typedef struct _TEqStatData{

        TDividendList     *divList;
        char              settleType; 
        TEquitySettlement *stm;

}TEqStatData;


/**-------------------------------------------------
 * Equaty static file data structure.
 * Equaty settlement data structure.
 * Copied from ALIB (eqsettle.c)
 * Structure with easier working
 * data for determination of settlement dates
 * for stock.
 *
 * \index{TEqStmPrivate}
 */

typedef struct _TEqStmPrivate {

    /*
     * For rolling settlement, ie T+n days
     * (US system, S&P500)
     *
     * if tradeDate > stmPeriodDates[i], use period[i]
     * if tradeDate > stmPeriodDates[last], use period[last+1]
     *
     * 'settlementHolName' is used to decide settlement dates with T+n
     */
    TDateList   *stmPeriodDates;  /* Date when stmPeriod switches */
    long        *stmPeriods;      /* settle periods for stmPeriodDates
                                   * Must be 1 more stmPeriods than
                                   * stmPeriodDates
                                   */
    long         rollingIndex;    /* Index into stmPeriodDates[] and
                                   * stmPeriods[] */

    char        *settlementHolName;

    /*
     * Allows more efficient searching
     * GTO_EQ_STM_DIRN_INC    Increasing dates
     * GTO_EQ_STM_DIRN_NONE   No pre-arranged direction
     * GTO_EQ_STM_DIRN_DEC    Decreasing dates
     *
     */
    long         direction;


    /*
     * For fixed settlement dates as published
     * by an exchange.
     * (French system, CAC40)
     * Before first ltd use rolling stm at stmPeriodBefore
     * After last ltd use rolling stm at stmPeriodAfter
     */
    TDateList   *lastTradingDates;
    TDateList   *settlementDates;
    long         settleIndex;
    long         stmPeriodBefore;
    long         stmPeriodAfter;


} TEqStmPrivate ;


int
CrxTEqStatGet(char         *pathdir,      /* (I) tmp direcory name (or NULL) */
              char         *eqStaFnam,    /* (I) file name (equity.sta)  */
              long         busDayConv,    /* (I) Business Day Conv    */
              char         *holidayFile,  /* (I) Holiday File       */
              TEqStatData  **eqStatData); /* (O) equaty static data */

TEqStatData *
CrxNewTEqStatData(TDividendList     *divList,
                  char              settleType,
                  TEquitySettlement *stm);


void CrxFreeTEqStatData(TEqStatData *eqStatData);


/**---------------------------------------------------------------------
 * Read data from equity.dyn in London format
 * Basis of vol curve is set to 4L (quarterly)
 *
 * Returns SUCCESS/FAILURE.
 */
int CrxEqDynDataGet(   
    char    *pathdir,               /* (I) tmp direcory name (or NULL) */
    char    *eqDynFnam,             /* (I) file name (equity.dyn)  */
    double  *spotValue,             /* (O)  */
    double  *correlation,           /* (O)  */
    TCurve **indxVolCurve);         /* (O) index volatility curve */


/**---------------------------------------------------------------
 * Compute forward price given the equity static data and funding curve.
 * Currently two underlying routines are available:
 * 1) the ALIB GtoEqForwardPrice function
 * 2) based on London code with proper modification for settlement and dividend
 *    dates adjustments.
 */
int
CrxForwardPriceGen(
    double      spotPrice,      /* (I) spot asset price   */
    TEqStatData *eqStatData,    /* (I) Equity static data */
    TCurve      *zcCoF,         /* (I) Funding zero curve */
    long        busDayConv,     /* (I) Business Day Conv    */
    char        *holidayFile,   /* (I) Holiday File       */
    long        numFwd,         /* (I) Number of fwd points */
    TDate       *fwdDate,       /* (I) Forward dates to compute fwdPrice
                                 *     fwdDate[0] must be the spotDate
                                 *     where spotPrice is known  */
    double      *fwdPrice);     /* (O) Forward price      */

/*f
***************************************************************************
** Reads from DR wrapper file for CrxTBondRepoCurve
***************************************************************************
*/
CrxTBondRepoCurve* CrxBondRepoCurveDRWRead(FILE *in_fp);

/*f
***************************************************************************
** Writes to DR wrapper file for CrxTBondRepoCurve
***************************************************************************
*/
void CrxBondRepoCurveDRWWrite(FILE *out_fp, CrxTBondRepoCurve *p);


/*f
***************************************************************************
** Reads from DR wrapper file for CrxTBondPrice
***************************************************************************
*/
CrxTBondPrice* CrxBondPriceDRWRead(FILE *in_fp);

/*f
***************************************************************************
** Writes to DR wrapper file for CrxTBondPrice
***************************************************************************
*/
void CrxBondPriceDRWWrite(FILE *out_fp, CrxTBondPrice *p);

/*f
***************************************************************************
** Reads from DR wrapper file for CrxTBondPriceVolCurve
***************************************************************************
*/
CrxTBondPriceVolCurve* CrxBondPriceVolCurveDRWRead(FILE *in_fp);

/*f
***************************************************************************
** Writes to DR wrapper file for CrxTBondPriceVolCurve
***************************************************************************
*/
void CrxBondPriceVolCurveDRWWrite(FILE *out_fp, CrxTBondPriceVolCurve *p);

/*f
***************************************************************************
** Reads from DR wrapper file for CrxTBondSpreadVolCurve
***************************************************************************
*/
CrxTBondSpreadVolCurve* CrxBondSpreadVolCurveDRWRead(FILE *in_fp);

/*f
***************************************************************************
** Writes to DR wrapper file for CrxTBondSpreadVolCurve
***************************************************************************
*/
void CrxBondSpreadVolCurveDRWWrite(FILE *out_fp, CrxTBondSpreadVolCurve *p);

/*f
***************************************************************************
** Writes to DR wrapper file for TCurve
***************************************************************************
*/
void CrxWriteTCurve (FILE* fp, TCurve *tc);

/* end of extern "C" scope */
#ifdef __cplusplus
}
#endif

#endif
