/* ======================================================
   FILENAME:  uttypes.h

   PURPOSE:   A few useful types to replace strings
   ====================================================== */

#ifndef UTTYPES_H
#define UTTYPES_H

/* ---------------------------------------------------------------------------
                                     BOOLEAN
   ---------------------------------------------------------------------------
 */

typedef enum SRT_Boolean_ {
  SRT_NO = 0,
  SRT_FALSE = 0,
  SRT_YES = 1,
  SRT_TRUE = 1
} SRT_Boolean;

typedef enum Forback_ { SRT_FOR = 0, SRT_BACK = 1, SRT_FORBACK = 2 } Forback;

/* ---------------------------------------------------------------------------
                                        DATES
   ---------------------------------------------------------------------------
 */

typedef enum WeekDay_ {
  SATURDAY = 0,
  SUNDAY,
  MONDAY,
  TUESDAY,
  WEDNESDAY,
  THURSDAY,
  FRIDAY,
  SAT = 0,
  SUN,
  MON,
  TUE,
  WED,
  THU,
  FRI
} WeekDay;

typedef enum SrtMonth_ {
  SRT_JAN = 1,
  SRT_FEB,
  SRT_MAR,
  SRT_APR,
  SRT_MAY,
  SRT_JUN,
  SRT_JUL,
  SRT_AUG,
  SRT_SEP,
  SRT_OCT,
  SRT_NOV,
  SRT_DEC
} SrtMonth;

typedef enum SrtUnit_ {
  SRT_DAY = 0,
  SRT_BDAY,
  SRT_WEEK,
  SRT_MONTH,
  SRT_YEAR
} SrtUnit;

typedef enum BusDayConv_ {
  NO_BUSDAY_CONVENTION = 0,
  MODIFIED_SUCCEEDING,
  SUCCEEDING,
  LASTBUSDAYCONV
} BusDayConv,
    SrtBusDayConv;

typedef enum HolidayConv_ { NO_HOLIDAYS, LASTHOLIDAYCONV } HolidayConv;

/* ---------------------------------------------------------------------------
             SWAP CONVENTIONS: BASIS      , COMPOUNDING
   ---------------------------------------------------------------------------
 */

typedef enum SrtCompounding_ {
  SRT_SIMPLE = 0,
  SRT_ANNUAL,
  SRT_SEMIANNUAL,
  SRT_QUARTERLY = 4,
  SRT_MONTHLY = 12
} SrtCompounding;

#define ISCompounding(x)                                                       \
  ((x) == SRT_ANNUAL || (x) == SRT_SEMIANNUAL || (x) == SRT_QUARTERLY ||       \
   (x) == SRT_MONTHLY)

typedef enum BasisCode_ {
  BASIS_ACT_ACT = 0,
  BASIS_ACT_365,
  BASIS_ACT_360,
  BASIS_30_360,
  BASIS_30_360E,
  BASIS_ACT_USD,
  LASTBASISCODE
} BasisCode,
    SrtBasisCode;

#define ISBasisCode(x) ((x) > -1 && (x) < LASTBASISCODE)

typedef enum SrtReceiverType_ {
  SRT_PAYER = 0,
  SRT_RECEIVER,
  SRT_FORWARD,
  SRT_LASTRECEIVERTYPE
} SrtReceiverType;

typedef enum InterpMethod_ {
  LIN_R,
  LIN_RT,
  LASTINTERPMETHOD
} InterpMethod,
    SrtInterpMethod;

typedef enum StructType_ {
  BOND_OPTION,
  SWAPTION,
  CAPFLOOR,
  RESETCAPFLOOR,
  RESETCMSOPTION,
  SWAP,
  BOND,
  ///// added by Albert Wang 08/25/03 - begin
  SIMPLEMIDAT,
  ///// added by Albert Wang 08/25/03 - end

  LASTSTRUCTTYPE
} StructType;

/* ---------------------------------------------------------------------------
             OPTIONS : CALL/PUT      , UP/DOWN...
   ---------------------------------------------------------------------------
 */

typedef enum SrtCallPutType_ {
  SRT_CALL = 0,
  SRT_PUT,
  SRT_STRADDLE,
  SRT_LASTCALLPUTTYPE
} SrtCallPutType;

typedef enum SrtDiffusionType_ {
  SRT_LOGNORMAL = 0,
  SRT_NORMAL,
  SRT_LOGNORMAL_ATS,
  SRT_NORMAL_ATS,
  SRT_BETAVOL,
  SRT_HESTONVOL,
  SRT_SABRVOL,
  SRT_BVMVOL,
  SRT_BVM2VOL,
  SRT_BVMHVOL,
  SRT_BVMH2VOL,
  SRT_BVMCVOL,
  SRT_QUADRAVOL,
  SRT_LOG_QUADRAVOL,
  SRT_SHIFTED_LOG_VOL,
  SRT_BMMVOL,
  SRT_BMM2VOL,
  SRT_BMM3VOL,
  SRT_BMMBVMVOL,
  SRT_BMMBVM2VOL,
  SRT_BMMBVMHVOL,
  SRT_BMMBVMH2VOL,
  SRT_BMMBVMCVOL,
  SRT_BMMCALVOL,
  SRT_BMM2CALVOL,
  SRT_BMM3CALVOL,
  SRT_LASTDIFFUSIONTYPE
} SrtDiffusionType;

typedef enum SrtBarrierType_ {
  SRT_DOWN = 0,
  SRT_UP,
  SRT_LASTBARRIERTYPE
} SrtBarrierType;

typedef enum SrtMinmaxType_ {
  SRT_MIN = 0,
  SRT_MAX,
  SRT_LASTMINMAXTYPE
} SrtMinmaxType;

typedef enum SrtLadderType_ {
  SRT_GRADUATED = 0,
  SRT_BEST,
  SRT_LASTLADDERTYPE
} SrtLadderType;

/* ----------------------------------------------------------------------------
                                RESET_OR_OPTIMISED
--------------------------------------------------------------------------------
*/
typedef enum SrtResOptType_ {
  SRT_AUTORESET = 0,
  SRT_OPTIMISED,
  SRT_LASTRESOPTTYPE
} SrtResOptType;

typedef enum SrtPriceType_ { SRT_PREMIUM = 0, SRT_VOLATILITY } SrtPriceType;

/*---------------------------------------------------------------------------------
                           BEST OR WORST OF TWO OPTIONS
  ---------------------------------------------------------------------------------*/
typedef enum SrtBestWorstType_ {
  SRT_BESTOFF = 0,
  SRT_WORSTOFF
} SrtBestWorstType;

/* ---------------------------------------------------------------------------
                               UNDERLYINGS
   ---------------------------------------------------------------------------
 */

typedef enum SrtDomForType_ {
  SRT_DOMESTIC = 0,
  SRT_FOREIGN,
  SRT_LASTDOMFORTYPE
} SrtDomForType;

/* ---------------------------------------------------------------------------
     Calibration Types for BGM
   ---------------------------------------------------------------------------
 */

typedef enum BGMCalibType_ {
  BGM_SDP = 0,
  BGM_CASC,
  BGM_LASTCALIBTYPE
} BGMCalibType;

typedef enum BGMCorrelType_ {
  BGM_SLIDING = 0,
  BGM_CONVERGING,
  BGM_LASTCORRELTYPE
} BGMCorrelType;

/* ---------------------------------------------------------------------------
     Component Types for BGMSABRGetVol
   ---------------------------------------------------------------------------
 */

typedef enum SABRVolComponent_ {
  SABR_ATMLOG = 0,
  SABR_ATMNORM,
  SABR_BETAVOL,
  SABR_LOGVOL,
  SABR_NORMVOL,
  SABR_ALPHA,
  SABR_BETA,
  SABR_RHO,
  SABR_ISMARKETSABR,
  SABR_ZETA
} SABRVolComponent;

/*---------------------------------------------------------------------------------
                                   FIRST_OR_SECOND DERIVATIVE
  ---------------------------------------------------------------------------------*/
typedef enum SrtDerType_ {
  SRT_FIRSTDER = 0,
  SRT_SECONDDER,
  SRT_TERDER,
  SRT_LASTDERTYPE
} SrtDerType;

typedef enum SrtShiftType_ {
  SRT_SHIFTVOL = 0,
  SRT_SHIFTYC,
  SRT_LASTSHIFTTYPE
} SrtShiftType;

/*========================FOR DELTA REPORT
 * ====================================*/

typedef enum SrtHedgeType_ { SRT_SWAP = 0, SRT_FRA, SRT_FUT } SrtHedgeType;

typedef enum SrtUndFRAType_ {
  SRT_UNDFRA_SWAP = 0,
  SRT_UNDFRA_FUT
} SrtUndFRAType;

/*---------------------------------------------------------------------------------*/
/*                   Volatility conversions */
/*---------------------------------------------------------------------------------*/

typedef enum SrtVolConversion_ {

  LOG_TO_SABR = 0,
  LOG_TO_HESTON,
  LOG_TO_NORMAL,
  LOG_TO_LOG,
  NORMAL_TO_LOG,
  NORMAL_TO_SABR,
  NORMAL_TO_HESTON,
  NORMAL_TO_NORMAL,
  HESTON_TO_LOG,
  HESTON_TO_NORMAL,
  HESTON_TO_SABR,
  HESTON_TO_HESTON,
  SABR_TO_LOG,
  SABR_TO_NORMAL,
  SABR_TO_HESTON,
  SABR_TO_SABR

} SrtVolConversion;

/*---------------------------------------------------------------------------------*/
/*                   Calibration Types */
/*---------------------------------------------------------------------------------*/

typedef enum SrtCalibrationType_ {

  MATCH_CONV,
  CHI2_MIN,
  MATCH_3STRIKES

} SrtCalibrationType;

typedef enum SrtInflationIdxBondAssetSwapType_ {
  SRT_FLAT_NOTIONAL = 0,
  SRT_FLAT_INDEXED_NOTIONAL,
  SRT_INDEXED_NOTIONAL

} SrtInflationIdxBondAssetSwapType;

#endif