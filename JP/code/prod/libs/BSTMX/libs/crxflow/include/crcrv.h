/***************************************************************
* Module:      cranalyst
* File:        crcrv.h
* Description: curve analytics
***************************************************************/
#ifndef _crcrv_H
#define _crcrv_H

#include <common/include/drmacros.h>
#include <common/include/drtypes.h>

#include <alib/tcurve.h>
#include <alib/ldate.h>
#include <alib/dcconvo.h>

#include <alib/datelist.h>
#include <alib/convert.h>
#include <alib/date_sup.h>
#include <alib/gtozc.h>
#include <string.h>

#define MAXPOINTS 500

#define INTERP_METHOD GTO_FLAT_FORWARDS

/*=============================================================================
* ENUMS
*============================================================================*/

/**\enum KProtPayConv
* This enum defines two possibilities: pay on default, and pay at maturity.*/
typedef enum 
{
    PAY_DEF, 
    PAY_MAT
} KProtPayConv;

/** \enum KAccrualConv
*  Flag for accrued fee upon default
*    - NONE
*    - ALL
*/
typedef enum
{
    ACCRUAL_PAY_NONE,
    ACCRUAL_PAY_ALL
} KAccrualConv;

/*=============================================================================
* CDS PROTECTION LEG STRUCT, CONSTRUCTOR & DESTRUCTOR
*============================================================================*/

/**\struct KProtLeg_D
* This structure defines a CDS protection leg.
*    - mStDate:        start date of the protection period.
*    - mNbPrd:         number of periods with different notional and recovery
*    - mDates:         date schedule
*    - mNotionals:     notional schedule.
*    - mRecoveries:    recovery rate schedule
*    - mPayType:       pay upon default (PAY_DEF) or at maturity (PAY_MAT).
*    - mDelay:         default payment delay interval
*    - mFrequency:     integration interval
*/
typedef struct
{
    TDate         mStDate;
    int           mNbPrd;
    TDate         *mDates;
    double        *mNotionals;
    double        *mRecoveries;
    KProtPayConv  mPayType;
    TDateInterval mDelay;
    TDateInterval mFrequency;
} KProtLeg_D;

/**ProtectionCreate
* Creates a protection leg and returns a pointer. Creates deep copies of arrays,
* so caller is responsible for freeing date, notionals & recoveries
* Return valid pointer, or NULL if unsuccessful.*/
KProtLeg_D* ProtectionCreate(
    TDate         stDate,        /**<(I)   Start date of protection          */
    int           nbPrd,         /**<(I)   Number of protection periods      */
    const TDate   *dates,        /**<(I)   End dates for protection schedules*/
    const double  *notionals,    /**<(I)   Notional for each interval        */
    const double  *recoveries,   /**<(I)   Recovery for each interval        */
    KProtPayConv  payType,       /**<(I)   Timing of protection payment      */
	TDateInterval delay,         /**<(I)   Payment delay (as offset from def)*/
	TDateInterval frequency      /**<(I)   Protection integration interval   */
);

/**Free memory for a KProtLeg_D*/
int CrxProtectionFree(
    KProtLeg_D *protLeg);

/*=============================================================================
* CDS FEE LEG STRUCT, CONSTRUCTORS & DESTRUCTOR
*============================================================================*/

/** \struct KFeeLeg_D
* This structure defines a CDS fee leg.
*    - mNbCF:          nb of fee payment periods.
*    - mAccStDates:    accrual start dates.
*    - mAccEndDates:   accrual end dates.
*    - mPayDates:      floating coupon payment dates.
*    - mNotionals:     notionals.
*    - mCoupons:       CDS fee
*    - mDCC:           day count convention for coupon accrual.
*        -# "30/360".
*        -# "ACT/360".
*        -# "ACT/365".
*    - mAccrualConv    accrued upon default payment convention
*    - mFrequency:     integration step (in days) for default accruals
*
*/
typedef struct
{
    int           mNbCF;
    TDate         *mAccStDates;
    TDate         *mAccEndDates;
    TDate         *mPayDates;
    double        *mNotionals;
    double        *mCoupons;
    TDayCountConv mDCC;
	KAccrualConv  mAccrualConv;
	TDateInterval mFrequency;

} KFeeLeg_D;

/**RiskyFeeCreate
* Creates a risky fee leg of a CDS from an explicit set of payment dates, 
* coupons, and notionals. The dates, notionals and coupon pointers are
* copied deeply into the resulting object, so the caller is responsible
* for freeing all array parameters to this constructor.
* Return NULL if failed, or a pointer to the object.*/
KFeeLeg_D* RiskyFeeCreate(
    int           nbCF,         /**<(I) No. of fee payments                  */
    const TDate   *accStDates,  /**<(I) Fee start dates                      */
    const TDate   *accEndDates, /**<(I) Fee end dates                        */
    const TDate   *payDates,    /**<(I) Fee payment dates                    */
    const double  *notionals,   /**<(I) Notionals                            */
    const double  *coupons,     /**<(I) Fee coupon rates                     */
    TDayCountConv dcc,          /**<(I) Fee daycount convention              */
    KAccrualConv  payAccr,      /**<(I) Pay accrued on default?              */
	TDateInterval frequency     /**<(I) Integration step for default accruals*/
);

/**CrxFeeLegCreateFromFreq
* Creates a CDS fee leg from a start date, maturity, and coupon frequency
* Return NULL if failed, or a pointer to the object.
*/
KFeeLeg_D* CrxFeeLegCreateFromFreq(
    TDate         sttDate,       
    TDate         endDate,
    TDateInterval couponFrequency,
    KStubLocation stub,
    double        notional,
    double        coupon,
    TDayCountConv dcc,
    KAccrualConv  payAccr,
	TDateInterval frequency);

/**Free memory for a KFeeLeg_D pointer.*/
int CrxFeeLegFree(
    KFeeLeg_D *feeLeg);

/*=============================================================================
* CDS STRUCT, CONSTRUCTOR & DESTRUCTOR
*============================================================================*/

/**Defines a CDS structure as a combination of a protection 
* leg and a fee payment leg.
*/
typedef struct
{
	KFeeLeg_D  *feeLeg;
	KProtLeg_D *protLeg;
} KCDS_D;

/**
 *  Create CDS structure from explicit payment dates and default terms.
 */
KCDS_D* CDSCreate(
    int           nbFeePeriods,    /**<(I) Number of fee payments            */
    const TDate   *accStDates,     /**<(I) Accrual start dates for payments  */
    const TDate   *accEndDates,    /**<(I) Accrual end dates for payments    */
    const TDate   *payDates,       /**<(I) Payment dates for fee payments    */
    const double  *notionals,      /**<(I) Notionals for each period         */
    const double  *coupons,        /**<(I) Coupon rates for each period      */
    const double  *recoveries,     /**<(I) recovery rates for coupon period  */
    TDayCountConv dcc,             /**<(I) Fee day-count convention          */
    KAccrualConv  payAccr,         /**<(I) Default treatment of fee accrual. */
    TDate         protStDate,      /**<(I) Protection start date             */
    KProtPayConv  payType,         /**<(I) Default payment convention        */
    TDateInterval delay,           /**<(I) Payment delay                     */
    TDateInterval frequency        /**<(I) Integration freq. for prot. leg   */
);

/**
 *  Free CDS structure
 */
int    CDSFree(
    KCDS_D        *cds);            /**<(I) Pointer to CDS structure         */

/*=============================================================================
* OTHER FUNCTIONS
*============================================================================*/

/**
 *  Risky discount factor
 * - startDate:      (I) start date 
 * - endDate:        (I) end   date
 * - irCurve:        (I) ir curve
 * - crCurve:        (I) cr curve
 */
double RiskyDiscountFactor(
    TDate           startDate,  /**<(I) Start date to discount to        */
    TDate           endDate,    /**<(I) End date to discount from        */
    const TCurve    *irCurve,   /**<(I) IR ZC Curve                      */
    const TCurve    *crCurve    /**<(I) Clean spread credit curve        */
);

/**
* Compute the value of a protection leg
* - pv:           (O) pv of protection leg
* - valDate:      (I) valuation date
* - stDate:       (I) start date of the protection period.
* - payType:      (I) pay upon default (PAY_DEF) or at maturity (PAY_MAT).
* - delay:        (I) payment delay for protection leg
* - nbPrd:        (I) number of periods with different notional and recovery
* - dates:        (I) date schedule
* - notionals:    (I) notional schedule.
* - recoveries:   (I) recovery rate schedule
* - frequency:    (I) integration step (in days) for protection leg
* Return SUCCESS/FAILURE.
*/
int ProtectionPV(
    double        *pv,              /**< (O) pv of protection leg            */
    TDate         today,            /**< (I) compute price as of date        */
    TDate         valDate,          /**< (I) return price as of date         */
    TDate         stDate,           /**< (I) start date of the protection    */
    int           nbPrd,            /**< (I) number of periods               */
    const TDate   *dates,           /**< (I) date schedule                   */
    const double  *notionals,       /**< (I) notional schedule.              */
    const double  *recoveries,      /**< (I) recovery rate schedule          */
    KProtPayConv  payType,          /**< (I) PAY_DEF or PAY_MAT.             */
    TDateInterval delay,            /**< (I) payment delay for protection leg*/
    TDateInterval frequency,        /**< (I) Integration interval            */
    const TCurve  *discZC,          /**< (I) Interest rate zero curve        */
    const TCurve  *sprdZC           /**< (I) Clean spread zero curve         */
);

/**
* This is the object version of protection leg valuation
*/
int ProtectionPV_O(
    double              *pv,        /**<(O) pv of protection leg             */
    TDate               today,      /**<(I) compute price as of date         */
    TDate               valDate,    /**<(I) return price as of date          */
    const KProtLeg_D    *protLeg,   /**<(I) protection leg object            */
    const TCurve	    *discZC,    /**<(I) Interest rate zero curve         */
    const TCurve        *sprdZC     /**<(I) Clean spread zero curve          */
);

/**
* Compute the value of a risky fee leg
* - pv:           (O) pv of risky fee leg
* - valDate:      (I) valuation date
* - nbCF:         (I) nb of fee payment periods.
* - accStDates:   (I) accrual start dates.
* - accEndDates:  (I) accrual end dates.
* - payDates:     (I) coupon payment dates.
* - notionals:    (I) notionals.
* - coupons:      (I) CDS fee schedule
* - dcc:          (I) day count convention for coupon accrual.
*        -# "30/360".
*        -# "ACT/360".
*        -# "ACT/365".
* - payAccr:      (I) accrued upon default payment convention
* - frequency:    (I) integration step (in days) for default accruals
* - priceConv     (I)
* - discBaseDate: (I) value date of IR discount curve
* - nbDiscZC:     (I) nb of IR discount zero curve points
* - discZCDates:  (I) IR discount zero date array
* - discZCRates:  (I) IR discount zero rate array
* - sprdBaseDate: (I) value dates of clean spread curve
* - nbSprd:       (I) nb of clean spread zero curve points
* - sprdZCDates:  (I) clean spread date array
* - sprdZCRates:  (I) clean spread array
*
* Return SUCCESS/FAILURE.
*
*/
int RiskyFeePV(
	double        *pv,          /**<(O) price of fee leg                     */
    TDate         today,        /**<(I) compute price as of date             */
	TDate         valDate,      /**<(I) return price as of date              */
    int           nbCF,         /**<(I) No. of fee payment periods.          */
    const TDate   *accStDates,  /**<(I) Fee accrual start dates.             */
    const TDate   *accEndDates, /**<(I) Fee accrual end dates                */
    const TDate   *payDates,    /**<(I) Fee payment dates                    */
    const double  *notionals,   /**<(I) Fee notionals                        */
    const double  *coupons,     /**<(I) Fee coupon rates                     */
    TDayCountConv dcc,          /**<(I) Fee daycount convention              */
    KAccrualConv  payAccr,      /**<(I) Pay accrued on default?              */
	KStubType     priceConv,    /**<(I) BOND, NONE, or SIMPLE                */
	TDateInterval frequency,    /**<(I) Integration frequency                */
    const TCurve  *discZC,      /**<(I) Interest rate zero curve             */
    const TCurve  *sprdZC       /**<(I) Clean spread zero curve              */
);

/**
* This is the object version of fee leg PV calculation
*/
int RiskyFeePV_O(
    double          *pv,        /**<(O) PV output                            */
    TDate           today,      /**<(I) compute price as of date             */
	TDate           valDate,    /**<(I) return price as of date              */
    const KFeeLeg_D *feeLeg,    /**<(I) input fee leg struct                 */
    KStubType       priceConv,  /**<(I) BOND, NONE, or SIMPLE                */
    const TCurve    *discZC,    /**<(I) Interest rate zero curve             */
    const TCurve    *sprdZC     /**<(I) Clean spread zero curve              */
);

/**
* Compute the adjustment factor to par spread for pay accrue upon default
* for a given accrual period.
* Return SUCCESS/FAILURE.
*/
int PayAccrOnDefaultAdj(
    double        *adjFactor,      /**<(O) adj. factor as % of fee coupon    */
    TDate         today,           /**<(I) compute price as of date          */
	TDate         valDate,         /**<(I) return price as of date           */
    TDate         accStDate,       /**<(I) accrual start date                */
    TDate         accEndDate,      /**<(I) accrual end date                  */
    TDayCountConv dcc,             /**<(I) DCC (30/360, ACT/360, ACT/365)    */
    double        recovery,        /**<(I) recovery rate                     */
	TDateInterval frequency,
    const TCurve  *discZC,         /**<(I) Interest rate zero curve          */
    const TCurve  *sprdZC          /**<(I) Clean spread zero curve           */
);

/**
 * Object version of CDS pricing
 * Return SUCCESS/FAILURE.
 */
int CDSPV_O(
    double		    *pv,        /**<(O) calculated PV, if SUCCESS            */
    TDate           today,      /**<(I) calculation date                     */
	TDate		    valDate,    /**<(I) Value date                           */
    const KCDS_D    *CDS,       /**<(I) Valid pointer to CDS to price        */ 
	KStubType       priceConv,  /**<(I) CDS stub type                        */
    const TCurve    *discZC,    /**<(I) Interest rate zero curve             */
    const TCurve    *sprdZC     /**<(I) Clean spread zero curve              */
);

/**
*  CreditZCBootstrap
*
*  This function creates a zero curve which prices to par a series of
*  CDS's with aligned cash-flows. It can also be used to bootstrap CDS's
*  with non-aligned cash-flows, by calling it for one such CDS at a time,
*  passing the curve constructed by the previous call as input.
*
*  Returns the bootstrapped curve, if successful, else NULL. The
*  caller is responsible for freeing the curve.
**/
TCurve* CreditZCBootstrap(
	TDate         today,            /**<(I) Trade date for new curve         */
	TDate		  valDate,          /**<(I) Spot-date for new curve          */
	TDayCountConv dcc,              /**<(I) Day-count for par spreads        */
	KAccrualConv  accrual,          /**<(I) Accrual convention               */
    KProtPayConv  payType,          /**<(I) Payment convention               */
	TDateInterval delay,            /**<(I) Payment delay for par spreads    */
    int           nbCouponDates,    /**<(I) No. of dates                     */
    const TDate   *accSttDates,     /**<(I) Accrual period start dates       */
    const TDate   *accEndDates,     /**<(I) Accrual period end dates         */
    const TDate   *payDates,        /**<(I) Payment dates                    */
    int           nbCDS,            /**<(I) No. of par CDS rates             */
    const int     *maturityIdx,     /**<(I) Index of maturity in cpn dates   */
	const double  *parFees,         /**<(I) Par fees                         */
	double        recovery,         /**<(I) Recovery rate (same for all)     */
    TDateInterval frequency,        /**<(I) Integration frequency for prot.  */
	const TCurve  *irCurve,         /**<(I) IR ZC Curve                      */
	const TCurve  *crCurveStub      /**<(I) Clean spread curve stub, or NULL */
    );

/**
* CrxBuildCurveFromFlatParSpreadOffset
* Builds a clean spread curve by offsetting the par-spreads by a fixed amount 
* relative to those implied by a given credit curve (crCurve). If crCurve is NULL,
* instead of offsetting by a fixed basis, a flat par-spread curve is built.
* Returns a pointer to the new curve if successful, else returns a NULL pointer; the
* caller is responsible for freeing the returned curve.
*/
TCurve* CrxBuildCurveFromFlatParSpreadOffset (
	TDate today,              /**<(I) Trade date for new curve               */
	TDate valueDate,          /**<(I) Spot-date for new curve                */
	TDayCountConv dcc,        /**<(I) Day-count convention for par spreads   */
	KAccrualConv accrual,     /**<(I) Accrual convention for par spreads     */
	KProtPayConv payType,     /**<(I) Payment convention for par spreads     */
	TDateInterval delay,      /**<(I) Payment delay for par spreads          */
	double parSpread,         /**<(I) Offfset to to crCurve or flat par-spd  */
	double recovery,          /**<(I) Recovery rate for par spreads          */
	TDateInterval feeFreq,    /**<(I) Fee frequency for par spreads          */
    TDateInterval integFreq,  /**<(I) Protection leg integration frequency   */
	const TCurve *irCurve,    /**<(I) IR discount curve                      */
	const TCurve *crCurve     /**<(I) Credit curve, or NULL for flat par spd */
);

/**
* Compute forward par CDS spread and annuity value
*
* Return SUCCESS/FAILURE.
*/
int  CrxFwdParCDSSpread(
    double        *parSprd,         /**<(O) par CDS spread                   */
    double        *annuity,         /**<(O) annuity                          */
    TDate         today,            /**<(I) Trade date for new curve         */
	TDate         valDate,          /**<(I) Spot-date for new curve          */
    TDate         cdsStDate,        /**<(I) Start date for CDS               */
    TDate         cdsMatDate,       /**<(I) CDS maturity                     */
    TDateInterval cdsFreq,          /**<(I) CDS fee payment frequency        */
    TDayCountConv cdsDCC,           /**<(I) CDS fee daycount                 */
	KStubLocation dateStub,         /**<(I) CDS fee stub type                */
    KAccrualConv  idxPayAccFlag,    /**<(I) Pay accrued on default or not    */
    KProtPayConv  protPayType,      /**<(I) Pay on default, or at maturity   */
	TDateInterval delay,            /**<(I) Payment delay                    */
    double        recovery,         /**<(I) Recovery rate                    */
	TDateInterval frequency,        /**<(I) Protection leg integration freq. */
    const TCurve  *discZC,          /**<(I) Interest rate zero curve         */
    const TCurve  *sprdZC);         /**<(I) Clean spread zero curve          */


/**
 * Change protection-leg integration frequency of a given credit curve. 
 * Returns a new curve if 
 * successful, or NULL if not. The caller is responsible for freeing the
 * new curve, as well as the old one.
 */
TCurve* CrxCleanCDSCurveConvert(
	TDateInterval newFrequency, /**<(I) Frequency for output curve           */
	TDateInterval frequency,    /**<(I) Frequency of input curve             */
    double        recoveryRate, /**<(I) Assumed recovery for both curves     */
    const TCurve  *discZC,      /**<(I) Interest rate zero curve             */
    const TCurve  *sprdZC       /**<(I) Clean spread zero curve              */
);

#endif    /* _crcrv_H */
