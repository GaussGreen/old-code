//----------------------------------------------------------------------------
//
//   Group       : Derivatives Research
//
//   Filename    : CCMPriceUtil.hpp
//
//   Description : CCM Price Utilities 
//
//   Date        : Aug 2004
//
//
//----------------------------------------------------------------------------

#ifndef EDRCCMPRICEUTIL_HPP
#define EDRCCMPRICEUTIL_HPP

#include "edginc/CashFlow.hpp"
#include "edginc/IDiscountCurveRisky.hpp"
DRLIB_BEGIN_NAMESPACE

class MARKET_DLL CCMPriceUtil{
public:
    enum ExpLossType{
        FLATFWD = 0,
        LINEAR
    };


/**
 * Calculate the PV of a accrual on default leg with delay
 *     PV(T_s,T_e,d)  = \int_{T_s}^{T_e} {(a+bt) Z_r(t+d) dZ_\Lambda(t)}
 * With the following notations
 * Z_r(t)       & = & \exp(-\int_0^t{r_s ds})
 * Z_\Lambda(t) & = & \exp(-\int_0^t{\lambda_s ds})
 *
 * The interval [0,T] is partitioned into I_i=[t_i, t_{i+1}].
 * An interpolation method is chosen within I_i:
 * - riskless rate r_t are piecewise constant (flat forward)
 * - either \lambda_t is piecewise constant (FLATFWD expected loss)
 *   or Z_\Lambda(t) is piecewise linear (LINEAR expected loss)
 */
    static double contingentAndAccrPrice(
        double             Ts,    /* (I) starting time                       */
        double             Te,    /* (I) end time                            */
        double             d,     /* (I) payment delay                       */
        const DoubleArray& tR,    /* (I) timeline for Zr                     */
        const DoubleArray& ZR,    /* (I) discount factor                     */
        const DoubleArray& tL,    /* (I) timeline for Zl                     */
        const DoubleArray& ZL,    /* (I) expected survival                   */
        double             a,     /* (I) ctg leg factor                      */
        double             b,     /* (I) acc payoff factor                   */
        ExpLossType        type); /* (I) FLATFWD or LINEAR                   */

    /**
     * Calculate the PV of a contingent leg with delay
     *     ctg(T_s,T_e,d)  = \int_{T_s}^{T_e} {Z_r(t+d) dZ_\Lambda(t)}
     * With the following notations
     * Z_r(t)       & = & \exp(-\int_0^t{r_s ds})
     * Z_\Lambda(t) & = & \exp(-\int_0^t{\lambda_s ds})
     *
     * The interval [0,T] is partitioned into I_i=[t_i, t_{i+1}].
     * An interpolation method is chosen within I_i:
     * - riskless rate r_t are piecewise constant (flat forward)
     * - either \lambda_t is piecewise constant (FLATFWD expected loss)
     *   or Z_\Lambda(t) is piecewise linear (LINEAR expected loss)
     */
    static double contingentPrice( 
        double             Ts,     /* (I) ctg leg begin                      */
        double             Te,     /* (I) ctg leg end                        */
        double             d,      /* (I) payment delay after default        */
        const DoubleArray& tR,     /* (I) timeline for Zr                    */
        const DoubleArray& ZR,     /* (I) discount factor                    */
        const DoubleArray& tL,     /* (I) timeline for Zl                    */
        const DoubleArray& ZL,     /* (I) expected survival                  */
        const ExpLossType  type);  /* (I) FLATFWD or LINEAR */

    /**
     * Calculate the PV of a accrual on default leg with delay
     *     acc(T_s,T_e,d)  = \int_{T_s}^{T_e} {t Z_r(t+d) dZ_\Lambda(t)}
     * With the following notations
     * Z_r(t)       & = & \exp(-\int_0^t{r_s ds})
     * Z_\Lambda(t) & = & \exp(-\int_0^t{\lambda_s ds})
     *
     * The interval [0,T] is partitioned into I_i=[t_i, t_{i+1}].
     * An interpolation method is chosen within I_i:
     * - riskless rate r_t are piecewise constant (flat forward)
     * - either \lambda_t is piecewise constant (FLATFWD expected loss)
     *   or Z_\Lambda(t) is piecewise linear (LINEAR expected loss)
     */
    static double accrualOnDefault( 
        double             Ts,     /* (I) ctg leg begin                   */
        double             Te,     /* (I) ctg leg end                     */
        double             d,      /* (I) payment delay after default     */
        const DoubleArray& tR,     /* (I) timeline for Zr                 */
        const DoubleArray& ZR,     /* (I) discount factor                 */
        const DoubleArray& tL,     /* (I) timeline for Zl                 */
        const DoubleArray& ZL,     /* (I) expected survival               */
        ExpLossType        type);  /* (I) FLATFWD or LINEAR*/

    /**
     * Calculate the expected average notional
     *     avgNtl(T_s,T_e)  = 
     \frac{1}{T_e-T_s} \int_{T_s}^{T_e} {t Z_\Lambda(t) dt}
     * With the following notations
     * Z_\Lambda(t) & = & \exp(-\int_0^t{\lambda_s ds})
     *
     * The interval [0,T] is partitioned into I_i=[t_i, t_{i+1}].
     * An interpolation method is chosen within I_i:
     * - either \lambda_t is piecewise constant (FLATFWD expected loss)
     *   or Z_\Lambda(t) is piecewise linear (LINEAR expected loss)
     */
    static double averageNotional( 
        double             Ts,     /* (I) ctg leg begin                     */
        double             Te,     /* (I) ctg leg end                       */
        const DoubleArray& tL,     /* (I) timeline for Zl                   */
        const DoubleArray& ZL,     /* (I) expected survival                 */
        ExpLossType        type);  /* (I) FLATFWD or LINEAR*/


    /** Linear Loss interpolation  */
    static double linearLossInterpolate(
        const DateTime&      maturity,
        const DateTimeArray& timeline,
        const DoubleArray&   expectedTrancheNotional);

    /* Returns the expected losses (in the fee leg) between today and maturity, 
       using the effectiveCurve */
    static double expectedFeeNotionalLoss(
        const double                outstandingNotional, /* initialTrancheSize -
                                                            pastTrancheLoss */
        const DateTime&             today,               /* reference date */
        const DateTime&             maturity,            // point to be returned
        const IDiscountCurveRiskySP effectiveCurve);     // effective curve to

    /* Factor out past outstanding from above method */
    static double outstandingNotionalAtMaturity(
        const double         initialSize,
        const DateTime&      maturity,            // point to be returned
        const CashFlowArray& pastTrancheLosses);

    static double pastLossesInPeriod(
        const CashFlowArray& pastTrancheLosses,
        const DateTime&      periodStart,         // point to be returned
        const DateTime&      periodEnd,
        const BoolArray&     payPastTrancheLosses);

    static double expectedLossesInPeriod(
        const double                outstandingNotional, /* initialTrancheSize -
                                                            pastTrancheLoss */
        const DateTime&             periodStart,
        const DateTime&             periodEnd,
        const IDiscountCurveRiskySP effectiveCurve);
};

DRLIB_END_NAMESPACE
#endif

