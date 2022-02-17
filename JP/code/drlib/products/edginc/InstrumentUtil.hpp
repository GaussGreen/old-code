//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : InstrumentUtil.hpp
//
//   Description : Utility functions for instruments
//
//   Author      : André Segger
//
//   Date        : 04 Jun 2001
//
//
//----------------------------------------------------------------------------

#ifndef EDG_IMNT_UTIL_HPP
#define EDG_IMNT_UTIL_HPP

#include "edginc/Control.hpp"
#include "edginc/Asset.hpp"
#include "edginc/InstrumentSettlement.hpp"

DRLIB_BEGIN_NAMESPACE

/** Utility functions for Instruments (all static methods). Separate from 
    CInstrument for clarity */
class PRODUCTS_DLL InstrumentUtil {
public:

    /** outputs all forwards at maturity */
    static void recordFwdAtMat(
                      Control*         control,
                      Results*         results,
                      const DateTime&  matDate,
                      const DateTime&  valueDate,
                      const CAsset*    asset);


    /** outputs indicative volatility */
    static void recordIndicativeVol(
                    Control*     control,
                    Results*     results,
                    double       indVol);

    /** outputs delay Price*/
    static void delayPriceHelper(Control*                    control,
                                 Results*                    results,
                                 const double&               fairValue,
                                 const DateTime&             valueDate,
                                 const YieldCurve*           discount,
                                 const Asset*                asset,
                                 const InstrumentSettlement* premiumSettlement);

    static void recordDelayPrice(Control*     control,
                                 Results*     results,
                                 double       delayPrice);

    static double calculateDelayPrice(const double&     fairValue,
                                      const DateTime&   valueDate,
                                      const DateTime&   payDate,
                                      const YieldCurve* discount);

    static double scalePremium(bool          oneContract,
                               bool          fwdStarting,
                               const double& notional,
                               const double& fwdAtStart,
                               const double& initialSpot);

    static void recordPricePCTpayoff(Control*     control,
                                     Results*     results,
                                     double       pricePCTpayoff);

    static void recordMaxPayoff(Control*     control,
                                Results*     results,
                                double       maxPayoff);

    /** Sorts the supplied performances and then takes the weight sum of them
        using supplied weights. weights.size() must equal perf.size(). The
        sorted performances are returned in perf.  */
    static double rainbowPerformance(const DoubleArray& weights,
                                     DoubleArray&       perf);

    /** Takes the weighted sum of the supplied performances.
        weights.size() must equal perf.size().  */
    static double weightedPerformance(const DoubleArray& weights,
                                      const DoubleArray& perf);

};

DRLIB_END_NAMESPACE
#endif
