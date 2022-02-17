//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : AssetUtil.hpp
//
//   Description : Utility functions for assets
//
//   Author      : Mark A Robson
//
//   Date        : 12 Mar 2001
//
//
//----------------------------------------------------------------------------

#ifndef EDG_ASSET_UTIL_HPP
#define EDG_ASSET_UTIL_HPP
#include "edginc/Asset.hpp"
#include "edginc/DividendCollector.hpp"
#include "edginc/EventAssetMove.hpp"
#include "edginc/Schedule.hpp"
#include "edginc/ObservationType.hpp"
#include "edginc/ObservationSource.hpp"

DRLIB_BEGIN_NAMESPACE

/** Utility functions for assets (all static methods). Separate from CAsset
    for clarity */
class MARKET_DLL AssetUtil{
    friend class ShutTheCompilerUp;
public:
    /** what ccy is an asset in ? */
    static string assetCcy(const CAsset*  asset);

    /** Ensures that weights are positive, that weights.size() =
        numAssets and that the weights sum to 1.0 to within reasonable
        tolerance. The difference between 1.0 and the sum of the
        weights is then distributed back among the weights. */
    static void checkWeights(DoubleArray&       weights, // note not const
                             int                numAssets);
    
    /** Calculates sum(assets[i]->fwdPrice() * weights[i]) */
    static void fwdValue(const CAssetWrapperArray& assets,
                         const CDoubleArray&       weights,
                         const DateTimeArray&      dates,
                         CDoubleArray&             result);

    /** Calculates sum(assets[i]->fwdPrice() * weights[i]) using provided
        algorithm */
    static void fwdValue(const CAssetWrapperArray&         assets,
                         const CDoubleArray&               weights,
                         const DateTimeArray&              dates,
                         const CAsset::FwdValueAlgorithm&  algo,
                         CDoubleArray&                     result);

    /** Calculate the settlement date associated with a given trade date for an
        array of assets. assets must contain an array of at least length 1 */
    static DateTime settleDate(const CAssetWrapperArray& assets,
                               const DateTime&           tradeDate);

    /** checks spotsAtStart.size() = numAssets and each spot at start > DBL_EPSILON 
        if not teies to get it from centralised history*/
    static void validateXCBSpotAtStart(const CAssetWrapperArray& assets,
                                       CDoubleArray&       spotsAtStart,
                                       const DateTime&     startDate,
                                       bool     hasSamplingInfo,
                                       const ObservationSourceArray& sources,
                                       const ObservationTypeArray& obsTypes);

    /** perform some validation for a unit weighted xcb. Checks that
        (i)   numAssets = assets->size() > 0
        (ii)  weights.size() = numAssets
        (iii) all assets != null
        (iv) all weights > 0 */
    static void validateUnitXCB(const CAssetWrapperArray& assets,
                                CDoubleArray&             weights);

    /** perform some validation for a percentage weighted xcb. Checks that
        (i)   numAssets = assets->size() > 0
        (ii)  weights.size() = numAssets
        (iii) all assets != null
        (iv)  all weights > 0 
        (v)   weights add up to one (note weights are 'renormalised' if almost
        equal to 1) */
    static void validatePercXCB(const CAssetWrapperArray& assets,
                                CDoubleArray&             weights);
        
    /** cross validate an asset against instrument.
        Validates:
        (i)  that start and value date are consistent if fwd starting
        (ii)  There are no duplicate assets within the main asset
        (iii) That fx rates are consistent within the assets
        (iv)  That the instrument has not already expired
        (v)   If forward starting, ensure that the option does not start before
              the asset does
        (vi)  If settling into futures, ensures instrument does not mature
              before corresponding future
        (vii) Ensures that the payoff currency matches the asset currency
        (viii)Validates that the value dates in the asset are consistent
        (ix)  Validates that the asset allows an instrument to be fwd starting
              if it is fwd starting
     */
    static void assetCrossValidate(const CAsset*         asset,
                                   bool                  fwdStarting,
                                   const DateTime&       startDate,
                                   const DateTime&       valueDate,
                                   const YieldCurve*     discCcy,
                                   const CInstrument*    instrument);
    //// as above but takes a YieldCurveWrapper
    static void assetCrossValidate(const CAsset*         asset,
                                   bool                  fwdStarting,
                                   const DateTime&       startDate,
                                   const DateTime&       valueDate,
                                   YieldCurveWrapper     discCcy,
                                   const CInstrument*    instrument);

    /** Returns the 'underlying' divs between (startDate, endDate] from
        the supplied asset in the currency of the 'underlying'. Any yield
        divs are left as is, whilst any continuous divs will cause a
        failure. Handles limited types of assets. Only provided for backwards
        compatibility - please use divsBetweenDates() instead */
    static DividendListSP getAllDivsBetweenDates(const CAsset*   asset, 
                                                 const DateTime& start, 
                                                 const DateTime& end);

    //// this one is to be retired
    static DividendListSP getAllDivsBetweenDates(const CAssetWrapper asset, 
                                                 const DateTime& start, 
                                                 const DateTime& end);

    
    /** Returns the divs between (startDate, endDate] from the supplied
        asset in the currency of the asset. Any yield divs are converted
        into dollar divs whilst any continuous divs will cause a
        failure. Handles all types of assets. For greater functionality
        use the DividendCollector class directly. */
    static DividendListSP divsBetweenDates(const CAsset*   asset, 
                                           const DateTime& valDate,
                                           const DateTime& start, 
                                           const DateTime& end);

    /** get only dollar divs from equity asset, works for basket */
    static DividendListSP getDollarDivsBetweenDates(const CAsset*   asset, 
                                                    const DateTime& start, 
                                                    const DateTime& end,
                                                    const bool      convertToDollars = true);

    /** get asset jump events, eturns false if nothing found */
    static bool getJumpEvents(const CAsset*       asset, 
                              const DateTime&     start, 
                              const DateTime&     end, 
                              int                 largestN,
                              EventAssetMove&     result);

    /** get asset jump events, eturns false if nothing found */
    static DividendListSP getDiscreteDivs(const CAsset*       asset,
                                          const DateTime&     valueDate,                        
                                          const DateTime&     start,
                                          const DateTime&     end,
                                          int                 largestN,
                                          DividendCollector::TDividendConversion conversion);

    /** checks if there is any trading time between two dates: (start, end]
    exclusive @start inclusive @end  */
    static bool hasTradingTime(CAssetConstSP asset,
                               const DateTime& start, 
                               const DateTime& end);

    /** checks whether an asset is a simple equity */
    static bool isSimpleEquity(CAssetConstSP asset);

    /** checks if an asset is a simple
        this is actually checking EquityBase */
    static bool isSimpleEquity(const CAsset* asset);

    /** checks whether an asset is a simple, struck or protected equity */
    static bool isSingleEquity(CAssetConstSP asset);

    /** checks whether an asset is a basket */
    static bool isBasket(CAssetConstSP asset);

    static HolidayConstSP getHoliday(const CAsset* asset);

    static MarketObjectSP surfaceSplined(const MarketData* market,
                                        const IModel* model, 
                                        const MarketObjectSP& mo);

    /* given exercise schedule, isAmer flag and underlying, set exercise flag on each time step */
    static void setStepExercise(vector<bool>& stepFlag, 
                                const DateTimeArray& stepDates, 
                                ScheduleConstSP sched, 
                                bool isAmerican, 
                                CAssetConstSP und);

private:
    AssetUtil();
    AssetUtil(const AssetUtil& rhs);
    AssetUtil& operator=(const AssetUtil& rhs);
};

DRLIB_END_NAMESPACE
#endif
