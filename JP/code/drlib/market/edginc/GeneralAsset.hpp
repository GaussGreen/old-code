//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : GeneralAsset.hpp
//
//   Description : GeneralAsset interface
//
//   Author      : Mark A Robson
//
//   Date        : 14 May 2004
//
//
//----------------------------------------------------------------------------

#ifndef EDR_GENERALASSET_HPP
#define EDR_GENERALASSET_HPP
#include "edginc/MarketFactor.hpp"
#include "edginc/DateTime.hpp"
#include "edginc/AtomicArray.hpp"

DRLIB_BEGIN_NAMESPACE
class IVolProcessed;
class CVolRequest;
class Theta;
class OutputRequest;
class Results;

/** A GeneralAsset covers equity/fx type assets as well as IR assets */
class MARKET_DLL IGeneralAsset: public virtual IMarketFactor{
public:
    static CClassConstSP const TYPE; // in MarketFactor.cpp

    IGeneralAsset();  // in MarketFactor.cpp

    virtual ~IGeneralAsset();  // in MarketFactor.cpp

    /** returns the asset name without the effect of any current treatment 
        ie the name that would have been obtained if currency treatment had 
        been CCY_TREATMENT_VANILLA. Default implementation returns getName() */
    virtual string getTrueName() const = 0;

    /** returns the spot price */
    virtual double getSpot() const = 0;

    /** Returns an processed vol - which combines the vol market data with the
        instrument data in the volRequest */
    virtual IVolProcessed* getProcessedVol(
        const CVolRequest* volRequest) const = 0;

    /** Calculates the expected spot price of the asset at the given date.
        Do not use this repeatedly to calculate values over a set of
        dates (poor performance) - instead use other fwdValue method */
    virtual double fwdValue(const DateTime& date) const = 0;

    /** Quicker version of above when you have multiple dates */
    virtual void fwdValue(const DateTimeArray& dateList,
                          CDoubleArray&        result) const = 0;

    /** Returns the name (not the ISO code) of the asset ccy */
    virtual string getYCName() const = 0;

    /** Calculates the expected spot price of the asset at the given date if
        the spot price had the given value spot on spotDate */
    virtual double fwdFwd(const DateTime& spotDate,
                          double          spot, 
                          const DateTime& fwdDate) const = 0;

    /** Array version of fwdFwd */
    virtual void fwdFwd(const DateTime&      spotDate,
                        double               spot, 
                        const DateTimeArray& fwdDates,
                        DoubleArray&         results) const = 0;

    /** record forwards at maturity*/
    virtual void recordFwdAtMat(OutputRequest*  request,
                                Results*        results,
                                const DateTime& maturityDate) const = 0;

    /** Returns the spot value to use when populating a sample during a
        theta shift */
    virtual double getThetaSpotOnDate(const Theta *shift,
                                      const DateTime &date) const = 0;

   /** Calculate the settlement date associated with a given trade date */
    virtual DateTime settleDate(const DateTime& tradeDate) const = 0;


private:
    static void load(CClassSP& clazz);
};

DRLIB_END_NAMESPACE
#endif
