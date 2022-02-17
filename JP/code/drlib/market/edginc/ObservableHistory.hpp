//----------------------------------------------------------------------------
//
//   Filename    : ObservableHistory.hpp
//
//   Description : Holds all the samples for a MarketObservable
//                 This could be a list of prices or a list of curves
//                 see AssetHistory QLib Design Proposal - Ian Stares Jan 2006
//
//   Author      : Ian Stares
//
//   Date        : January 31 2006
//
//
//----------------------------------------------------------------------------

#ifndef OBSHISTORY_HPP
#define OBSHISTORY_HPP

#include "edginc/config.hpp"
#include "edginc/ObservationType.hpp"
#include "edginc/SamplingConvention.hpp"
#include "edginc/FixingType.hpp"

DRLIB_BEGIN_NAMESPACE

// interface that observable history needs
// so ability to retrieve a single sample (possibly with context and
// adjusting ro holidays, and ability to retrieve actual sample dates
// and judge whether a given date is a holiday
class MARKET_DLL IObservableHistory {
public:
    // retrieve a single sample
    virtual double pastValue(const DateTime&           sampleDate,
                             const ObservationType*    obsType,
                             const FixingType*         fixType,
                             const SamplingConvention* sampleRule) const = 0;

    // retrieve a single observation date
    // Returns false if obs is to be omitted
    virtual bool observationDate(const DateTime&           sampleDate,
                                 const SamplingConvention* sampleRule,
                                 DateTime*                 obsDate) const = 0;

    // is the given date a holiday for the relevant source
    virtual bool isHoliday(const DateTime& sampleDate) const = 0;

    // returns the name of the market observable for which this is
    // the history
    virtual string getMarketObservableName() const = 0;

    // get the name of the source of this history
    virtual string getSource() const = 0;
};

DRLIB_END_NAMESPACE

#endif




