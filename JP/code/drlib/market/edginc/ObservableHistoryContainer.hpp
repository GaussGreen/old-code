//----------------------------------------------------------------------------
//
//   Filename    : ObservableHistoryContainer.hpp
//
//   Description : Interface between MarketObservable and ObservableHistory
//                 stores all histories for a MarketObservable and
//                 handles look up based on source before delegating
//                 down to ObservableHistory
//                 see AssetHistory QLib Design Proposal - Ian Stares Jan 2006
//
//   Author      : Ian Stares
//
//   Date        : January 31 2006
//
//
//----------------------------------------------------------------------------

#ifndef OBSHISTORYCONT_HPP
#define OBSHISTORYCONT_HPP

#include "edginc/config.hpp"
#include "edginc/ObservableHistory.hpp"
#include "edginc/ObservationType.hpp"
#include "edginc/SamplingConvention.hpp"
#include "edginc/FixingType.hpp"

DRLIB_BEGIN_NAMESPACE

class MARKET_DLL IObservableHistoryContainer {
public:
    // retrieve a single sample
    virtual double pastValue(const DateTime&            sampleDate,
                             const ObservationType*     obsType,
                             const string&              source,
                             const FixingType*          fixType,
                             const SamplingConvention*  sampleRule) const = 0;

    // retrieve a single observation date
    // Returns false if obs is to be omitted
    virtual bool observationDate(const DateTime&           sampleDate,
                                const string&             source,
                                const SamplingConvention* sampleRule,
                                DateTime*                 obsDate) const = 0;

    // is the given date a holiday for the relevant source
    virtual bool isHoliday(const DateTime& sampleDate,
                           const string&   source) const = 0;
};

DRLIB_END_NAMESPACE

#endif




