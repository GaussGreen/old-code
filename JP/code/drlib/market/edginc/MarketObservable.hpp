//----------------------------------------------------------------------------
//
//   Filename    : MarketObservable.hpp
//
//   Description : A market observable 
//                 This could be a price based asset (equity, commodity, fx)
//                 or a curve/strip of rates for energy/rates
//                 see AssetHistory QLib Design Proposal - Ian Stares Jan 2006
//
//   Author      : Ian Stares
//
//   Date        : January 31 2006
//

//
//----------------------------------------------------------------------------

#ifndef MKTOBSERVABLE_HPP
#define MKTOBSERVABLE_HPP

#include "edginc/config.hpp"
#include "edginc/ObservableHistory.hpp"
#include "edginc/ObservationType.hpp"
#include "edginc/ObservationSource.hpp"
#include "edginc/ObservationOverride.hpp"
#include "edginc/SamplingConvention.hpp"
#include "edginc/FixingType.hpp"
#include "edginc/PastSamplesEvent.hpp"

DRLIB_BEGIN_NAMESPACE

class MARKET_DLL IMarketObservable : virtual public IObject {
public:
    static CClassConstSP const TYPE;

    static const string DEFAULT_SOURCE;

    static ObservationSourceSP getDefaultObsSource();

    // retrieve a single sample
    virtual double pastValue(const DateTime&             sampleDate,
                             const ObservationType*      obsType,
                             const ObservationSource*    source,
                             const FixingType*           fixType,
                             const IObservationOverride* overrides,
                             const SamplingConvention*   sampleRule) const = 0;

    // retrieve single observation date - returns false if obs is to be omitted
    virtual bool observationDate(const DateTime&           sampleDate,
                                 const ObservationSource*  source,
                                 const SamplingConvention* sampleRule,
                                 DateTime*                 obsDate) const = 0;

    // is the given date a holiday for the relevant source
    virtual bool isHoliday(const DateTime&            sampleDate,
                           const ObservationSource*   source) const = 0;

    // adds to events for past samples
    virtual double addPastSampleEvent(const DateTime&             sampleDate,
                                const ObservationType*      obsType,
                                const ObservationSource*    source,
                                const FixingType*           fixType,
                                const IObservationOverride* overrides,
                                const SamplingConvention*   sampleRule,
                                PastSamplesCollector*        collector) const = 0;

    static void load(CClassSP& clazz);
};

DRLIB_END_NAMESPACE

#endif




