//----------------------------------------------------------------------------
//
//   Filename    : ObservationOverride.hpp
//
//   Description : Interface for classes which provide overrides for
//                 observations for IMarketObservables 
//                 e.g. PastValues provides such overrides
//
//   Author      : Ian Stares
//
//   Date        : July 4 2006
//
//
//----------------------------------------------------------------------------

#ifndef OBSOVERRIDE_HPP
#define OBSOVERRIDE_HPP

#include "edginc/config.hpp"
#include "edginc/DateTime.hpp"
#include "edginc/ObservationType.hpp"
#include "edginc/PastSamplesEvent.hpp"

DRLIB_BEGIN_NAMESPACE

class MARKET_DLL IObservationOverride : public virtual IObject {
public:
    static CClassConstSP const TYPE;

    // retrieve an override for a given sample date
    // note sample dates are dates we ATTEMPT to sample
    // returns true if override is present and if so returns result in value
    virtual bool value(const DateTime&           sampleDate,
                       const ObservationType*    obsType,
                       double* value) const = 0;

    // add a past sample to an event collector
    // returns false if it fails to find an override
    // if it succeeds the reuslt is stuck in level
    virtual bool addPastSampleEvent(const DateTime&          sampleDate,
                                    const ObservationType*   obsType,
                                    const ObservationSource* source,
                                    const FixingType*        fixType,
                                    PastSamplesCollector*    collector,
                                    double&                  level) const = 0;
private:
    static void load(CClassSP& clazz);
};

typedef smartPtr<IObservationOverride> IObservationOverrideSP;
typedef array<IObservationOverrideSP, IObservationOverride> IObservationOverrideArray;

DRLIB_END_NAMESPACE

#endif




