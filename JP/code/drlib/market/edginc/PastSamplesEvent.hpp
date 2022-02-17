
#ifndef PAST_SAMPLES_H
#define PAST_SAMPLES_H

#include "edginc/Object.hpp"
#include "edginc/Events.hpp"
#include "edginc/FixingType.hpp"
#include "edginc/ObservationType.hpp"
#include "edginc/ObservationSource.hpp"

DRLIB_BEGIN_NAMESPACE

class IModel;
class CInstrument;
class EventResults;

class MARKET_DLL PastSamplesEvent: public Event, public virtual IEvent{
public:
    static CClassConstSP const TYPE;
    friend class PastSamplesCollector;

    static void load(CClassSP& clazz);

    // interface instrument must satisfy to handle this event
    class MARKET_DLL IEventHandler {
    public:
        static CClassConstSP const TYPE;

        virtual void getEvents(const PastSamplesEvent*  samples,
                               IModel*                  model,
                               const DateTime&          eDate,
                               EventResults*            events) const = 0;
    };

    PastSamplesEvent(const DateTime& eDate, const string source, 
                     const string obsType, const FixingType* context);

    void addSample(const DateTime& sampleDate, const DateTime& obsDate, 
                   double level, bool isOverride);
    
    // once we're done we convert IntArray to a BoolArray
    void buildBoolArray();
        
    CClassConstSP getEventInterface() const;

    virtual void getEvents(IObjectConstSP       obj,
                           IModel*              model,
                           const DateTime&      eDate,
                           EventResults*        events) const;

    // if we get risk and legal events we need to know whether a risk
    // event will just be a duplicated legal event (and hence removable)
    virtual bool isValidRiskEvent(const DateTime& today) const;

private:
    PastSamplesEvent();
    PastSamplesEvent(const PastSamplesEvent &rhs);
    PastSamplesEvent& operator=(const PastSamplesEvent& rhs);
    static IObject* defaultPastSamplesEvent();

    string          observationType;
    string          source;
    DateTimeArraySP samplingDates;   
    DateTimeArraySP observationDates;   
    DoubleArraySP   levels;
    BoolArraySP     overridden;
    FixingTypeSP    context;
    // you can't iterate on a BoolArray so a lame transient field
    IntArraySP      overrideInts; 
};
// support for smart ptrs and arrays
typedef smartPtr<PastSamplesEvent> PastSamplesEventSP;
typedef array<PastSamplesEventSP, PastSamplesEvent> PastSamplesEventArray;
typedef smartPtr<PastSamplesEventArray> PastSamplesEventArraySP;

// a collector class to get samples and put the info together
class MARKET_DLL PastSamplesCollector : public CObject {
public:
    static CClassConstSP const TYPE;

    static void load(CClassSP& clazz);

    // add a single sample
    void addSample(const DateTime&          sampleDate,
                   const DateTime&          obsDate,
                   const ObservationSource* source,
                   const ObservationType*   obsType,
                   double                   level,
                   const FixingType*        fixType,
                   bool                     isOverride);
    
    // put the finished event results in the results structure
    void finaliseResults(EventResults* events);

    PastSamplesCollector(const DateTime& today);
 
private:
    PastSamplesCollector();
    PastSamplesCollector(const PastSamplesCollector &rhs);
    PastSamplesCollector& operator=(const PastSamplesCollector& rhs);
    static IObject* defaultPastSamplesCollector();

    PastSamplesEventArraySP samples;
    int                     numResultSets;
    DateTime                today;
};
typedef smartPtr<PastSamplesCollector> PastSamplesCollectorSP;

DRLIB_END_NAMESPACE

#endif
