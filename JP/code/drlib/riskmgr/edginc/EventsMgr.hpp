
#ifndef EDG_EVENTS_MGR_H
#define EDG_EVENTS_MGR_H
#include "edginc/Instrument.hpp"
#include "edginc/Model.hpp"
#include "edginc/EventResults.hpp"
#include "edginc/MarketData.hpp"

DRLIB_BEGIN_NAMESPACE

class OutputRequestCalculator;

/** Directs retrieval of MATRIX events */
class RISKMGR_DLL EventsMgr {
public:  
    static CClassConstSP const TYPE;
    
    // generate events for given imnt/model/market
    // mode drives the interpretation of the filter
    // see static strings below
    static EventResults* run(const IModel*      model,
                             const CInstrument* inst,
                             const MarketData*  market,
                             const string&      mode,
                             const StringArray* filter);


    /** Creates an OutputRequestCalculator which can be used for a default
        implementation of the Events output request */
    static OutputRequestCalculator* createCalculator();

    // recurse through and get events
    // market data already extracted
    static EventResults* getEvents(IModel*            model,
                                   CInstrument*       inst,
                                   const DateTime&    eventDate,
                                   const string&      mode,
                                   const StringArray* filter);

    // how to drive events:
    // 1 - only do what was asked for
    // 2 - everything but what was asked for
    // 3 - everything (default)
    static string SELECTION_MODE;
    static string EXCLUSION_MODE;
    static string COMPLETE_MODE;

private:
    class Calculator;
    class GetEventAction;

    // gets one set of events - either risk or legal
    static void getEventsSet(IModel*             model,
                             CInstrument*        inst,
                             const DateTime&     eventDate,
                             const string&       mode,
                             const StringArray*  filter,
                             EventResults*       events);
    EventsMgr();
    EventsMgr(const EventsMgr& rhs);
    EventsMgr& operator=(const EventsMgr& rhs);
};

DRLIB_END_NAMESPACE

#endif
