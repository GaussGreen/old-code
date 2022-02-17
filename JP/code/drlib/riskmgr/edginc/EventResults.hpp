
#ifndef EDG_EVENT_RESULTS_H
#define EDG_EVENT_RESULTS_H

#include "edginc/Events.hpp"
#include "edginc/Object.hpp"

DRLIB_BEGIN_NAMESPACE

// Stores MATRIX event results
class RISKMGR_DLL EventResults: public CObject{
public:
    static CClassConstSP const TYPE;

    /** Stores a generic event */ 
    void addEvent(Event* event);

    static void load(CClassSP& clazz);

    // how many events have we got so far
    int numEvents();

    // flag telling us we're currently retrieveing legal terms events
    void setToLegal();

    // if we only get risk events we need to copy them to legal
    // and delete the risk ones
    void copyRiskToLegal();

    // if we get risk and legal events we need remove duplicates
    void removeDuplicateRiskEvents(const DateTime& eventDate);

    /** Creates an empty results set */
    EventResults();

private:
    static IObject* defaultEventResults();

    EventArray  riskEvents;
    EventArray  legalEvents;
    bool        doingLegal;
};

typedef smartPtr<EventResults> EventResultsSP;

DRLIB_END_NAMESPACE

#endif
