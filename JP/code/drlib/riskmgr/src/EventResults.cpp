
#include "edginc/config.hpp"
#include "edginc/EventResults.hpp"

DRLIB_BEGIN_NAMESPACE

// for reflection 
EventResults::EventResults(): CObject(TYPE), doingLegal(false) {
    riskEvents = EventArray(0);
    legalEvents = EventArray(0);
}

// Invoked when Class is 'loaded'
void EventResults::load(CClassSP& clazz){
    clazz->setPublic(); // make visible to EAS/spreadsheet
    REGISTER(EventResults, clazz);
    SUPERCLASS(CObject);
    EMPTY_SHELL_METHOD(defaultEventResults);
    FIELD(riskEvents, "The array of Event object pointers for risk terms");
    FIELD(legalEvents, "The array of Event object pointers for risk terms");
    FIELD(doingLegal, "Which set of results are we currently filling in");
    FIELD_MAKE_TRANSIENT(doingLegal);
}

IObject* EventResults::defaultEventResults(){
    return new EventResults();
}

CClassConstSP const EventResults::TYPE = CClass::registerClassLoadMethod(
    "EventResults", typeid(EventResults), load);

void EventResults::addEvent(Event* event){
    if (doingLegal){
        legalEvents.append(EventSP(event));
    } else {
        riskEvents.append(EventSP(event));
    }
}

// how many events have we got so far
int EventResults::numEvents() {
    int size = riskEvents.empty() ? 0 : riskEvents.size();
    return size + (legalEvents.empty() ? 0 : legalEvents.size());
}

// flag telling us we're currently retrieving legal terms events
void EventResults::setToLegal() {
    doingLegal = true;
}

// if we only get risk events we need to copy them to legal
// and delete the risk ones
void EventResults::copyRiskToLegal() {
    for (int i = 0; i < riskEvents.size(); ++i) {
        legalEvents.append(riskEvents[i]);
    }
    riskEvents.clear();
}

// if we get risk and legal events we need remove duplicates
void EventResults::removeDuplicateRiskEvents(const DateTime& eventDate) {
    EventArray realRiskEvents = EventArray(0);

    for (int i = 0; i < riskEvents.size(); i++) {
        IEvent* event = dynamic_cast<IEvent*>(riskEvents[i].get());
        if (event->isValidRiskEvent(eventDate)) {
            // it's a good event stash it for later
            realRiskEvents.append(riskEvents[i]);
        }
    }
    riskEvents.clear();
    riskEvents = realRiskEvents;
}

DRLIB_END_NAMESPACE
