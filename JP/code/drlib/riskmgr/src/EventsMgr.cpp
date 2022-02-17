
#include "edginc/config.hpp"
#include "edginc/EventsMgr.hpp"
#include "edginc/LegalTerms.hpp"
#include "edginc/Modifier.hpp"
#include "edginc/OutputRequest.hpp"
#include "edginc/Holiday.hpp"
#include "edginc/Results.hpp"
#include "edginc/ObjectIteration.hpp"

DRLIB_BEGIN_NAMESPACE

string EventsMgr::SELECTION_MODE = "SELECTION";
string EventsMgr::EXCLUSION_MODE = "EXCLUSION";
string EventsMgr::COMPLETE_MODE  = "COMPLETE";


EventResults* EventsMgr::run(const IModel*      model,
                             const CInstrument* inst,
                             const MarketData*  market,
                             const string&      mode,
                             const StringArray* filter){
    static const string method("EventsMgr::run");

    try {
        // work on a copy of the instrument just to be safe
        CInstrumentSP instrument(copy(inst));
        // copy model as well since we might be filling it with market data
        IModelSP mdl(copy(model));

        mdl->getInstrumentAndModelMarket(market, instrument.get());

        // better do this before we do anything
        instrument->Validate();

        // roll instrument to now so that past closings etc get updated
        // very lame. Might need rethink when we do asset history properly
        HolidaySP holidays(Holiday::noHolidays());
        ThetaSP thetaShift(new Theta(0, holidays));
        thetaShift->applyScenario(instrument);

        // now go and round up those events!
        EventResultsSP results(getEvents(mdl.get(), instrument.get(), 
                                         market->GetReferenceDate(),
                                         mode,
                                         filter));
        return results.release();
    }
    catch (exception &e) {
        throw ModelException(e, method, "event retrieval failed");
    }
}

//// EVENTS output request is not handled by the instrument
//// but at the top level through the events mgr interface. 
//// This calculator class sits inside the OutputRequestHelper
class EventsMgr::Calculator : public OutputRequestCalculator{
public:
    void calculate(OutputRequest*     request, 
                   const IModel*      model, 
                   const CInstrument* instrument,
                   Results*           results) {
        static const string method("EventsMgr::Calculator::calculate");
        try {        
            // work on a copy of the instrument just to be safe
            CInstrumentSP inst(copy(instrument));
            // copy model as well 
            IModelSP      mdl(copy(model));

            // better do this before we do anything
            inst->Validate();

            EventResultsSP events(getEvents(mdl.get(), inst.get(), 
                                            inst->getValueDate(),
                                            COMPLETE_MODE,
                                            0));

            if (events->numEvents() > 0) {
                results->storeRequestResult(request, events);
            } else {
                results->storeNotApplicable(request);
            }
        }
        catch (exception &e) {
            throw ModelException(e, method, "event retrieval failed");
        }
    }
    Calculator() {};
};

/** Creates an OutputRequestCalculator which can be used for a default
    implementation of the Events output request */
OutputRequestCalculator* EventsMgr::createCalculator(){
    return new Calculator();
}

// Support for getting events from an invidual object. The invoke method below
// is called for each object of the right type.
class EventsMgr::GetEventAction : virtual public ObjectIteration::IActionConst{
public:
    ////// fields /////
    const DateTime&      eventDate;
    IModel*              model;
    EventResults*        events;
    IEvent*              eventType;

    // Constructor - builds the EventResults structure
    GetEventAction(const DateTime& eDate, IModel* mdl, 
                   IEvent* eventType, EventResults* events): 
        eventDate(eDate), model(mdl), events(events), eventType(eventType) {}

    // method invoked by recurse routine to get events
    bool invoke(const ObjectIteration::State& state, IObjectConstSP obj){
        // call handler method on the event type which 
        // will delegate to the implemented interface
        eventType->getEvents(obj, model, eventDate, events);
        return true;
    }
};

// gets one set of events - either risk or legal
// note the flag for which one we're doing lives in the EventResults
// and is set outside of this function
void EventsMgr::getEventsSet(IModel*             model,
                             CInstrument*        inst,
                             const DateTime&     eventDate,
                             const string&       mode,
                             const StringArray*  filter,
                             EventResults*       events) {
    static const string method("EventsMgr::getEventsSet");
    try {
        CClassVec eventIFaces;

        if (mode != COMPLETE_MODE && !filter) {
            throw ModelException(method, 
                                 "Running in " + mode + 
                                 " mode but no filter supplied");
        }
            
        if (mode == COMPLETE_MODE) {
            // do everything possible
            eventIFaces = IEvent::TYPE->getAssignableClasses();
        }
        else if (mode == SELECTION_MODE) {
            // only do what was asked for
            for (int i = 0; i < filter->size(); i++) {
                eventIFaces.push_back(CClass::forName((*filter)[i]));
            }
        }
        else if (mode == EXCLUSION_MODE) {
            // don't do what's in the filter
            // this is a weak way of doing the filter - rubbish performance
            CClassVec allEvents = IEvent::TYPE->getAssignableClasses();
            CClassVec dontDoMe;
            int i;
            for (i = 0; i < filter->size(); i++) {
                dontDoMe.push_back(CClass::forName((*filter)[i]));
            }
            for (i = 0; i < (int)allEvents.size(); i++) {
                bool skip = false;
                for (unsigned int j = 0; !skip && j < dontDoMe.size(); j++) {
                    skip = (allEvents[i] == dontDoMe[j]);
                }
                if (!skip) {
                    eventIFaces.push_back(allEvents[i]);
                }
            }           
        } 
        else {
            throw ModelException(method, "Unknown filter mode " + mode);
        }

        // create a map of names to classes and a sorted vector of names
        map<string, CClassConstSP> classMap;
        vector<string> classNames;
        for(unsigned int i=0; i<eventIFaces.size(); i++) {
            CClassConstSP iface = eventIFaces[i];
            classNames.push_back(iface->getName());
            classMap[iface->getName()] = iface;
        }
        sort(classNames.begin(), classNames.end());

        // process the interfaces in alphabetical order
        for (vector<string>::iterator i = classNames.begin();
                                      i != classNames.end();
                                      i++ ){
            CClassConstSP eventType = classMap[*i];
            if (!Modifier::isAbstract(eventType->getModifiers())){
                IObjectSP eventObj(eventType->newInstance());
                IEvent* eInst = dynamic_cast<IEvent*>(eventObj.get());
                CClassConstSP eventHandler = eInst->getEventInterface();
                
                if (ObjectIteration::dependsUpon(IObjectConstSP(inst),
                                                 eventHandler, false)) {
                    // copy inst and model every time through        
                    IModelSP modelCopy(copy(model));
                    CInstrumentSP instCopy(copy(inst));

                    // create our Action class
                    GetEventAction action(eventDate, modelCopy.get(), 
                                          eInst, events);
                    // create the ObjectIteration class
                    ObjectIteration iteration(eventHandler);
                    // then recurse over all the relevant objects
                    iteration.recurse(action, instCopy);
                }
            }
        }
    } catch (exception& e) {
        throw ModelException(e, method);
    }
}

// recurse through and get events
EventResults* EventsMgr::getEvents(IModel*            model,
                                   CInstrument*       inst,
                                   const DateTime&    eventDate,
                                   const string&      mode,
                                   const StringArray* filter){
    try{
        EventResultsSP events(new EventResults());

        //get the 'Risk' events
        getEventsSet(model, inst, eventDate, mode, filter, events.get());

        // and repeat for Legal Terms if necessary
        IScenarioShiftSP legalTermsScenario(new LegalTerms());
        CInstrumentSP instLT(copy(inst));
        // we don't want to look for events again if there are no legal terms
        if (legalTermsScenario->applyScenario(instLT)){
            // switch flag so events are put in the legal bucket
            events->setToLegal();

            getEventsSet(model, instLT.get(), eventDate, 
                         mode, filter,
                         events.get());
            
            // also given that for dates in the past legal/risk coincide
            // we need to remove any spurious risk events
            events->removeDuplicateRiskEvents(eventDate);                
        } else {
            // no legal terms so we only found one set of events
            // we should move them to the legal result set
            // as, in this case, that is where they belong
            events->copyRiskToLegal();                
        }

        return events.release();
    } catch (exception& e){
        throw ModelException(e, "EventsMgr::getEvents");
    }
}

DRLIB_END_NAMESPACE

