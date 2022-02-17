
#include "edginc/config.hpp"
#include "edginc/EventsManager.hpp"
#include "edginc/EventsMgr.hpp"
#include "edginc/Addin.hpp"
#include "edginc/DRWrapper.hpp"

DRLIB_BEGIN_NAMESPACE

// for reflection
EventsManager::EventsManager() : 
    CObject(TYPE), model(0), inst(0), market(0), mode(EventsMgr::COMPLETE_MODE) {
    // empty
}

/** Runs 'regression test' for model and inst in this class */
IObjectSP EventsManager::runTest() const{
    return run();
}

// EdrAction 
IObjectSP EventsManager::run(){
    const EventsManager* emr = this; // getting into a mess with names
    return emr->run();
}

IObjectSP EventsManager::run() const{
    static const string method("EventsManager::run");

    try {
        if (!market){
            throw ModelException(method, "NULL market");
        }
        if (!model){
            throw ModelException(method, "NULL model");
        }
        if (!inst){
            throw ModelException(method, "NULL instrument");
        }

        if (mode != EventsMgr::COMPLETE_MODE && !filter.get()) {
            throw ModelException(method, "Running in " + mode + 
                                 " mode but no filter supplied");
        }
            
        // first convert the inst/model from drwrapper if necessary
        IModelSP      modelSP(IModelSP::dynamicCast(
            DRWrapper::drWrapperToObject(model, IModel::TYPE)));
        CInstrumentSP instSP(CInstrumentSP::dynamicCast(
            DRWrapper::drWrapperToObject(inst, CInstrument::TYPE)));

        // now go and round up those events!
        EventResultsSP events(EventsMgr::run(modelSP.get(), instSP.get(), 
                                             market.get(),
                                             mode,
                                             filter.get()));
        return events;

    }
    catch (exception &e) {
        throw ModelException(e, method, "event retrieval failed");
    }
}

/** addin function wrapper for EventsManager */
static IObjectSP addinEventsMgr(EventsManager* params){
    return params->run();
}

class EventsManagerHelper{
public:
    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        clazz->setPublic(); // make visible to EAS/spreadsheet
        REGISTER(EventsManager, clazz);
        SUPERCLASS(CObject);
        IMPLEMENTS(ClientRunnable);
        EMPTY_SHELL_METHOD(defaultInterface);
        FIELD(model, "model");
        FIELD(inst, "instrument");
        FIELD(market, "market");
        FIELD(mode, "mode: filter events or not");
        FIELD_MAKE_OPTIONAL(mode);
        FIELD(filter, "filter");
        FIELD_MAKE_OPTIONAL(filter);

        // registration for addin function
        Addin::registerClassObjectMethod("EventsManager",
                                         Addin::RISK,
                                         "Retrieves events",
                                         EventsManager::TYPE,
                                         true,
                                         Addin::returnHandle,
                                         (Addin::ObjMethod*)addinEventsMgr);
        
    }

    static IObject* defaultInterface(){
        return new EventsManager();
    }
};

CClassConstSP const EventsManager::TYPE = CClass::registerClassLoadMethod(
    "EventsManager", typeid(EventsManager), EventsManagerHelper::load); 

DRLIB_END_NAMESPACE

