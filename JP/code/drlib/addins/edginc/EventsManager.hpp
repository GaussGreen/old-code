
#ifndef EDG_EVENTS_MANAGER_H
#define EDG_EVENTS_MANAGER_H
#include "edginc/RegressionTest.hpp"
#include "edginc/ClientRunnable.hpp"
#include "edginc/MarketData.hpp"
#include "edginc/Instrument.hpp"
#include "edginc/Model.hpp"

DRLIB_BEGIN_NAMESPACE

class OutputRequestCalculator;

/** Directs retrieval of MATRIX events */
class ADDINS_DLL EventsManager : public CObject,
                      virtual public IRegressionTest,
                      virtual public ClientRunnable {
public:  
    static CClassConstSP const TYPE;
    friend class EventsManagerHelper;

    /** Runs 'regression test' for model, inst and control in this class */
    virtual  IObjectSP runTest() const;

    // EdrAction 
    virtual IObjectSP run();

    IObjectSP run() const;

private:
    // model and inst may be DRWrappers
    IObjectSP model;
    IObjectSP inst;

    CMarketDataSP market;
    string        mode;   // do everything, selection or exclude
    StringArraySP filter; // specific events to do/ignore
    
    EventsManager();
    EventsManager(const EventsManager& rhs);
    EventsManager& operator=(const EventsManager& rhs);

    /** for addin - runs inputs through events manager */
    static IObjectSP addinEventsMgr(EventsManager* addinParams);
};

DRLIB_END_NAMESPACE

#endif
