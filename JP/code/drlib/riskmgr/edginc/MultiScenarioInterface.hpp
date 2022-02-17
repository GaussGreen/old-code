//----------------------------------------------------------------------------
//
//   Group       : QR
//
//   Description : Defines interface to Scenario when pricing 
//                 IInstrumentCollection
//
//   Author      : Mark A Robson
//
//   Date        : 10 Feb 2005
//
//
//----------------------------------------------------------------------------

#ifndef MULTISCENARIOINTERFACE_HPP
#define MULTISCENARIOINTERFACE_HPP

#include "edginc/Scenario.hpp"
#include "edginc/RegressionTest.hpp"
#include "edginc/ClientRunnable.hpp"
#include "edginc/RiskMgr.hpp"

DRLIB_BEGIN_NAMESPACE

/** Defines interface to Scenario when pricing IInstrumentCollection */
class RISKMGR_DLL MultiScenarioInterface: public CObject,
                              virtual public IRegressionTest,
                              virtual public ClientRunnable {
public:
    static CClassConstSP const TYPE;

    virtual ~MultiScenarioInterface();

    /** Does not clone supplied parameters */
    MultiScenarioInterface(ScenarioSP              scenario,
                           IModelSP                model, 
                           IInstrumentCollectionSP insts, 
                           CControlSP              ctrl,
                           CMarketDataSP           market);

    /** Runs 'regression test' */
    virtual IObjectSP runTest() const;

    // EdrAction 
    virtual IObjectSP run();

private:
    IObjectSP go() const;
    static void load(CClassSP& clazz);
    static IObject* defaultConstructor();

    /** for addin - runs inputs through risk manager */
    static IObjectSP addinScenario(MultiScenarioInterface* addinParams);

    MultiScenarioInterface();
    MultiScenarioInterface(const MultiScenarioInterface& rhs);
    MultiScenarioInterface& operator=(const MultiScenarioInterface& rhs);
    /// fields
    ScenarioSP              scenario;
    IModelSP                model;
    IInstrumentCollectionSP insts;
    CControlSP              ctrl;
    CMarketDataSP           market;
};

// typedef for smart pointers to MultiScenarioInterface
typedef smartConstPtr<MultiScenarioInterface> MultiScenarioInterfaceConstSP;
typedef smartPtr<MultiScenarioInterface> MultiScenarioInterfaceSP;

DRLIB_END_NAMESPACE

#endif
