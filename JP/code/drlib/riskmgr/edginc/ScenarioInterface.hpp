//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : ScenarioInterface.hpp
//
//   Description : Defines interface to Scenario
//
//   Author      : Andrew J Swain
//
//   Date        : 25 April 2001
//
//
//----------------------------------------------------------------------------

#ifndef SCENARIOINTERFACE_HPP
#define SCENARIOINTERFACE_HPP

#include "edginc/Scenario.hpp"
#include "edginc/RegressionTest.hpp"
#include "edginc/ClientRunnable.hpp"
#include "edginc/RiskMgr.hpp"

DRLIB_BEGIN_NAMESPACE

class RISKMGR_DLL ScenarioInterface: public CObject,
                         virtual public IRegressionTest,
                         virtual public ClientRunnable {
public:
    static CClassConstSP const TYPE;

    virtual ~ScenarioInterface();

    /// fields
    ScenarioSP    scenario;
    IModelSP      model;
    CInstrumentSP inst;
    CControlSP    ctrl;
    CMarketDataSP market;

    /** Does not clone supplied parameters */
    ScenarioInterface(ScenarioSP    scenario,
                      IModelSP      model, 
                      CInstrumentSP inst, 
                      CControlSP ctrl,
                      CMarketDataSP market);

    /** Runs 'regression test' */
    virtual IObjectSP runTest() const;

    // EdrAction 
    virtual IObjectSP run();

private:
    friend class ScenarioInterfaceHelper;
    IObjectSP go() const;
    static void load(CClassSP& clazz);
    static IObject* defaultInterface();
    
    ScenarioInterface();
    ScenarioInterface(const ScenarioInterface& rhs);
    ScenarioInterface& operator=(const ScenarioInterface& rhs);
};

// typedef for smart pointers to ScenarioInterface
typedef smartConstPtr<ScenarioInterface> ScenarioInterfaceConstSP;
typedef smartPtr<ScenarioInterface> ScenarioInterfaceSP;

DRLIB_END_NAMESPACE

#endif
