//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : Scenario.hpp
//
//   Description : Defines a collection of scenario shifts
//
//   Author      : Andrew J Swain
//
//   Date        : 24 April 2001
//
//
//----------------------------------------------------------------------------


#ifndef SCENARIO_HPP
#define SCENARIO_HPP
#include <string>
#include "edginc/Object.hpp"
#include "edginc/ScenarioShift.hpp"
#include "edginc/Control.hpp"
#include "edginc/Instrument.hpp"
#include "edginc/Model.hpp"
#include "edginc/Results.hpp"
#include "edginc/RiskMgr.hpp"

DRLIB_BEGIN_NAMESPACE
/** Defines a collection of scenario shifts */
class RISKMGR_DLL Scenario: public CObject{
public:
    static CClassConstSP const TYPE;

    /** Takes copy of shift array */
    Scenario(IScenarioShiftArrayConstSP shifts);

    /** Apply this scenario to supplied instrument/market */
    void apply(IModelSP                model, // (M)
               IInstrumentCollectionSP instruments) const; // (M)

    /** Apply this scenario to supplied instrument/market, before market data is retrieved */
    void preapply(IModelSP                model, // (M)
                  IInstrumentCollectionSP instruments) const; // (M)

    /** Writes a ScenarioInterface object to file which consists of this and the
        supplied parameters. A ScenarioInterface object is capable of being run
        through the regression tester */
    void writeInputs(IModelSP      model,
                     CInstrumentSP instrument,
                     CControlSP    control,
                     MarketDataSP  market);
    /** Writes a ScenarioInterface object to file which consists of this and the
        supplied parameters. A ScenarioInterface object is capable of being run
        through the regression tester */
    void writeInputs(IModelSP                model,
                     IInstrumentCollectionSP instruments,
                     CControlSP              control,
                     MarketDataSP            market);
    /** Validation */
    virtual void validatePop2Object();      
private:
    Scenario();
    Scenario(const Scenario &rhs);
    Scenario& operator=(const Scenario& rhs);
    static void load(CClassSP& clazz);
    static IObject* defaultScenario();

    IScenarioShiftArraySP shifts;

};

// typedef for smart pointers to Scenario
typedef smartConstPtr<Scenario> ScenarioConstSP;
typedef smartPtr<Scenario> ScenarioSP;

DRLIB_END_NAMESPACE
#endif
