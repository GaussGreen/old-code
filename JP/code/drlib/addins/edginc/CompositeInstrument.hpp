//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : CompositeInstrument.hpp
//
//   Description : Defines a composite instrument
//
//   Author      : Andrew J Swain
//
//   Date        : 1 May 2001
//
//
//----------------------------------------------------------------------------


#ifndef COMPOSITEINSTRUMENT_HPP
#define COMPOSITEINSTRUMENT_HPP
#include "edginc/AtomicArray.hpp"
#include "edginc/Scenario.hpp"
#include "edginc/Control.hpp"
#include "edginc/Instrument.hpp"
#include "edginc/Model.hpp"
#include "edginc/Results.hpp"
#include "edginc/RegressionTest.hpp"
#include "edginc/ResultsSet.hpp"
#include "edginc/ClientRunnable.hpp"

DRLIB_BEGIN_NAMESPACE
/**  Defines a composite instrument */
class ADDINS_DLL CompositeInstrument: public CObject,
                           virtual public IRegressionTest,
                           virtual public ClientRunnable {
public:
    static CClassConstSP const TYPE;
   
    /** Does not clone supplied parameters */
    CompositeInstrument(ScenarioSP    scenario,
                        ObjectArraySP model, 
                        ObjectArraySP inst, 
                        ObjectArraySP ctrl,
                        DoubleArraySP multipliers,
                        DoubleArraySP weight,
                        CMarketDataSP market);

    /** Calculates price and sensitivities for given instrument and
        context together with supplied market data + any scenario shifts.
        Results are returned in either a ResultSet object. ( >1 instrument) 
        or a Results object (1 instrument) */
    IObject* run() const;
    
    // EdrAction 
    virtual IObjectSP run();

    /** Runs 'regression test' */
    virtual IObjectSP runTest() const;

private:
    // create a de-wrappered version of input file
    void fileConvert() const;

    // fields
    ScenarioSP     scenario;

    // any number of items in the ObjectArraySPs may be DRWrappers
    ObjectArraySP model;
    ObjectArraySP inst;
    ObjectArraySP ctrl;

    DoubleArraySP  multiplier;
    DoubleArraySP  weight;
    CMarketDataSP  market;

    friend class CompositeInstrumentHelper;
    friend class XLGetWrappers;
    friend class XLGetYieldCurveNames;
    friend class XLGetMarket;
    friend class FileInstTypeAddin;
    CompositeInstrument();
    CompositeInstrument(const CompositeInstrument &rhs);
    CompositeInstrument& operator=(const CompositeInstrument& rhs);
};

// typedef for smart pointers to CompositeInstrument
typedef smartConstPtr<CompositeInstrument> CompositeInstrumentConstSP;
typedef smartPtr<CompositeInstrument> CompositeInstrumentSP;

DRLIB_END_NAMESPACE
#endif
