/**
 * @file MultiRiskMgrInterface.hpp
 */

#ifndef DRLIB_MultiRiskMgrInterface_H
#define DRLIB_MultiRiskMgrInterface_H

#include "edginc/RegressionTest.hpp"
#include "edginc/ClientRunnable.hpp"
#include "edginc/Model.hpp"
#include "edginc/MarketData.hpp"

using namespace std;

DRLIB_BEGIN_NAMESPACE

/**
 * Client entry point for pricing/greeking multiple instruments in one run
 *
 * This is the external interface to QLib's "vectorized pricing" facility (see
 * IInstrumentCollection for an overview).
 *
 * MultiRiskMgrInterface closely parallels RiskMgrInterface, but accepts an
 * "array" of instruments rather than just one, and returns an array of
 * Results, one for each.
 *
 * The input IInstrumentCollection can be an actual array
 * (ArrayInstrumentCollection) or something more specialised (such as
 * VanillaGridInstrumentCollection).  Using a specialised collection
 * may be more convenient for clients, and allows QLib to leverage
 * knowledge of the structure of the collection to obtain performance
 * benefits.
 *
 * Note that the single-instrument RiskMgrInterface is now a wrapper
 * round MultiRiskMgrInterface.
 */

class RISKMGR_DLL MultiRiskMgrInterface : public CObject,
                              virtual public IRegressionTest,
                              virtual public ClientRunnable {
    static void load(CClassSP& clazz);

public:
    static CClassConstSP const TYPE;
    
    /** Does not clone supplied parameters */
    MultiRiskMgrInterface(IModelSP model,
                          IInstrumentCollectionSP insts,
                          CControlSP ctrl,
                          CMarketDataSP market);

    /** Runs 'regression test' for model, insts and control in this class */
    virtual IObjectSP runTest() const;

    // EdrAction 
    virtual IObjectSP run();

private:
    MultiRiskMgrInterface();
    static IObject *defaultMultiRiskMgrInterface();
    MultiRiskMgrInterface(const MultiRiskMgrInterface& rhs);
    MultiRiskMgrInterface& operator=(const MultiRiskMgrInterface& rhs);
    IObjectSP go() const;

    //// fields
    IModelSP                model;
    IInstrumentCollectionSP insts;
    CControlSP              ctrl;
    CMarketDataSP           market;
};

// typedef for smart pointers to MultiRiskMgrInterface
typedef smartConstPtr<MultiRiskMgrInterface> MultiRiskMgrInterfaceConstSP;
typedef smartPtr<MultiRiskMgrInterface> MultiRiskMgrInterfaceSP;

DRLIB_END_NAMESPACE

#endif
