//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : RiskMgrInterface.hpp
//
//   Description : Defines interface to RiskMgr
//
//   Author      : Andrew J Swain
//
//   Date        : 19 January 2001
//
//
//----------------------------------------------------------------------------

#ifndef RISKMGRINTERFACE_HPP
#define RISKMGRINTERFACE_HPP

#include "edginc/RegressionTest.hpp"
#include "edginc/ClientRunnable.hpp"
#include "edginc/Model.hpp"
#include "edginc/MarketData.hpp"
#include "edginc/Instrument_forward.hpp"

DRLIB_BEGIN_NAMESPACE

class RISKMGR_DLL RiskMgrInterface : public CObject,
                         virtual public IRegressionTest,
                         virtual public ClientRunnable {
public:
    static CClassConstSP const TYPE;
    friend class RiskMgrInterfaceHelper;
    
    /** Does not clone supplied parameters */
    RiskMgrInterface(IModelSP      model, CInstrumentSP inst, CControlSP ctrl,
                     CMarketDataSP market);

    IModelSP      model;
    CInstrumentSP inst;
    CControlSP    ctrl;
    CMarketDataSP market;

    /** Runs 'regression test' for model, inst and control in this class */
    virtual  IObjectSP runTest() const;

    // EdrAction 
    virtual IObjectSP run();


private:
    /** for addin - runs inputs through risk manager */
    static IObjectSP addinRiskMgr(RiskMgrInterface* addinParams);
    IObjectSP go() const;

    RiskMgrInterface();
    RiskMgrInterface(const RiskMgrInterface& rhs);
    RiskMgrInterface& operator=(const RiskMgrInterface& rhs);
};

// typedef for smart pointers to RiskMgrInterface
typedef smartConstPtr<RiskMgrInterface> RiskMgrInterfaceConstSP;
typedef smartPtr<RiskMgrInterface> RiskMgrInterfaceSP;

DRLIB_END_NAMESPACE

#endif
