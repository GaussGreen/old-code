//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : DRWrapperInterface.hpp
//
//   Description : Defines DR Wrapper interface to RiskMgr
//
//   Author      : Andrew J Swain
//
//   Date        : 20 April 2001
//
//
//----------------------------------------------------------------------------

#ifndef DRWRAPPERINTERFACE_HPP
#define DRWRAPPERINTERFACE_HPP

#include "edginc/RegressionTest.hpp"
#include "edginc/Scenario.hpp"

DRLIB_BEGIN_NAMESPACE

class ADDINS_DLL DRWrapperInterface : public CObject,
                           public IRegressionTest{
public:
    static CClassConstSP const TYPE;
    
    /** Does not clone supplied parameters */
    DRWrapperInterface(ScenarioSP    scenario,
                       IObjectSP     model, 
                       IObjectSP     inst,
                       IObjectSP     ctrl,
                       CMarketDataSP market);

    ScenarioSP    scenario;
    IObjectSP     model;
    IObjectSP     inst;
    IObjectSP     ctrl;
    CMarketDataSP market;

    /** Runs 'regression test' for model, inst and control in this class */
    virtual  IObjectSP runTest() const;

private:
    friend class DRWrapperInterfaceHelper;

    /** for addin - runs inputs through risk manager */
    static IObjectSP addinDRWrapper(DRWrapperInterface* addinParams);

    DRWrapperInterface();
    DRWrapperInterface(const DRWrapperInterface& rhs);
    DRWrapperInterface& operator=(const DRWrapperInterface& rhs);
};

// typedef for smart pointers to DRWrapperInterface
typedef smartConstPtr<DRWrapperInterface> DRWrapperInterfaceConstSP;
typedef smartPtr<DRWrapperInterface> DRWrapperInterfaceSP;

DRLIB_END_NAMESPACE

#endif
