//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : MegaShuffleInterface.hpp
//
//   Description : Decorates RiskMgr to enable a shuffled mega test.
//
//   Author      : Jon Dee
//
//   Date        : 22 November 2006
//
//
//----------------------------------------------------------------------------

#ifndef EDR_MEGA_SHUFFLE_INTERFACE
#define EDR_MEGA_SHUFFLE_INTERFACE

#include "edginc/RegressionTest.hpp"
#include "edginc/ClientRunnable.hpp"
#include "edginc/RiskMgrInterface.hpp"
#include "edginc/Control.hpp"
#include "edginc/Instrument.hpp"
#include "edginc/Model.hpp"
#include "edginc/Results.hpp"

DRLIB_BEGIN_NAMESPACE

FORWARD_DECLARE(Control);


/** decorates RiskMgrInterface, allows shuffled mega sens tests */
class RISKMGR_DLL MegaShuffleInterface : public CObject,
    virtual public IRegressionTest,
    virtual public ClientRunnable
{
    RiskMgrInterfaceSP riskMgr;
    int seed;
public:
    static CClassConstSP const TYPE;

    MegaShuffleInterface(IObjectSP _riskMgr, int _seed);
    virtual ~MegaShuffleInterface(void) {}

    static void load(CClassSP& clazz);
    static IObject* defaultInterface();

    //IRegressionTest
    virtual IObjectSP runTest() const;

    //Client runnable
    virtual IObjectSP run();

    // Returns a seed based on today's date. 
    // The seed is constant during a given day.
    static int createShuffleSeed();

private:
    MegaShuffleInterface();
    IObjectSP go() const;
    ControlSP makeControl(SensitivitySP sensToUse) const;
    void removeUnpredictableOutput(CControlSP control) const;

};

DECLARE(MegaShuffleInterface);

DRLIB_END_NAMESPACE

#endif


