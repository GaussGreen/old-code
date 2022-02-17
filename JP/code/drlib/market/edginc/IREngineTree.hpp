//----------------------------------------------------------------------------
//
//   Group       : QR Rates
//
//   Filename    : IREngineTree.hpp
//
//   Description : Tree Engine
//
//   Author      : Anwar E Sidat
//
//   Date        : 18-Aug-2006
//
//----------------------------------------------------------------------------

#ifndef QLIB_IREngineTree_HPP
#define QLIB_IREngineTree_HPP

#include "edginc/AtomicArray.hpp"
#include "edginc/IRModelConfig.hpp"


DRLIB_BEGIN_NAMESPACE

/** VNFM Model Class.
 */
class MARKET_DLL IREngineTree : public IRModelConfig
{
public:
    static CClassConstSP const TYPE;

    IREngineTree();
    virtual ~IREngineTree();

    /** Returns name of model. */
    virtual string getName() const;

    /** overrides default */
    virtual void validatePop2Object();

    //Fields
    string       Key;               // Handle name
    int          nbPPY;             // Number of Points Per Year
    int          nbCET;             // Number of Calibration Enhancement Tool iterations
    int          nbStdDevs;         // Number of Standard Deviations
    int          nbStateVars;       // Number of state variables
    int          nbStateVarStdDevs; // Number of state variable std deviation
    double       nbStdDevsMQ;       // Number of std deviations MultiQ
    double       MQNCK;             // MQ's numerical configuration key

protected:

    IREngineTree(const CClassConstSP& clazz);
    IREngineTree(const IREngineTree& irv);
    IREngineTree& operator=(const IREngineTree& irv);

private:
    static void load(CClassSP& clazz);
    static IObject* defaultConstructor(void) { return new IREngineTree(); }
};

typedef smartConstPtr<IREngineTree> IREngineTreeConstSP;
typedef smartPtr<IREngineTree>      IREngineTreeSP;
typedef MarketWrapper<IREngineTree> IREngineTreeWrapper;

#ifndef QLIB_IREngineTree_CPP
EXTERN_TEMPLATE(class MARKET_DLL_SP smartConstPtr<IREngineTree>);
EXTERN_TEMPLATE(class MARKET_DLL_SP smartPtr<IREngineTree>);
#else
INSTANTIATE_TEMPLATE(class MARKET_DLL smartConstPtr<IREngineTree>);
INSTANTIATE_TEMPLATE(class MARKET_DLL smartPtr<IREngineTree>);
#endif

// support for wrapper class
#ifndef QLIB_IREngineTree_CPP
EXTERN_TEMPLATE(class MARKET_DLL MarketWrapper<IREngineTree>);
#else
INSTANTIATE_TEMPLATE(class MARKET_DLL MarketWrapper<IREngineTree>);
#endif

DRLIB_END_NAMESPACE

#endif
