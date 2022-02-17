//----------------------------------------------------------------------------
//
//   Group       : QR Rates
//
//   Filename    : IREngineMC.hpp
//
//   Description : Monte Carlo (European) Engine
//
//   Author      : Anwar E Sidat
//
//   Date        : 18-Aug-2006
//
//----------------------------------------------------------------------------

#ifndef QLIB_IREngineMC_HPP
#define QLIB_IREngineMC_HPP

//#include "edginc/AtomicArray.hpp"
#include "edginc/IRModelConfig.hpp"


DRLIB_BEGIN_NAMESPACE

/** IREngineMC Engine Class.
 */
class MARKET_DLL IREngineMC : public IRModelConfig
{
public:
    static CClassConstSP const TYPE;

    IREngineMC();
    virtual ~IREngineMC();

    /** Returns name of model. */
    virtual string getName() const;

    /** overrides default */
    virtual void validatePop2Object();

    //Fields
    string       Key;       // Handle name
    int          nbPPY;     // Number of Points per Year
    int          nbCET;     // Number of Calibration Enhancement Tool iterations
    double       nbPaths;   // Number of paths
    int          Dice;      // Dice value
    double       Seed;      // Seed value
    int          nbStdDevs; // Number of Standard Deviations
    int          nbRuns;    // Number of runs (not paths)

protected:

    IREngineMC(const CClassConstSP& clazz);
    IREngineMC(const IREngineMC& irv);
    IREngineMC& operator=(const IREngineMC& irv);

private:
    static void load(CClassSP& clazz);
    static IObject* defaultConstructor(void) { return new IREngineMC(); }
};

typedef smartConstPtr<IREngineMC> IREngineMCConstSP;
typedef smartPtr<IREngineMC>      IREngineMCSP;
typedef MarketWrapper<IREngineMC> IREngineMCWrapper;

#ifndef QLIB_IREngineMC_CPP
EXTERN_TEMPLATE(class MARKET_DLL_SP smartConstPtr<IREngineMC>);
EXTERN_TEMPLATE(class MARKET_DLL_SP smartPtr<IREngineMC>);
#else
INSTANTIATE_TEMPLATE(class MARKET_DLL smartConstPtr<IREngineMC>);
INSTANTIATE_TEMPLATE(class MARKET_DLL smartPtr<IREngineMC>);
#endif

// support for wrapper class
#ifndef QLIB_IREngineMC_CPP
EXTERN_TEMPLATE(class MARKET_DLL MarketWrapper<IREngineMC>);
#else
INSTANTIATE_TEMPLATE(class MARKET_DLL MarketWrapper<IREngineMC>);
#endif

DRLIB_END_NAMESPACE

#endif
