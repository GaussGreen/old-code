//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : SpreadCurve.hpp
//
//   Description : a liquidity spread curve
//
//   Author      : André Segger
//
//
//----------------------------------------------------------------------------

#ifndef EDR_LIQUIDITY_SPREAD_CURVE_H
#define EDR_LIQUIDITY_SPREAD_CURVE_H

#include "edginc/Class.hpp"
#include "edginc/MarketObject.hpp"
#include "edginc/Expiry.hpp"
#include "edginc/LiquiditySpreadRhoParallel.hpp"
#include "edginc/LiquiditySpreadPointwise.hpp"
#include "edginc/IRestorableWithRespectTo.hpp"

DRLIB_BEGIN_NAMESPACE

class MARKET_DLL LiquiditySpreadCurve : public MarketObject {
public:
// not done yet
    static CClassConstSP const TYPE;
    static void load(CClassSP& );

    LiquiditySpreadCurve();

    LiquiditySpreadCurve(const string&      name,
                         const ExpiryArray* expiries,
                         const DoubleArray& rates);

    LiquiditySpreadCurve(const string&      name,
                         const ExpiryArray* expiries,
                         const DoubleArray& rates,
                         const int          parSwapFreq,
                         const string&      parDCC,
                         const double       parRecovery,
                         const bool         parAccrueFee);

    LiquiditySpreadCurve(const CClassConstSP& clazz);

    virtual string getName() const
    {
        return name;
    }

    /** SpreadRhoParallel support */
    string sensName(LiquiditySpreadRhoParallel* shift) const;
    bool sensShift(LiquiditySpreadRhoParallel* shift);
    void sensRestore(LiquiditySpreadRhoParallel* shift);

    /** LiquiditySpreadRhoPointwise support */
    string sensName(const LiquiditySpreadPointwise*) const;
    ExpiryWindowArrayConstSP sensQualifiers(const LiquiditySpreadPointwise*) const;
    TweakOutcome sensShift(const PropertyTweak<LiquiditySpreadPointwise>&);
    void sensRestore(const PropertyTweak<LiquiditySpreadPointwise>&);
    
    ExpiryArrayConstSP getExpiries()const;

    DateTimeArray getExpiryDates(const DateTime& today)const;

    void validate();

protected:
    string          name;

    ExpiryArraySP   expiries;
    CDoubleArray    liquiditySpreads;

    int             parSwapFreq;
    string          parDCC;
    double          parRecovery;
    bool            parAccrueFee;
};

typedef smartConstPtr<LiquiditySpreadCurve> LiquiditySpreadCurveConstSP;
typedef smartPtr<LiquiditySpreadCurve>      LiquiditySpreadCurveSP;
#ifndef QLIB_LIQUIDITYSPREADCURVE_CPP
EXTERN_TEMPLATE(class MARKET_DLL_SP smartConstPtr<LiquiditySpreadCurve>);
EXTERN_TEMPLATE(class MARKET_DLL_SP smartPtr<LiquiditySpreadCurve>);
#else
INSTANTIATE_TEMPLATE(class MARKET_DLL smartConstPtr<LiquiditySpreadCurve>);
INSTANTIATE_TEMPLATE(class MARKET_DLL smartPtr<LiquiditySpreadCurve>);
#endif
// support for wrapper class
typedef MarketWrapper<LiquiditySpreadCurve> LiquiditySpreadCurveWrapper;
#ifndef QLIB_LIQUIDITYSPREADCURVE_CPP
EXTERN_TEMPLATE(class MARKET_DLL MarketWrapper<LiquiditySpreadCurve>);
#else
INSTANTIATE_TEMPLATE(class MARKET_DLL MarketWrapper<LiquiditySpreadCurve>);
#endif
    
DRLIB_END_NAMESPACE
#endif
