//----------------------------------------------------------------------------
//
//   Group       : QR Rates
//
//   Filename    : IRSmile2Q.hpp
//
//   Description : Smile 2Q
//
//   Author      : Anwar E Sidat
//
//   Date        : 18-Aug-2006
//
//----------------------------------------------------------------------------

#ifndef QLIB_IRSmile2Q_HPP
#define QLIB_IRSmile2Q_HPP

#include "edginc/IRExoticParam.hpp"
#include "edginc/IRSmile2QTweak.hpp"
#include "edginc/TweakOutcome.hpp"
#include "edginc/PropertyTweak.hpp"
#include "edginc/ITweakableWithRespectTo.hpp"

DRLIB_BEGIN_NAMESPACE

/** Smile2Q Class.
 */
class MARKET_DLL IRSmile2Q : public IRExoticParam,
    virtual public ITweakableWithRespectTo<IRSmile2QTweak::Property>
{
public:
    static CClassConstSP const TYPE;

    IRSmile2Q();
    virtual ~IRSmile2Q();

    /** Get methods for input parameters. */
    virtual double qLeft() const    { return QLeft; }
    virtual double qRight() const   { return QRight; }
    virtual double fwdShift() const { return FwdShift; }

    /** Returns name of model. */
    virtual string getName() const;

    /** overrides default */
    virtual void validatePop2Object();

    /************ Smile2Q Tweak **************/
    string sensName(const IRSmile2QTweak::Property*) const;
    TweakOutcome sensShift (const PropertyTweak<IRSmile2QTweak::Property>& shift);


protected:

    IRSmile2Q(const CClassConstSP& clazz);
    IRSmile2Q(const IRSmile2Q& irv);
    IRSmile2Q& operator=(const IRSmile2Q& irv);

    //Fields
    string       Key;        // Handle name
    double       QLeft;
    double       QRight;
    double       FwdShift;

private:
    static void load(CClassSP& clazz);
    static IObject* defaultConstructor(void) { return new IRSmile2Q(); }
};

typedef smartConstPtr<IRSmile2Q> IRSmile2QConstSP;
typedef smartPtr<IRSmile2Q>      IRSmile2QSP;
typedef MarketWrapper<IRSmile2Q> IRSmile2QWrapper;

#ifndef QLIB_IRSmile2Q_CPP
EXTERN_TEMPLATE(class MARKET_DLL_SP smartConstPtr<IRSmile2Q>);
EXTERN_TEMPLATE(class MARKET_DLL_SP smartPtr<IRSmile2Q>);
#else
INSTANTIATE_TEMPLATE(class MARKET_DLL smartConstPtr<IRSmile2Q>);
INSTANTIATE_TEMPLATE(class MARKET_DLL smartPtr<IRSmile2Q>);
#endif

// support for wrapper class
#ifndef QLIB_IRSmile2Q_CPP
EXTERN_TEMPLATE(class MARKET_DLL MarketWrapper<IRSmile2Q>);
#else
INSTANTIATE_TEMPLATE(class MARKET_DLL MarketWrapper<IRSmile2Q>);
#endif

DRLIB_END_NAMESPACE

#endif
