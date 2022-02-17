//----------------------------------------------------------------------------
//
//   Group       : QR Rates
//
//   Filename    : IRSmile2Q.cpp
//
//   Description : Smile2Q.
//
//   Author      : Anwar E Sidat
//
//   Date        : 21-Aug-2006
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#define QLIB_IRSmile2Q_CPP
#include "edginc/IRSmile2Q.hpp"
#include "edginc/Addin.hpp"
#include "edginc/Maths.hpp"

DRLIB_BEGIN_NAMESPACE

// Make sure the class links
bool IRSmile2QLoad() { return (IRSmile2Q::TYPE != 0); }

IRSmile2Q::IRSmile2Q()
    : IRExoticParam(TYPE)
{
}

IRSmile2Q::~IRSmile2Q() {}

IRSmile2Q::IRSmile2Q(const CClassConstSP& clazz): 
    IRExoticParam(clazz) {}

void IRSmile2Q::validatePop2Object()
{
    static const string method = "IRSmile2Q::validatePop2Object";
    try
    {
    }
    catch (exception& e)
    {
        throw ModelException(e, method);
    }
}


/************ Smile2Q Tweak **************/

string IRSmile2Q::sensName(const IRSmile2QTweak::Property*) const {
    return Key;
}

TweakOutcome IRSmile2Q::sensShift (const PropertyTweak<IRSmile2QTweak::Property>& shift) {
    try {
        if (!Maths::isZero(shift.coefficient)) {
            QLeft += shift.coefficient;
            QRight += shift.coefficient;
            FwdShift += shift.coefficient;
        }
        return TweakOutcome(shift.coefficient, false);
    } 
    catch (exception& e) {
        throw ModelException(e, __FUNCTION__,
            "Failed for " + getName());
    }
}


void IRSmile2Q::load(CClassSP& clazz)
{
    clazz->setPublic();
    clazz->setDescription("Smile2Q");
    REGISTER(IRSmile2Q, clazz);
    SUPERCLASS(IRExoticParam);
    EMPTY_SHELL_METHOD(defaultConstructor);
    FIELD(Key,      "Identifier for this object in market data cache.");
    FIELD(QLeft,    "Q-Left parameter (0..1).");
    FIELD(QRight,   "Q-Right parameter (0..1).");
    FIELD(FwdShift, "Forward Shift.");

    Addin::registerConstructor("IRSmile2Q",
                                Addin::MARKET,
                                "Creates a handle for a Smile2Q ie.Lognormal (QLeft=QRight=FwdShift=0), Normal (QLeft=QRight=1, FwdShift=0).",
                                IRSmile2Q::TYPE);
}

CClassConstSP const IRSmile2Q::TYPE = CClass::registerClassLoadMethod(
    "IRSmile2Q", typeid(IRSmile2Q), IRSmile2Q::load);

/** Returns name of struct */
string IRSmile2Q::getName() const
{
    return Key;
}
// definition of TYPE for MarketWrapper template class
DEFINE_TEMPLATE_TYPE(IRSmile2QWrapper);

DRLIB_END_NAMESPACE
