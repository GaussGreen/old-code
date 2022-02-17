/**
 * @file RiskQuantity.cpp
 */

#include "edginc/config.hpp"
#include "edginc/TRACE.hpp"
#include "edginc/DefaultConstructor.hpp"
#include "edginc/RiskQuantity.hpp"
#include "edginc/HypotheticalQuantity.hpp"
#include "edginc/NotApplicableException.hpp"

DRLIB_BEGIN_NAMESPACE

RiskQuantity::RiskQuantity(CClassConstSP type,
                           HypotheticalQuantityArrayConstSP parameters):
    CObject(type),
    _parameters(parameters)
{}

RiskQuantity::RiskQuantity(CClassConstSP type):
    CObject(type)
{}

RiskQuantity::~RiskQuantity() {}

HypotheticalQuantityArrayConstSP RiskQuantity::parameters() const {
    return _parameters;
}

bool RiskQuantity::isExceptional() const {
    return false;
}

bool RiskQuantity::isNotApplicable() const {
    return false;
}

class ConstantRiskQuantity: public RiskQuantity {

    static void load(CClassSP& clazz) {
        REGISTER(ConstantRiskQuantity, clazz);
        SUPERCLASS(RiskQuantity);
        FIELD(c, "c");
        EMPTY_SHELL_METHOD(DefaultConstructor<ConstantRiskQuantity>::iObject);
    }

public:
    static CClassConstSP const TYPE;

    double c;

    ConstantRiskQuantity(double c = 0):
        RiskQuantity(TYPE, HypotheticalQuantityArray::SP()),
        c(c)
    {}

    double value(const CDoubleArray&, const CDoubleArray&) const {
        TRACE_METHOD;
        TRACE("Value of quantity is just constant " << c);
        return c;
    }

    virtual bool isConstant() const { return true; }
};

DECLARE(ConstantRiskQuantity)

CClassConstSP const ConstantRiskQuantity::TYPE = CClass::registerClassLoadMethod(
    "ConstantRiskQuantity", typeid(ConstantRiskQuantity), load);

RiskQuantityConstSP RiskQuantity::constant(double c) {
    return RiskQuantityConstSP(new ConstantRiskQuantity(c));
}

class NotApplicableRiskQuantity: public RiskQuantity {

    static void load(CClassSP& clazz) {
        REGISTER(NotApplicableRiskQuantity, clazz);
        SUPERCLASS(RiskQuantity);
        EMPTY_SHELL_METHOD(DefaultConstructor<NotApplicableRiskQuantity>::iObject);
    }

public:
    static CClassConstSP const TYPE;

    NotApplicableRiskQuantity():
        RiskQuantity(TYPE, HypotheticalQuantityArray::SP())
    {}

    double value(const CDoubleArray&, const CDoubleArray&) const {
        TRACE_METHOD;
        TRACE("Value of quantity is 'not applicable'");
        throw NotApplicableException();
    }

    bool isExceptional() const { return true; }

    bool isNotApplicable() const { return true; }
};

DECLARE(NotApplicableRiskQuantity)

CClassConstSP const NotApplicableRiskQuantity::TYPE = CClass::registerClassLoadMethod(
    "NotApplicableRiskQuantity", typeid(NotApplicableRiskQuantity), load);

RiskQuantityConstSP RiskQuantity::notApplicable() {
    return RiskQuantityConstSP(new NotApplicableRiskQuantity());
}

class UntweakableRiskQuantity: public RiskQuantity {

    static void load(CClassSP& clazz) {
        REGISTER(UntweakableRiskQuantity, clazz);
        SUPERCLASS(RiskQuantity);
        FIELD(message, "message");
        EMPTY_SHELL_METHOD(DefaultConstructor<UntweakableRiskQuantity>::iObject);
    }

public:
    static CClassConstSP const TYPE;

    string message;

    UntweakableRiskQuantity(): RiskQuantity(TYPE) {}

    UntweakableRiskQuantity(const ModelException& e):
        RiskQuantity(TYPE, HypotheticalQuantityArray::SP())
    {
        char* m = e.stackTrace();
        message = m;
        free(m);
    }

    double value(const CDoubleArray&, const CDoubleArray&) const {
        TRACE_METHOD;
        TRACE("Value of quantity is 'untweakable'");
        throw ModelException(message);
    }

	virtual bool isExceptional() const { return true; }
};

DECLARE(UntweakableRiskQuantity)

CClassConstSP const UntweakableRiskQuantity::TYPE = CClass::registerClassLoadMethod(
    "UntweakableRiskQuantity", typeid(UntweakableRiskQuantity), load);

RiskQuantityConstSP RiskQuantity::untweakable(const ModelException& e) {
    return RiskQuantityConstSP(new UntweakableRiskQuantity(e));
}

void RiskQuantity::load(CClassSP& clazz) {
    REGISTER(RiskQuantity, clazz);
    SUPERCLASS(CObject);
    FIELD(_parameters, "parameters");
}

CClassConstSP const RiskQuantity::TYPE =
    CClass::registerClassLoadMethod(
        "RiskQuantity", typeid(RiskQuantity), load);

DRLIB_END_NAMESPACE
