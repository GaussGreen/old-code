/**
 * @file RiskAxis.cpp
 */

#include "edginc/config.hpp"
#include "edginc/Void.hpp"
#include "edginc/Expiry.hpp"
#include "edginc/ExpiryWindow.hpp"
#include "edginc/OutputName.hpp"
#include "edginc/IRiskAxis.hpp"
#include "edginc/RiskAxis.hpp"
#include <map>

DRLIB_BEGIN_NAMESPACE

using std::map;

class SimpleFreezer: public CObject {
    static CClassConstSP const TYPE;

    static void load(CClassSP& clazz) {
        REGISTER(SimpleFreezer, clazz);
        SUPERCLASS(CObject);
    }

public:

    string className, marketNameFieldName, qualifierFieldName; // $unregistered

    SimpleFreezer(string className, string marketNameFieldName,
                  string qualifierFieldName):
        CObject(TYPE),
        className(className),
        marketNameFieldName(marketNameFieldName),
        qualifierFieldName(qualifierFieldName)
    {
        try {
            if (CClass::forName(className)->
                    getDeclaredField(qualifierFieldName)->getType() ==
                        Void::TYPE) {
                this->qualifierFieldName = "";
            }
        }
        catch (exception& e) {
            throw ModelException(e, "SimpleFreezer::SimpleFreezer()");
        }
    }

    RiskAxisConstSP frozen(IRiskAxisConstSP it) const;
};

CClassConstSP const SimpleFreezer::TYPE = CClass::registerClassLoadMethod(
    "SimpleFreezer", typeid(SimpleFreezer), load);

DECLARE(SimpleFreezer)

static map<string, SimpleFreezerSP> reg;

static IRiskAxisConstSP reg_thawed(string className,
                                   string marketName,
                                   ExpiryWindowSP expiry) {
    try {
        map<string, SimpleFreezerSP>::iterator i = reg.find(className);
        if (i == reg.end()) {
            throw ModelException("Class " + className + " " +
                                 "has not been registerSimple'd");
        }

        CClassConstSP clazz = CClass::forName(className);

        IObject* obj = clazz->newInstance();
        IRiskAxisSP axis(dynamic_cast<IRiskAxis*>(obj));
        if (!axis) {
            delete obj;
            throw ModelException("Class " + className + " is not an IRiskAxis "
                                 "implementation");
        }

        SimpleFreezer& freezer = *i->second;

        clazz->getDeclaredField(freezer.marketNameFieldName)->
            set(axis, OutputName::SP(marketName));

        if (!expiry) {
            if (freezer.qualifierFieldName != "") {
                throw ModelException(
                    "No expiry given but class " + className +
                    " has a qualifier field " + freezer.qualifierFieldName +
                    " registered through registerSimple()");
            }
        }
        else {
            clazz->getDeclaredField(freezer.qualifierFieldName)->
                set(axis, expiry);
        }

        axis->validatePop2Object();

        return axis;
    }
    catch (exception& e) {
        throw ModelException(e, "reg_thawed");
    }
}

RiskAxis::RiskAxis()
{}

RiskAxis::~RiskAxis() {}

class ScalarRiskAxis: public CObject,
                      public virtual RiskAxis {
public:

    static CClassConstSP const TYPE;

private:

    ScalarRiskAxis(): CObject(TYPE) {}

    static IObject* defaultScalarRiskAxis() {
        return new ScalarRiskAxis();
    }

    static void load(CClassSP& clazz) {
        clazz->setPublic();
        REGISTER(ScalarRiskAxis, clazz);
        SUPERCLASS(CObject);
        IMPLEMENTS(RiskAxis);
        EMPTY_SHELL_METHOD(defaultScalarRiskAxis);
        FIELD(className, "className");
        FIELD(marketDataName, "marketDataName");
    }

    string className, marketDataName;

public:

    ScalarRiskAxis(string className, string marketDataName):
        CObject(TYPE),
        className(className),
        marketDataName(marketDataName)
    {}

    IRiskAxisConstSP thawed() const {
        try {
            return reg_thawed(className, marketDataName, ExpiryWindowSP());
        }
        catch (exception& e) {
            throw ModelException(e, "ScalarRiskAxis::thawed()");
        }
    }
};

CClassConstSP const ScalarRiskAxis::TYPE = CClass::registerClassLoadMethod(
    "ScalarRiskAxis", typeid(ScalarRiskAxis), load);

class ExpiryRiskAxis: public CObject,
                      public virtual RiskAxis {
public:

    static CClassConstSP const TYPE;

private:

    ExpiryRiskAxis(): CObject(TYPE) {}

    static IObject* defaultExpiryRiskAxis() {
        return new ExpiryRiskAxis();
    }

    static void load(CClassSP& clazz) {
        clazz->setPublic();
        REGISTER(ExpiryRiskAxis, clazz);
        SUPERCLASS(CObject);
        IMPLEMENTS(RiskAxis);
        EMPTY_SHELL_METHOD(defaultExpiryRiskAxis);
        FIELD(className, "className");
        FIELD(marketDataName, "marketDataName");
        FIELD(previous, "previous");
        FIELD_MAKE_OPTIONAL(previous);
        FIELD(expiry, "expiry");
        FIELD(next, "next");
        FIELD_MAKE_OPTIONAL(next);
    }

    string className;
    string marketDataName;
    ExpiryConstSP previous, expiry, next;

public:

    ExpiryRiskAxis(string className,
                   string marketDataName,
                   ExpiryWindowConstSP expiry):
        CObject(TYPE),
        className(className),
        marketDataName(marketDataName),
        previous(expiry->previous),
        expiry(expiry->expiry),
        next(expiry->next)
    {}

    IRiskAxisConstSP thawed() const {
        try {
            return reg_thawed(className, marketDataName,
                              ExpiryWindow::SP(previous, expiry, next));
        }
        catch (exception& e) {
            throw ModelException(e, "ExpiryRiskAxis::thawed()");
        }
    }
};

CClassConstSP const ExpiryRiskAxis::TYPE = CClass::registerClassLoadMethod(
    "ExpiryRiskAxis", typeid(ExpiryRiskAxis), load);

void RiskAxis::registerSimple(string className,
                              string marketNameFieldName,
                              string qualifierFieldName) {
    reg[className] = SimpleFreezerSP(new SimpleFreezer(
        className, marketNameFieldName, qualifierFieldName));
}

void RiskAxis::registerSimple(string className,
                              string marketNameFieldName) {
    reg[className] = SimpleFreezerSP(new SimpleFreezer(
        className, marketNameFieldName, ""));
}

RiskAxisConstSP SimpleFreezer::frozen(IRiskAxisConstSP it) const {
    try {
        string marketName = OutputNameConstSP::dynamicCast(
                it->getClass()->getDeclaredField(marketNameFieldName)->
                    get(IRiskAxisSP::constCast(it))
            )->toString();

        if (qualifierFieldName == "") {
            return RiskAxisConstSP(new ScalarRiskAxis(className, marketName));
        }
        else {
            CFieldConstSP q = it->getClass()->getDeclaredField(
                                                  qualifierFieldName);
            if (q->getType() == Void::TYPE) {
                return RiskAxisConstSP(new ScalarRiskAxis(
                    className, marketName));
            }
            else if (q->getType() == ExpiryWindow::TYPE) {
                return RiskAxisConstSP(new ExpiryRiskAxis(
                    className,
                    marketName,
                    ExpiryWindowSP::dynamicCast(q->get(
                        IObjectSP::constCast(it)))));
            }
            else {
                throw ModelException(
                    "SimpleFreezer can't encode a 'qualifier' of type " +
                    q->getType()->getName());
            }
        }
    }
    catch (exception& e) {
        throw ModelException(e, "SimpleFreezer::frozen()",
                             "Freezing a '" + it->getClass()->getName() + "'");
    }
}

RiskAxisConstSP RiskAxis::frozenSimple(const IRiskAxis* it) {
    try {
        map<string, SimpleFreezerSP>::iterator i =
            reg.find(it->getClass()->getName());
        if (i == reg.end()) {
            throw ModelException("Class " + it->getClass()->getName() + " " +
                                 "has not been registerSimple'd");
        }

        return i->second->frozen(IRiskAxisConstSP::attachToRef(it));
    }
    catch (exception& e) {
        throw ModelException(e, "RiskAxis::frozenSimple()");
    }
}

static void loadRiskAxis(CClassSP& clazz) {
    REGISTER_INTERFACE(RiskAxis, clazz);
    clazz->setPublic();
    EXTENDS(IObject);
}

CClassConstSP const RiskAxis::TYPE = CClass::registerInterfaceLoadMethod(
    "RiskAxis", typeid(RiskAxis), loadRiskAxis);

DEFINE_TEMPLATE_TYPE(RiskAxisArray);

DRLIB_END_NAMESPACE

// Local Variables: ***
// c-basic-offset:4 ***
// End: ***
