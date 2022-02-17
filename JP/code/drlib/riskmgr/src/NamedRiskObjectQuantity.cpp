//----------------------------------------------------------------------------
//
//   Group       : QR - Credit Hybrids
//
//   Filename    : NamedRiskObjectQuantity.cpp
//
//   Description : Like NamedRiskQuantity, but for returning IObjects.
//
//   Author      : Linus Thand
//
//   Date        : 26 July 2006
//
//
//----------------------------------------------------------------------------


#include "edginc/config.hpp"
#include "edginc/TRACE.hpp"
#include "edginc/Void.hpp"
#include "edginc/Untweakable.hpp"
#include "edginc/Results.hpp"
#include "edginc/IResultsIdentifier.hpp"
#include "edginc/NamedRiskObjectQuantity.hpp"
#include "edginc/NotApplicableException.hpp"
#include "edginc/DefaultConstructor.hpp"
#include "edginc/RiskQuantity.hpp"
#include "edginc/HypotheticalQuantity.hpp"

DRLIB_BEGIN_NAMESPACE

//// IObjectConstantRiskQuantity

class IObjectConstantRiskQuantity: public RiskQuantity {

    static void load(CClassSP& clazz) {
        REGISTER(IObjectConstantRiskQuantity, clazz);
        SUPERCLASS(RiskQuantity);
        EMPTY_SHELL_METHOD(DefaultConstructor<IObjectConstantRiskQuantity>::iObject);
    }

public:
    static CClassConstSP const TYPE;

    IObjectConstantRiskQuantity():
        RiskQuantity(TYPE, HypotheticalQuantityArray::SP())
    {}

    double value(const CDoubleArray&, const CDoubleArray&) const {
        throw ModelException("IObjectConstantRiskQuantity::value shouldn't be used this way");
    }

    virtual bool isConstant() const { return true; }
};

CClassConstSP const IObjectConstantRiskQuantity::TYPE = CClass::registerClassLoadMethod(
    "IObjectConstantRiskQuantity", typeid(IObjectConstantRiskQuantity), load);






////NamedRiskObjectQuantity

NamedRiskObjectQuantity::NamedRiskObjectQuantity(
        IObjectSP value,
        IResultsIdentifierConstSP resultsName):
        value(value),
    NamedRiskQuantity(RiskQuantityConstSP(new IObjectConstantRiskQuantity()), resultsName, 1)
{}

NamedRiskObjectQuantitySP NamedRiskObjectQuantity::SP(
        IObjectSP value,
        IResultsIdentifierConstSP resultsName) {
    return NamedRiskObjectQuantitySP(new NamedRiskObjectQuantity(
        value, resultsName));
}

NamedRiskObjectQuantity::~NamedRiskObjectQuantity() {}

void NamedRiskObjectQuantity::storeResult(const CDoubleArray& vals,
                                    const CDoubleArray& dsts,
                                    CResultsSP results) const {
    TRACE_METHOD;

    ASSERT(vals.size() == dsts.size());

    if (resultsName->exists(results)) {
        TRACE(*resultsName << " is already in the Results: just leave it");
    }
    else {
        TRACE("Calculating " << *resultsName << " from previously "
              "computed values");

        for (int w = 0; w < vals.size(); ++w) {
            TRACE("Price [or whatever] in world " << w << " = " << vals[w]);
            TRACE("Distance of world " << w << " from base case = " << dsts[w]);
        }

        try {
            resultsName->store(results,  value);
        }
        catch (NotApplicableException&) {
            // ... thrown by RiskQuantity::notApplicable()

            TRACE("So the result is NotApplicable");
            resultsName->storeNotApplicableToName(results);
        }
        catch (ModelException& e) {
            // ... thrown by RiskQuantity::untweakable() and also e.g.
            // IScalarDerivative::oneSided() on zero divisor

            TRACE("So the result is Untweakable");
            resultsName->storeUntweakable(
                results, UntweakableConstSP(new Untweakable(e)));
        }
        catch (exception& e) {
            TRACE("Oops, something went wrong => Untweakable");
            resultsName->storeUntweakable(
                results, UntweakableConstSP(new Untweakable(string(e.what()))));
        }
    }
                                      
    
}

static IObject* defaultNamedRiskObjectQuantity() {
    return new NamedRiskObjectQuantity(IObjectSP(),
                                 IResultsIdentifierConstSP());
}

void NamedRiskObjectQuantity::load(CClassSP& clazz) {
    REGISTER(NamedRiskQuantity, clazz);
    SUPERCLASS(NamedRiskQuantity);
    EMPTY_SHELL_METHOD(defaultNamedRiskObjectQuantity);
}

CClassConstSP const NamedRiskObjectQuantity::TYPE = CClass::registerClassLoadMethod(
    "NamedRiskObjectQuantity", typeid(NamedRiskObjectQuantity), load);

DEFINE_TEMPLATE_TYPE(NamedRiskObjectQuantityArray);

DRLIB_END_NAMESPACE
