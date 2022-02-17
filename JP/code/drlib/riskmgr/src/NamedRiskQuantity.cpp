/**
 * @file NamedRiskQuantity.cpp
 */

#include "edginc/config.hpp"
#include "edginc/TRACE.hpp"
#include "edginc/Void.hpp"
#include "edginc/Untweakable.hpp"
#include "edginc/Results.hpp"
#include "edginc/IResultsIdentifier.hpp"
#include "edginc/NamedRiskQuantity.hpp"
#include "edginc/NotApplicableException.hpp"

DRLIB_BEGIN_NAMESPACE

NamedRiskQuantity::NamedRiskQuantity(
        RiskQuantityConstSP riskQuantity,
        IResultsIdentifierConstSP resultsName,
        double unit):
    CObject(TYPE),
    riskQuantity(riskQuantity),
    resultsName(resultsName),
    unit(unit)
{}

NamedRiskQuantitySP NamedRiskQuantity::SP(
        RiskQuantityConstSP riskQuantity,
        IResultsIdentifierConstSP resultsName,
        double unit) {
    return NamedRiskQuantitySP(new NamedRiskQuantity(
        riskQuantity, resultsName, unit));
}

NamedRiskQuantity::~NamedRiskQuantity() {}

void NamedRiskQuantity::storeResult(const CDoubleArray& vals,
                                    const CDoubleArray& dsts,
                                    CResultsSP results) const {
    TRACE_METHOD;

    ASSERT(vals.size() == dsts.size());

    if (resultsName->exists(results)) {
        TRACE(*resultsName << " is already in the Results: just leave it");
    }
    else if (riskQuantity->isNotApplicable()) {
        // Typically, RiskQuantityFactorySensitivity added a
        // RiskQuantity::notApplicable() to signal that no sensitive market
        // data were found.  We don't strictly have to handle this case here,
        // since the "catch NotApplicableException" below would handle it; but
        // throwing and catching the exception takes rather a long time.

        TRACE(*resultsName << " was marked as NotApplicable");
        resultsName->storeNotApplicableToName(results);
    }
    else {
        TRACE("Calculating " << *resultsName << " from previously "
              "computed values");

        for (int w = 0; w < vals.size(); ++w) {
            TRACE("Price [or whatever] in world " << w << " = " << vals[w]);
            TRACE("Distance of world " << w << " from base case = " << dsts[w]);
        }

        try {
            double value = riskQuantity->value(vals, dsts);
            TRACE("So the result is " << value << " * " << unit << " = " <<
                  unit * value);
            resultsName->store(results, unit * value);
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

static IObject* defaultNamedRiskQuantity() {
    return new NamedRiskQuantity(RiskQuantityConstSP(),
                                 IResultsIdentifierConstSP());
}

void NamedRiskQuantity::load(CClassSP& clazz) {
    REGISTER(NamedRiskQuantity, clazz);
    SUPERCLASS(CObject);
    EMPTY_SHELL_METHOD(defaultNamedRiskQuantity);
    FIELD(riskQuantity, "riskQuantity");
    FIELD(resultsName, "resultsName");
    FIELD(unit, "unit");
}

CClassConstSP const NamedRiskQuantity::TYPE = CClass::registerClassLoadMethod(
    "NamedRiskQuantity", typeid(NamedRiskQuantity), load);

DEFINE_TEMPLATE_TYPE(NamedRiskQuantityArray);

DRLIB_END_NAMESPACE
