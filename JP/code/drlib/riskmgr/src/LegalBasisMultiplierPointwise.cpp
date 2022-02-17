/**
 * @file LegalBasisMultiplierPointwise.cpp
 */

#include "edginc/config.hpp"
#include "edginc/DefaultConstructor.hpp"
#include "edginc/MultiTweakGroup.hpp"
#include "edginc/IResultsFunction.hpp"
#include "edginc/RiskProperty.hpp"
#include "edginc/IScalarDerivative.hpp"
#include "edginc/PerNameRiskPropertySensitivity.hpp"
#include "edginc/LegalBasisRelativePointwise.hpp"
#include "edginc/GenericSensitivityFactory.hpp"
#include "edginc/ParSpreadMaxTweakSize.hpp"

DRLIB_BEGIN_NAMESPACE

static const char* NAME = "LEGAL_BASIS_MULTIPLIER_POINTWISE";

/**
 * A greek calculated by tweaking each market name's par spreads at all its
 * defined expiries.
 *
 * LegalBasisMultiplierPointwise has been rewritten to use the "declarative" sensitivities
 * framework: see IRiskQuantityFactory for an overview.
 */

class LegalBasisMultiplierPointwise: public PerNameRiskPropertySensitivity<ExpiryWindow> {

    static void load(CClassSP& clazz) {
        clazz->setPublic();
        REGISTER(LegalBasisMultiplierPointwise, clazz);
        SUPERCLASS(PerNameRiskPropertySensitivity<ExpiryWindow>);
        EMPTY_SHELL_METHOD(&DefaultConstructor<LegalBasisMultiplierPointwise>::iObject);
        SensitivityFactory::addSens(NAME,
                                    new GenericSensitivityFactory<LegalBasisMultiplierPointwise>(), 
                                    new LegalBasisMultiplierPointwise(DEFAULT_SHIFT),
                                    ITweakableWithRespectTo<LegalBasisRelativePointwise>::TYPE);
    }

public:

    static CClassConstSP const TYPE;

    static const double DEFAULT_SHIFT;

    LegalBasisMultiplierPointwise(double shiftSize = DEFAULT_SHIFT):
        PerNameRiskPropertySensitivity<ExpiryWindow>(TYPE, shiftSize, NAME)
    {}

    PerNameRiskPropertySensitivity<ExpiryWindow>::Deriv deriv() const {
        return Deriv(IResultsFunction::price(),
                     ParSpreadMaxTweakSize::adapted(
                         RiskProperty<LegalBasisRelativePointwise>::SP()),
                     IScalarDerivative::oneSided(),
                     0.01);
    }
};

const double LegalBasisMultiplierPointwise::DEFAULT_SHIFT = 0.01; // Requirement: 1%

CClassConstSP const LegalBasisMultiplierPointwise::TYPE = CClass::registerClassLoadMethod(
    "LegalBasisMultiplierPointwise", typeid(LegalBasisMultiplierPointwise), load);

bool LegalBasisMultiplierPointwiseLinkIn() {
    return LegalBasisMultiplierPointwise::TYPE != NULL;
}

DRLIB_END_NAMESPACE

// Local Variables: ***
// c-basic-offset:4 ***
// End: ***
