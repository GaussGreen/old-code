/**
 * @file LegalBasisAdditivePointwise.cpp
 */

#include "edginc/config.hpp"
#include "edginc/DefaultConstructor.hpp"
#include "edginc/MultiTweakGroup.hpp"
#include "edginc/IResultsFunction.hpp"
#include "edginc/RiskProperty.hpp"
#include "edginc/IScalarDerivative.hpp"
#include "edginc/Additive.hpp"
#include "edginc/PerNameRiskPropertySensitivity.hpp"
#include "edginc/CreditTweak.hpp"
#include "edginc/LegalBasisAbsolutePointwise.hpp"
#include "edginc/GenericSensitivityFactory.hpp"
#include "edginc/ParSpreadMaxTweakSize.hpp"

DRLIB_BEGIN_NAMESPACE

static const char* NAME = "LEGAL_BASIS_ADDITIVE_POINTWISE";
const double ONE_BASIS_POINT = 0.0001;

/**
 * A greek calculated by tweaking each market name's par spreads at all its
 * defined expiries.
 *
 * LegalBasisAdditivePointwise has been rewritten to use the "declarative" sensitivities
 * framework: see IRiskQuantityFactory for an overview.
 */

class LegalBasisAdditivePointwise: public PerNameRiskPropertySensitivity<ExpiryWindow>,
                                   public virtual Additive,
                                   public CreditTweak {

    static void load(CClassSP& clazz) {
        clazz->setPublic();
        REGISTER(LegalBasisAdditivePointwise, clazz);
        IMPLEMENTS(Additive);
        SUPERCLASS(PerNameRiskPropertySensitivity<ExpiryWindow>);
        EMPTY_SHELL_METHOD(DefaultConstructor<LegalBasisAdditivePointwise>::iObject);
        SensitivityFactory::addSens(NAME,
                                    new GenericSensitivityFactory<LegalBasisAdditivePointwise>(), 
                                    new LegalBasisAdditivePointwise(DEFAULT_SHIFT),
                                    ITweakableWithRespectTo<LegalBasisAbsolutePointwise>::TYPE);
    }

public:

    static CClassConstSP const TYPE;

    static const double DEFAULT_SHIFT;

    LegalBasisAdditivePointwise(double shiftSize = DEFAULT_SHIFT):
        PerNameRiskPropertySensitivity<ExpiryWindow>(TYPE, shiftSize, NAME)
    {}

    PerNameRiskPropertySensitivity<ExpiryWindow>::Deriv deriv() const {
        return Deriv(IResultsFunction::price(),
                     ParSpreadMaxTweakSize::adapted(
                         RiskProperty<LegalBasisAbsolutePointwise>::SP()),
                     IScalarDerivative::oneSided(),
                     ONE_BASIS_POINT);
    }
};

const double LegalBasisAdditivePointwise::DEFAULT_SHIFT = 0.1 * ONE_BASIS_POINT;

CClassConstSP const LegalBasisAdditivePointwise::TYPE = CClass::registerClassLoadMethod(
    "LegalBasisAdditivePointwise", typeid(LegalBasisAdditivePointwise), load);

bool LegalBasisAdditivePointwiseLinkIn() {
    return LegalBasisAdditivePointwise::TYPE != NULL;
}

DRLIB_END_NAMESPACE

// Local Variables: ***
// c-basic-offset:4 ***
// End: ***
