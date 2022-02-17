/**
 * @file LegalBasisAdditiveRelativeParallel.cpp
 */

#include "edginc/config.hpp"
#include "edginc/DefaultConstructor.hpp"
#include "edginc/MultiTweakGroup.hpp"
#include "edginc/IResultsFunction.hpp"
#include "edginc/RiskProperty.hpp"
#include "edginc/IScalarDerivative.hpp"
#include "edginc/ScalarRiskPropertySensitivity.hpp"
#include "edginc/LegalBasisAdditiveRelativeParallelTweak.hpp"
#include "edginc/GenericSensitivityFactory.hpp"
#include "edginc/ParSpreadMaxTweakSize.hpp"

DRLIB_BEGIN_NAMESPACE

static const char* NAME = "LEGAL_BASIS_ADDITIVE_RELATIVE_PARALLEL";


/**
 * Calculates derivative of instrument price w.r.t. relative shifts in the additive
 * coefficients of the Legal Basis adjustment done to the CDSParSpreadCurve
 */

class LegalBasisAdditiveRelativeParallel: public ScalarRiskPropertySensitivity {

    static void load(CClassSP& clazz) {
        clazz->setPublic();
        REGISTER(LegalBasisAdditiveRelativeParallel, clazz);
        SUPERCLASS(ScalarRiskPropertySensitivity);
        EMPTY_SHELL_METHOD(
            &DefaultConstructor<LegalBasisAdditiveRelativeParallel>::iObject);
        SensitivityFactory::addSens(
            NAME,
            new GenericSensitivityFactory<LegalBasisAdditiveRelativeParallel>(), 
            new LegalBasisAdditiveRelativeParallel(DEFAULT_SHIFT),
            ITweakableWithRespectTo<LegalBasisAdditiveRelativeParallelTweak>::TYPE);
    }

public:

    static CClassConstSP const TYPE;

    static const double DEFAULT_SHIFT;

    LegalBasisAdditiveRelativeParallel(double shiftSize = DEFAULT_SHIFT):
        ScalarRiskPropertySensitivity(TYPE, NAME, shiftSize)
    {}

    ScalarRiskPropertySensitivity::Deriv deriv() const {
        return Deriv(IResultsFunction::price(),
                     RiskProperty<LegalBasisAdditiveRelativeParallelTweak>::SP(),
                     IScalarDerivative::oneSided(),
                     0.01);
    }
};

const double LegalBasisAdditiveRelativeParallel::DEFAULT_SHIFT = 0.01; // Requirement: 1%

CClassConstSP const LegalBasisAdditiveRelativeParallel::TYPE = 
    CClass::registerClassLoadMethod("LegalBasisAdditiveRelativeParallel", 
                                    typeid(LegalBasisAdditiveRelativeParallel), 
                                    load);

bool LegalBasisAdditiveRelativeParallelLinkIn() {
    return LegalBasisAdditiveRelativeParallel::TYPE != NULL;
}

DRLIB_END_NAMESPACE

// Local Variables: ***
// c-basic-offset:4 ***
// End: ***
