/**
 * @file RhoPointwiseTwoSided.cpp
 */

#include "edginc/config.hpp"
#include "edginc/Maths.hpp"
#include "edginc/DefaultConstructor.hpp"
#include "edginc/MultiTweakGroup.hpp"
#include "edginc/IInstrumentCollection.hpp"
#include "edginc/IResultsFunction.hpp"
#include "edginc/PropertyRiskAxis.hpp"
#include "edginc/RiskProperty.hpp"
#include "edginc/PropertyTweakHypothesis.hpp"
#include "edginc/RiskQuantity.hpp"
#include "edginc/IScalarDerivative.hpp"
#include "edginc/NamedRiskQuantity.hpp"
#include "edginc/RatePointwise.hpp"
#include "edginc/RhoPointwise.hpp"
#include "edginc/IResultsIdentifier.hpp"
#include "edginc/Additive.hpp"
#include "edginc/PerNameRiskPropertySensitivity.hpp"
#include "edginc/SensitivityFactory.hpp"
#include "edginc/GenericSensitivityFactory.hpp"

DRLIB_BEGIN_NAMESPACE

FORWARD_DECLARE(RhoPointwiseTwoSided)

class RhoPointwiseTwoSided: public PerNameRiskPropertySensitivity<ExpiryWindow>,
                            public virtual Additive {
    
    static void load(CClassSP& clazz) {
        clazz->setPublic();
        REGISTER(RhoPointwiseTwoSided, clazz);
        IMPLEMENTS(Additive);
        SUPERCLASS(PerNameRiskPropertySensitivity<ExpiryWindow>);
        EMPTY_SHELL_METHOD(&DefaultConstructor<RhoPointwiseTwoSided>::iObject);
        SensitivityFactory::addSens(NAME,
                                    new GenericSensitivityFactory<RhoPointwiseTwoSided>(), 
                                    new RhoPointwiseTwoSided(DEFAULT_SHIFT),
                                    ITweakableWithRespectTo<RatePointwise>::TYPE);
        SensitivityFactory::addSens(NAME2,
                                    new GenericSensitivityFactory<RhoPointwiseTwoSided>(), 
                                    new RhoPointwiseTwoSided(DEFAULT_SHIFT),
                                    ITweakableWithRespectTo<RatePointwise>::TYPE);
    }

    Deriv deriv() const {
        return Deriv(IResultsFunction::price(),
                     RiskProperty<RatePointwise>::SP(),
                     IScalarDerivative::twoSided(),
                     SENSITIVITY_UNIT,
                     IScalarDerivative::second(),
                     SENSITIVITY_UNIT);
    }

public:

    static CClassConstSP const TYPE;
    
    static const double DEFAULT_SHIFT;
    static const double SENSITIVITY_UNIT;
    static const string NAME;
    static const string NAME2;

    RhoPointwiseTwoSided(double shiftSize = DEFAULT_SHIFT):
        PerNameRiskPropertySensitivity<ExpiryWindow>(TYPE, shiftSize, NAME, NAME2)
    {}

    ~RhoPointwiseTwoSided() {}
};

const string RhoPointwiseTwoSided::NAME = "RHO_POINTWISE_TWO_SIDED";
const string RhoPointwiseTwoSided::NAME2 = "IR_GAMMA_POINTWISE";
const double RhoPointwiseTwoSided::DEFAULT_SHIFT = RhoPointwise::DEFAULT_SHIFT;
const double RhoPointwiseTwoSided::SENSITIVITY_UNIT = RhoPointwise::SENSITIVITY_UNIT;

CClassConstSP const RhoPointwiseTwoSided::TYPE = CClass::registerClassLoadMethod(
    "RhoPointwiseTwoSided", typeid(RhoPointwiseTwoSided), load);

bool RhoPointwiseTwoSidedLinkIn() {
    return RhoPointwiseTwoSided::TYPE != NULL;
}

DRLIB_END_NAMESPACE

// Local Variables: ***
// c-basic-offset:4 ***
// End: ***
