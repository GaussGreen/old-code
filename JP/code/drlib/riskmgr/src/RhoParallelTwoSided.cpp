/**
 * @file RhoParallelTwoSided.cpp
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
#include "edginc/RateParallel.hpp"
#include "edginc/RhoParallel.hpp"
#include "edginc/IResultsIdentifier.hpp"
#include "edginc/Additive.hpp"
#include "edginc/ScalarRiskPropertySensitivity.hpp"
#include "edginc/SensitivityFactory.hpp"
#include "edginc/GenericSensitivityFactory.hpp"

DRLIB_BEGIN_NAMESPACE

FORWARD_DECLARE(RhoParallelTwoSided)

class RhoParallelTwoSided: public ScalarRiskPropertySensitivity,
                    public virtual Additive {

    static void load(CClassSP& clazz) {
        clazz->setPublic();
        REGISTER(RhoParallelTwoSided, clazz);
        IMPLEMENTS(Additive);
        SUPERCLASS(ScalarRiskPropertySensitivity);
        EMPTY_SHELL_METHOD(&DefaultConstructor<RhoParallelTwoSided>::iObject);
        SensitivityFactory::addSens(NAME,
                                    new GenericSensitivityFactory<RhoParallelTwoSided>(), 
                                    new RhoParallelTwoSided(DEFAULT_SHIFT),
                                    ITweakableWithRespectTo<RateParallel>::TYPE);
        SensitivityFactory::addSens(NAME2,
                                    new GenericSensitivityFactory<RhoParallelTwoSided>(), 
                                    new RhoParallelTwoSided(DEFAULT_SHIFT),
                                    ITweakableWithRespectTo<RateParallel>::TYPE);
    }

    Deriv deriv() const {
        return Deriv(IResultsFunction::price(),
                     RiskProperty<RateParallel>::SP(),
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

    RhoParallelTwoSided(double shiftSize = DEFAULT_SHIFT):
        ScalarRiskPropertySensitivity(TYPE, NAME, NAME2, shiftSize)
    {}

    ~RhoParallelTwoSided() {}
};

const string RhoParallelTwoSided::NAME = "RHO_PARALLEL_TWO_SIDED";
const string RhoParallelTwoSided::NAME2 = "IR_GAMMA_PARALLEL";
const double RhoParallelTwoSided::DEFAULT_SHIFT = RhoParallel::DEFAULT_SHIFT;
const double RhoParallelTwoSided::SENSITIVITY_UNIT = RhoParallel::SENSITIVITY_UNIT;

CClassConstSP const RhoParallelTwoSided::TYPE = CClass::registerClassLoadMethod(
    "RhoParallelTwoSided", typeid(RhoParallelTwoSided), load);

bool RhoParallelTwoSidedLinkIn() {
    return RhoParallelTwoSided::TYPE != NULL;
}

DRLIB_END_NAMESPACE

// Local Variables: ***
// c-basic-offset:4 ***
// End: ***
