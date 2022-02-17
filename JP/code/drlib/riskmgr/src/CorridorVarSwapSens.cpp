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
#include "edginc/IResultsIdentifier.hpp"
#include "edginc/NakedBondRhoParallel.hpp"
#include "edginc/SensitivityFactory.hpp"
#include "edginc/GenericSensitivityFactory.hpp"
#include "edginc/OutputRequest.hpp"
#include "edginc/Delta.hpp"
#include "edginc/VegaParallel.hpp"
#include "edginc/VolParallel.hpp"

DRLIB_BEGIN_NAMESPACE

/** Sensitivity of OutputRequest::CORRIDOR_VARIANCE_VALUE when tweaking the name's spot price  */
class RISKMGR_DLL CorridorVarianceLegDelta: public ScalarRiskPropertySensitivity,
                                            public virtual Additive {
    static void load(CClassSP& clazz) {
        clazz->setPublic();
        REGISTER(CorridorVarianceLegDelta, clazz);
        IMPLEMENTS(Additive);
        SUPERCLASS(ScalarRiskPropertySensitivity);
        EMPTY_SHELL_METHOD(&DefaultConstructor<CorridorVarianceLegDelta>::iObject);
        SensitivityFactory::addSens(NAME,
                                    new GenericSensitivityFactory<CorridorVarianceLegDelta>(), 
                                    new CorridorVarianceLegDelta(DEFAULT_SHIFT),
                                    ITweakableWithRespectTo<Spot>::TYPE);
    }

    Deriv deriv() const {
        return Deriv(IResultsFunction::outputRequest(OutputRequest::CORRIDOR_VARIANCE_VALUE),
                     RiskProperty<Spot>::SP(),
                     IScalarDerivative::twoSided(),
                     1.0,
                     IScalarDerivative::second(),
                     1.0);
    }

public:
    static CClassConstSP const TYPE;

    static const double DEFAULT_SHIFT;
    static const string NAME;
    static const string SECOND_ORDER_NAME;

    CorridorVarianceLegDelta(double shiftSize = DEFAULT_SHIFT): 
    ScalarRiskPropertySensitivity(TYPE, NAME, SECOND_ORDER_NAME, shiftSize) {}
    
    ~CorridorVarianceLegDelta() {};
};

const string CorridorVarianceLegDelta::NAME = "CORRIDOR_VARIANCE_DELTA";
const string CorridorVarianceLegDelta::SECOND_ORDER_NAME = "CORRIDOR_VARIANCE_GAMMA";
const double CorridorVarianceLegDelta::DEFAULT_SHIFT = Delta::DEFAULT_SHIFT;

CClassConstSP const CorridorVarianceLegDelta::TYPE = CClass::registerClassLoadMethod(
    "CorridorVarianceLegDelta", typeid(CorridorVarianceLegDelta), load);


////////////////////////////////////////////////////////////////////////////////////////////


/** Sensitivity of OutputRequest::CORRIDOR_ACCRUAL_VALUE when tweaking the name's spot price  */
class RISKMGR_DLL CorridorAccrualLegDelta: public ScalarRiskPropertySensitivity,
                                           public virtual Additive {
    static void load(CClassSP& clazz) {
        clazz->setPublic();
        REGISTER(CorridorAccrualLegDelta, clazz);
        IMPLEMENTS(Additive);
        SUPERCLASS(ScalarRiskPropertySensitivity);
        EMPTY_SHELL_METHOD(&DefaultConstructor<CorridorAccrualLegDelta>::iObject);
        SensitivityFactory::addSens(NAME,
                                    new GenericSensitivityFactory<CorridorAccrualLegDelta>(), 
                                    new CorridorAccrualLegDelta(DEFAULT_SHIFT),
                                    ITweakableWithRespectTo<Spot>::TYPE);
    }

    Deriv deriv() const {
        return Deriv(IResultsFunction::outputRequest(OutputRequest::CORRIDOR_ACCRUAL_VALUE),
                     RiskProperty<Spot>::SP(),
                     IScalarDerivative::twoSided(),
                     1.0,
                     IScalarDerivative::second(),
                     1.0);
    }

public:
    static CClassConstSP const TYPE;

    static const double DEFAULT_SHIFT;
    static const string NAME;
    static const string SECOND_ORDER_NAME;

    CorridorAccrualLegDelta(double shiftSize = DEFAULT_SHIFT): 
    ScalarRiskPropertySensitivity(TYPE, NAME, SECOND_ORDER_NAME, shiftSize) {}
    
    ~CorridorAccrualLegDelta() {};
};

const string CorridorAccrualLegDelta::NAME = "CORRIDOR_ACCRUAL_DELTA";
const string CorridorAccrualLegDelta::SECOND_ORDER_NAME = "CORRIDOR_ACCRUAL_GAMMA";
const double CorridorAccrualLegDelta::DEFAULT_SHIFT = Delta::DEFAULT_SHIFT;

CClassConstSP const CorridorAccrualLegDelta::TYPE = CClass::registerClassLoadMethod(
    "CorridorAccrualLegDelta", typeid(CorridorAccrualLegDelta), load);


////////////////////////////////////////////////////////////////////////////////////////////


/** Sensitivity of OutputRequest::CORRIDOR_VARIANCE_VALUE when tweaking the name's vol (parallel)  */
class RISKMGR_DLL CorridorVarianceLegVega: public ScalarRiskPropertySensitivity,
                                           public virtual Additive {
    static void load(CClassSP& clazz) {
        clazz->setPublic();
        REGISTER(CorridorVarianceLegVega, clazz);
        IMPLEMENTS(Additive);
        SUPERCLASS(ScalarRiskPropertySensitivity);
        EMPTY_SHELL_METHOD(&DefaultConstructor<CorridorVarianceLegVega>::iObject);
        SensitivityFactory::addSens(NAME,
                                    new GenericSensitivityFactory<CorridorVarianceLegVega>(), 
                                    new CorridorVarianceLegVega(DEFAULT_SHIFT),
                                    ITweakableWithRespectTo<VolParallel>::TYPE);
    }

    Deriv deriv() const {
        return Deriv(IResultsFunction::outputRequest(OutputRequest::CORRIDOR_VARIANCE_VALUE),
                     RiskProperty<VolParallel>::SP(),
                     IScalarDerivative::oneSided(),
                     SENSITIVITY_UNIT);
    }

public:
    static CClassConstSP const TYPE;

    static const double DEFAULT_SHIFT;
    static const double SENSITIVITY_UNIT;
    static const string NAME;

    CorridorVarianceLegVega(double shiftSize = DEFAULT_SHIFT): 
    ScalarRiskPropertySensitivity(TYPE, NAME, shiftSize) {}
    
    ~CorridorVarianceLegVega() {};
};

const string CorridorVarianceLegVega::NAME = "CORRIDOR_VARIANCE_VEGA";
const double CorridorVarianceLegVega::SENSITIVITY_UNIT = VegaParallel::SENSITIVITY_UNIT;
const double CorridorVarianceLegVega::DEFAULT_SHIFT = VegaParallel::DEFAULT_SHIFT;

CClassConstSP const CorridorVarianceLegVega::TYPE = CClass::registerClassLoadMethod(
    "CorridorVarianceLegVega", typeid(CorridorVarianceLegVega), load);


////////////////////////////////////////////////////////////////////////////////////////////


/** Sensitivity of OutputRequest::CORRIDOR_VARIANCE_VALUE when tweaking the name's vol (parallel)  */
class RISKMGR_DLL CorridorAccrualLegVega: public ScalarRiskPropertySensitivity,
                                          public virtual Additive {
    static void load(CClassSP& clazz) {
        clazz->setPublic();
        REGISTER(CorridorAccrualLegVega, clazz);
        IMPLEMENTS(Additive);
        SUPERCLASS(ScalarRiskPropertySensitivity);
        EMPTY_SHELL_METHOD(&DefaultConstructor<CorridorAccrualLegVega>::iObject);
        SensitivityFactory::addSens(NAME,
                                    new GenericSensitivityFactory<CorridorAccrualLegVega>(), 
                                    new CorridorAccrualLegVega(DEFAULT_SHIFT),
                                    ITweakableWithRespectTo<VolParallel>::TYPE);
    }

    Deriv deriv() const {
        return Deriv(IResultsFunction::outputRequest(OutputRequest::CORRIDOR_ACCRUAL_VALUE),
                     RiskProperty<VolParallel>::SP(),
                     IScalarDerivative::oneSided(),
                     SENSITIVITY_UNIT);
    }

public:
    static CClassConstSP const TYPE;

    static const double DEFAULT_SHIFT;
    static const double SENSITIVITY_UNIT;
    static const string NAME;

    CorridorAccrualLegVega(double shiftSize = DEFAULT_SHIFT): 
    ScalarRiskPropertySensitivity(TYPE, NAME, shiftSize) {}
    
    ~CorridorAccrualLegVega() {};
};

const string CorridorAccrualLegVega::NAME = "CORRIDOR_ACCRUAL_VEGA";
const double CorridorAccrualLegVega::SENSITIVITY_UNIT = VegaParallel::SENSITIVITY_UNIT;
const double CorridorAccrualLegVega::DEFAULT_SHIFT = VegaParallel::DEFAULT_SHIFT;

CClassConstSP const CorridorAccrualLegVega::TYPE = CClass::registerClassLoadMethod(
    "CorridorAccrualLegVega", typeid(CorridorAccrualLegVega), load);


////////////////////////////////////////////////////////////////////////////////////////////


bool CorridorVarSwapSensLinkIn() {
    bool success = 
        CorridorVarianceLegDelta::TYPE != NULL && 
        CorridorAccrualLegDelta::TYPE != NULL &&
        CorridorVarianceLegVega::TYPE != NULL && 
        CorridorAccrualLegVega::TYPE != NULL;
    return success;
}


DRLIB_END_NAMESPACE
