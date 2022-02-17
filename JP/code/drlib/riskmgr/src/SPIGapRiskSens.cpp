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

/** Sensitivity of OutputRequest::SPI_GAP_RISK when tweaking the name's spot price  */
class RISKMGR_DLL SPIGapRiskDelta: public ScalarRiskPropertySensitivity,
                                            public virtual Additive {
    static void load(CClassSP& clazz) {
        clazz->setPublic();
        REGISTER(SPIGapRiskDelta, clazz);
        IMPLEMENTS(Additive);
        SUPERCLASS(ScalarRiskPropertySensitivity);
        EMPTY_SHELL_METHOD(&DefaultConstructor<SPIGapRiskDelta>::iObject);
        SensitivityFactory::addSens(NAME,
                                    new GenericSensitivityFactory<SPIGapRiskDelta>(), 
                                    new SPIGapRiskDelta(DEFAULT_SHIFT),
                                    ITweakableWithRespectTo<Spot>::TYPE);
    }

    Deriv deriv() const {
        return Deriv(IResultsFunction::outputRequest(OutputRequest::SPI_GAP_RISK),
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

    SPIGapRiskDelta(double shiftSize = DEFAULT_SHIFT): 
    ScalarRiskPropertySensitivity(TYPE, NAME, SECOND_ORDER_NAME, shiftSize) {}
    
    ~SPIGapRiskDelta() {};
};

const string SPIGapRiskDelta::NAME = "SPI_GAP_RISK_DELTA";
const string SPIGapRiskDelta::SECOND_ORDER_NAME = "SPI_GAP_RISK_GAMMA";
const double SPIGapRiskDelta::DEFAULT_SHIFT = Delta::DEFAULT_SHIFT;

CClassConstSP const SPIGapRiskDelta::TYPE = CClass::registerClassLoadMethod(
    "SPIGapRiskDelta", typeid(SPIGapRiskDelta), load);


////////////////////////////////////////////////////////////////////////////////////////////

/** Sensitivity of OutputRequest::SPI_GAP_RISK when tweaking the name's vol (parallel)  */
class RISKMGR_DLL SPIGapRiskVega: public ScalarRiskPropertySensitivity,
                                           public virtual Additive {
    static void load(CClassSP& clazz) {
        clazz->setPublic();
        REGISTER(SPIGapRiskVega, clazz);
        IMPLEMENTS(Additive);
        SUPERCLASS(ScalarRiskPropertySensitivity);
        EMPTY_SHELL_METHOD(&DefaultConstructor<SPIGapRiskVega>::iObject);
        SensitivityFactory::addSens(NAME,
                                    new GenericSensitivityFactory<SPIGapRiskVega>(), 
                                    new SPIGapRiskVega(DEFAULT_SHIFT),
                                    ITweakableWithRespectTo<VolParallel>::TYPE);
    }

    Deriv deriv() const {
        return Deriv(IResultsFunction::outputRequest(OutputRequest::SPI_GAP_RISK),
                     RiskProperty<VolParallel>::SP(),
                     IScalarDerivative::oneSided(),
                     SENSITIVITY_UNIT);
    }

public:
    static CClassConstSP const TYPE;

    static const double DEFAULT_SHIFT;
    static const double SENSITIVITY_UNIT;
    static const string NAME;

    SPIGapRiskVega(double shiftSize = DEFAULT_SHIFT): 
    ScalarRiskPropertySensitivity(TYPE, NAME, shiftSize) {}
    
    ~SPIGapRiskVega() {};
};

const string SPIGapRiskVega::NAME = "SPI_GAP_RISK_VEGA";
const double SPIGapRiskVega::SENSITIVITY_UNIT = VegaParallel::SENSITIVITY_UNIT;
const double SPIGapRiskVega::DEFAULT_SHIFT = VegaParallel::DEFAULT_SHIFT;

CClassConstSP const SPIGapRiskVega::TYPE = CClass::registerClassLoadMethod(
    "SPIGapRiskVega", typeid(SPIGapRiskVega), load);


////////////////////////////////////////////////////////////////////////////////////////////


bool SPIGapRiskSensLoad() {
    bool success = 
        SPIGapRiskDelta::TYPE != NULL && 
        SPIGapRiskVega::TYPE != NULL;
    return success;
}


DRLIB_END_NAMESPACE
