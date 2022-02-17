/**
 * @file EnergyDeltaPointwiseWithShift.cpp
 */

#include "edginc/config.hpp"
#include "edginc/Maths.hpp"
#include "edginc/Expiry.hpp"
#include "edginc/ExpiryWindow.hpp"
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
#include "edginc/GenericSensitivityFactory.hpp"
#include "edginc/DefaultConstructor.hpp"
#include "edginc/PerNameRiskPropertySensitivity.hpp"
#include "edginc/EnergyFuturesCurvePointwise.hpp"
#include "edginc/EnergyFuturesCurveParallel.hpp"


DRLIB_BEGIN_NAMESPACE



class EnergyDeltaPointwiseWithShift: public PerNameRiskPropertySensitivity<ExpiryWindow>,
                                     public virtual Additive 
{

public:
    static CClassConstSP const TYPE;

    static const string NAME;
    static const double DEFAULT_SHIFT;

    EnergyDeltaPointwiseWithShift(double shiftSize = DEFAULT_SHIFT,
                                  double parallelShiftSize = 0.020 /* 10% */):
	    PerNameRiskPropertySensitivity<ExpiryWindow>(TYPE, shiftSize, NAME, "ENERGY_GAMMA_POINTWISE"),
        //VectorRiskPropertySensitivity(TYPE, NAME, shiftSize), // added a ENERGY_GAMMA_POINTWISE
        parallelShiftSize(parallelShiftSize)
    {}

    ~EnergyDeltaPointwiseWithShift() {}

private:

    double parallelShiftSize;

    PerNameRiskPropertySensitivity<ExpiryWindow>::Deriv deriv() const {
		return Deriv(IResultsFunction::price(),
                     RiskProperty<EnergyFuturesCurvePointwise>::SP(),
                     IScalarDerivative::twoSided()->averagedOverScenarios(
                         RiskProperty<EnergyFuturesCurveParallel>::SP(),
                         parallelShiftSize, -parallelShiftSize),
                     1.,// adding for ENERGY_GAMMA_POINTWISE
					 IScalarDerivative::second()->averagedOverScenarios(
                         RiskProperty<EnergyFuturesCurveParallel>::SP(),
                         parallelShiftSize, -parallelShiftSize),
                     1.);
        
    }

    static void load(CClassSP& clazz) {
        clazz->setPublic();
        REGISTER(EnergyDeltaPointwiseWithShift, clazz);
        IMPLEMENTS(Additive);
        SUPERCLASS(PerNameRiskPropertySensitivity<ExpiryWindow>);
        EMPTY_SHELL_METHOD(&DefaultConstructor<EnergyDeltaPointwiseWithShift>::iObject);
        FIELD(parallelShiftSize, "Shift size for parallel scenario");
        FIELD_MAKE_OPTIONAL(parallelShiftSize);
        SensitivityFactory::addSens(NAME,
                                    new GenericSensitivityFactory<EnergyDeltaPointwiseWithShift>(), 
                                    new EnergyDeltaPointwiseWithShift(),
                                    ITweakableWithRespectTo<EnergyFuturesCurvePointwise>::TYPE);
    }
};

const string EnergyDeltaPointwiseWithShift::NAME = "ENERGY_DELTA_POINTWISE";
const double EnergyDeltaPointwiseWithShift::DEFAULT_SHIFT = 0.1;

CClassConstSP const EnergyDeltaPointwiseWithShift::TYPE = 
    CClass::registerClassLoadMethod("EnergyDeltaPointwiseWithShift", 
                                    typeid(EnergyDeltaPointwiseWithShift), 
                                    load);

bool EnergyDeltaPointwiseWithShiftLinkIn() {
    return EnergyDeltaPointwiseWithShift::TYPE != NULL;
}

DRLIB_END_NAMESPACE
