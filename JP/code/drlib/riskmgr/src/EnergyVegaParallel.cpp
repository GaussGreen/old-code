//----------------------------------------------------------------------------
//
//   Group       : GCCG QR&D
//
//   Filename    : EnergyVega.cpp
//
//   Description : Energy Vega tweaking
//
//   Author      : Sean Chen
//
//   Date        : May 11 2005
//
//
//----------------------------------------------------------

/**
 * @file EnergyVegaParallel.cpp
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
#include "edginc/ScalarRiskPropertySensitivity.hpp"

#include "edginc/VolParallel.hpp"


DRLIB_BEGIN_NAMESPACE




class EnergyVegaParallel: public ScalarRiskPropertySensitivity,
                                     public virtual Additive 
{

public:

    static CClassConstSP const TYPE;
    static const string NAME;
	static const double SENSITIVITY_UNIT;
	static const double DEFAULT_SHIFT;


    EnergyVegaParallel(double shiftSize = DEFAULT_SHIFT) :
	             ScalarRiskPropertySensitivity(TYPE, NAME, shiftSize)
    {}

    ~EnergyVegaParallel() {}

private:

   
    ScalarRiskPropertySensitivity::Deriv deriv() const 
	{
		return Deriv(IResultsFunction::price(),
                     RiskProperty<VolParallel>::SP(),
                     IScalarDerivative::oneSided(),
                     SENSITIVITY_UNIT);
        
    }

    static void load(CClassSP& clazz) 
	{
        clazz->setPublic();
        REGISTER(EnergyVegaParallel, clazz);
        IMPLEMENTS(Additive);
        SUPERCLASS(ScalarRiskPropertySensitivity);
        EMPTY_SHELL_METHOD(&DefaultConstructor<EnergyVegaParallel>::iObject);
        SensitivityFactory::addSens(NAME,
                                    new GenericSensitivityFactory<EnergyVegaParallel>(), 
                                    new EnergyVegaParallel(),
                                    ITweakableWithRespectTo<VolParallel>::TYPE);
    }
};

const string EnergyVegaParallel::NAME = "ENERGY_VEGA_PARALLEL";
const double EnergyVegaParallel::SENSITIVITY_UNIT = 0.01;
const double EnergyVegaParallel::DEFAULT_SHIFT = 0.001;

CClassConstSP const EnergyVegaParallel::TYPE = 
                CClass::registerClassLoadMethod("EnergyVegaParallel", 
                                                typeid(EnergyVegaParallel), 
                                                load);
bool EnergyVegaParallelLinkIn() {
    return EnergyVegaParallel::TYPE != NULL;
}

DRLIB_END_NAMESPACE
