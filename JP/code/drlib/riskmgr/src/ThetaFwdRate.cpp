//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : ThetaFwdRate.cpp
//
//   Description : Theta Forward Rate - preserve forward zero rates (df's)
//
//   Author      : Ning Shen
//
//   Date        : 5 Nov 2002
//
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/DateTime.hpp"
#include "edginc/ThetaFwdRate.hpp"
#include "edginc/Holiday.hpp"
#include "edginc/SensitivityFactory.hpp"

DRLIB_BEGIN_NAMESPACE

/** Sens Control for theta */
const string ThetaFwdRate::NAME = "THETA_FWD_RATE";

/** for reflection */
ThetaFwdRate::ThetaFwdRate(): Theta(TYPE, NAME){}

ThetaFwdRate::ThetaFwdRate(int offset, HolidaySP hols):
Theta(offset, hols, TYPE, NAME){}

/** Factory class dictates what methods of building this
    sensitivity are supported. Here we support a default */
class ThetaFwdRateFactory: public SensitivityFactory::IDefault {
public:
    virtual Sensitivity* createDefault(){
        HolidaySP hols(Holiday::weekendsOnly());
        return new ThetaFwdRate(Theta::DEFAULT_SHIFT, hols);
    }
};

/** Invoked when ThetaFwdRate is 'loaded' */
void ThetaFwdRate::load(CClassSP& clazz){
    clazz->setPublic(); // make visible to EAS/spreadsheet
    REGISTER(ThetaFwdRate, clazz);
    SUPERCLASS(Theta);
    EMPTY_SHELL_METHOD(defaultThetaFwdRate);
    // register how to build our sensitivity
    SensitivityFactory::addSens(ThetaFwdRate::NAME, 
                                new ThetaFwdRateFactory(), 
                                new ThetaFwdRate(),
                                ThetaFwdRate::Shift::TYPE);
}

IObject* ThetaFwdRate::defaultThetaFwdRate(){
    return new ThetaFwdRate();
}

CClassConstSP const ThetaFwdRate::TYPE = CClass::registerClassLoadMethod(
    "ThetaFwdRate", typeid(ThetaFwdRate), load);


DRLIB_END_NAMESPACE
