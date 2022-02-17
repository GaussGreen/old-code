//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : ThetaFwdSpot.cpp
//
//   Description : Theta Forward Spot shift
//
//   Author      : Mark A Robson
//
//   Date        : 26 June 2001
//
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/Holiday.hpp"
#include "edginc/ThetaFwdSpot.hpp"
#include "edginc/SensitivityFactory.hpp"

DRLIB_BEGIN_NAMESPACE

/** Sens Control for theta */
const string ThetaFwdSpot::NAME = "THETA_FWD_SPOT";


/** for reflection */
ThetaFwdSpot::ThetaFwdSpot(): Theta(TYPE, NAME){}

ThetaFwdSpot::ThetaFwdSpot(int offset, HolidaySP hols):
              Theta(offset, hols, TYPE, NAME) {}

/** Returns true */
bool ThetaFwdSpot::useAssetFwds() const{
    return true;
}

/** Factory class dictates what methods of building this
    sensitivity are supported. Here we support a default */
class TFSFactory: public SensitivityFactory::IDefault {
public:
    virtual Sensitivity* createDefault(){
        HolidaySP hols(Holiday::weekendsOnly());
        return new ThetaFwdSpot(Theta::DEFAULT_SHIFT, hols);
    }
};

/** Invoked when ThetaFwdSpot is 'loaded' */
void ThetaFwdSpot::load(CClassSP& clazz){
    clazz->setPublic(); // make visible to EAS/spreadsheet
    REGISTER(ThetaFwdSpot, clazz);
    SUPERCLASS(Theta);
    EMPTY_SHELL_METHOD(defaultThetaFwdSpot);
    // register how to build our sensitivity
    SensitivityFactory::addSens(ThetaFwdSpot::NAME, 
                                new TFSFactory(), 
                                new ThetaFwdSpot(),
                                ThetaFwdSpot::Shift::TYPE);
}

IObject* ThetaFwdSpot::defaultThetaFwdSpot(){
    return new ThetaFwdSpot();
}

CClassConstSP const ThetaFwdSpot::TYPE = CClass::registerClassLoadMethod(
    "ThetaFwdSpot", typeid(ThetaFwdSpot), load);


DRLIB_END_NAMESPACE
