//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : TweakGroup.cpp
//
//   Description : Holds objects to be tweaked for greeks - basically
//                 instrument and model
//
//   Author      : Mark A Robson
//
//   Date        : 9 May 2002
//
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#define QLIB_TWEAKGROUP_CPP
#include "edginc/TweakGroup.hpp"
DRLIB_BEGIN_NAMESPACE

/** Constructor: takes references to inputs, not copies */
TweakGroup::TweakGroup(const InstrumentSP& inst,
                       const IModelSP&      model): 
    CObject(TYPE), inst(inst), model(model){}

/** Returns the instrument */
Instrument* TweakGroup::getInstrument() const{
    return inst.get();
}

/** Returns the instrument as SP - todo: merge this with getInstrument */
InstrumentSP TweakGroup::getInstrumentSP() const{
    return inst;
}

/** Returns the IModel */
IModel* TweakGroup::getModel() const{
    return model.get();
}
/** Returns the Model as SP */
IModelSP TweakGroup::getModelSP() const{
    return model;
}

IObject* TweakGroup::defaultTweakGroup(){
    return new TweakGroup();
}

TweakGroup::TweakGroup(): CObject(TYPE){}

void TweakGroup::load(CClassSP& clazz){
    REGISTER(TweakGroup, clazz);
    SUPERCLASS(CObject);
    EMPTY_SHELL_METHOD(defaultTweakGroup);
    FIELD(inst, "instrument");
    FIELD(model, "algorithm");
}

CClassConstSP const TweakGroup::TYPE = CClass::registerClassLoadMethod(
    "TweakGroup", typeid(TweakGroup), load);

DRLIB_END_NAMESPACE

