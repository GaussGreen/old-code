/**
 * @file MultiTweakGroup.cpp
 */

#include "edginc/config.hpp"
#include "edginc/IInstrumentCollection.hpp"
#include "edginc/MultiTweakGroup.hpp"

DRLIB_BEGIN_NAMESPACE

MultiTweakGroup::~MultiTweakGroup(){}

/** Constructor: takes references to inputs, not copies */
MultiTweakGroup::MultiTweakGroup(IInstrumentCollectionSP insts,
                                 IModelSP model): 
    CObject(TYPE), insts(insts), model(model){}

/** Constructor: takes references to inputs, not copies */
MultiTweakGroupSP MultiTweakGroup::SP(IInstrumentCollectionSP inst,
                                      IModelSP model) {
    return MultiTweakGroupSP(new MultiTweakGroup(inst, model));
}

/** Returns the instruments */
IInstrumentCollectionSP MultiTweakGroup::getInstruments() const{
    return insts;
}

/** Returns the IModel */
IModel* MultiTweakGroup::getModel() const{
    return model.get();
}
/** Returns the Model as SP */
IModelSP MultiTweakGroup::getModelSP() const{
    return model;
}

IObject* MultiTweakGroup::defaultMultiTweakGroup(){
    return new MultiTweakGroup();
}

MultiTweakGroup::MultiTweakGroup(): CObject(TYPE) {}

void MultiTweakGroup::load(CClassSP& clazz){
    clazz->setPublic();
    REGISTER(MultiTweakGroup, clazz);
    SUPERCLASS(CObject);
    EMPTY_SHELL_METHOD(defaultMultiTweakGroup);
    FIELD(insts, "instruments");
    FIELD(model, "algorithm");
}

CClassConstSP const MultiTweakGroup::TYPE = CClass::registerClassLoadMethod(
    "MultiTweakGroup", typeid(MultiTweakGroup), load);

DRLIB_END_NAMESPACE

// Local Variables: ***
// c-basic-offset:4 ***
// End: ***
