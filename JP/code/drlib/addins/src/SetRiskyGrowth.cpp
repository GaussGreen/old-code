/**
 * @file SetRiskyGrowth.cpp
 */

#include "edginc/config.hpp"
#include "edginc/XLConvert.hpp"
#include "edginc/ConvBond.hpp"

DRLIB_BEGIN_NAMESPACE

struct SetRiskyGrowth: public CObject {
    static CClassConstSP const TYPE;

    ConvBondSP          convBond;
    bool                 isRisky;

    static IObjectSP setRiskyGrowth(SetRiskyGrowth* params){
        ConvBondSP inst(copy(params->convBond.get()));
        inst->setRisky(params->isRisky);
        return inst;
    }
 
    SetRiskyGrowth(): CObject(TYPE){}
    
    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        REGISTER(SetRiskyGrowth, clazz);
        SUPERCLASS(CObject);
        EMPTY_SHELL_METHOD(defaultSetRiskyGrowth);
        FIELD(convBond,      "Convertible Bond");
        FIELD(isRisky, "if true, set to use risky growth.  false, set to not use risky growth.");
        Addin::registerClassObjectMethod("SET_USE_RISKY_GROWTH_FOR_CONVERT",
                                         Addin::UTILITIES,
                                         "Sets risky field in a convert",
                                         TYPE,
                                         true,
                                         Addin::returnHandle,
                                         (Addin::ObjMethod*)setRiskyGrowth);

    }
    
    static IObject* defaultSetRiskyGrowth(){
        return new SetRiskyGrowth();
    }
};

CClassConstSP const SetRiskyGrowth::TYPE= CClass::registerClassLoadMethod(
    "SetRiskyGrowth", typeid(SetRiskyGrowth), SetRiskyGrowth::load);

bool SetRiskyGrowthLoad() {
    return SetRiskyGrowth::TYPE != NULL;
}

DRLIB_END_NAMESPACE
