//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : VolBase.cpp
//
//   Description : Abstract vol interface
//
//   Author      : Mark A Robson
//
//   Date        : 13 Jan 2001
//
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#define QLIB_VOLBASE_CPP
#include "edginc/FXVol.hpp"

DRLIB_BEGIN_NAMESPACE

/** populate from market cache - default implementation provided */
void CVolBase::getMarket(const IModel* model, const MarketData* market) {
    // do nothing
}

/** Override the one in MarketObject to allow for automatic conversion
    of vols to fx vols */
void CVolBase::convert(IObjectSP&    object,
                       CClassConstSP requiredType) const{
    if (FXVolBaseWrapper::TYPE->isAssignableFrom(requiredType) &&
        !FXVolBase::TYPE->isInstance(this)){
        // have an ordinary [non fx] vol, want a fxvol wrapper
        // so create empty wrapper
        FXVolBaseWrapper& wrapper = 
            dynamic_cast<FXVolBaseWrapper&>(*requiredType->newInstance());
        object = IObjectSP(&wrapper); // return the wrapper as our output
        // create the fx vol
        MarketObjectSP fxVol(new FXVol(this));
        wrapper.setObject(fxVol);
        wrapper.setCacheUse(false);
    } else if (FXVolBase::TYPE->isAssignableFrom(requiredType)) {
        // have an ordinary [non fx] vol, want a fxvolbase
        object = IObjectSP(new FXVol(this)); // create an fx vol
    } else {
        // default to original method
        MarketObject::convert(object, requiredType);
    }

    // check we got what we asked for on the FX side
    if (!requiredType->isAssignableFrom(object->getClass())) {
        throw ModelException("CVolBase::convert",
                             "Required " + requiredType->getName() + 
                             " but conversion yielded " + 
                             object->getClass()->getName() + " instead"); 
    }
}

CVolBase::~CVolBase(){}

CVolBase::CVolBase(const CClassConstSP& clazz): 
    MarketObject(clazz){}

class CVolBaseHelper{
public:
    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        clazz->setPublic(); // make visible to EAS/spreadsheet
        REGISTER(CVolBase, clazz);
        SUPERCLASS(MarketObject);
        IMPLEMENTS(IGetMarket);
    }
};

CClassConstSP const CVolBase::TYPE = CClass::registerClassLoadMethod(
    "VolBase", typeid(CVolBase), CVolBaseHelper::load);

// definition of TYPE for MarketWrapper template class
DEFINE_TEMPLATE_TYPE_WITH_NAME("VolBaseWrapper", CVolBaseWrapper);

DRLIB_END_NAMESPACE
