//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : IRGridPoint.cpp
//
//   Description : Identifies a point on an IR vol surface
//
//   Author      : Andrew J Swain
//
//   Date        : 22 February 2002
//
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/IRGridPoint.hpp"

DRLIB_BEGIN_NAMESPACE

IRGridPoint::~IRGridPoint(){}

IRGridPoint::IRGridPoint(): CObject(TYPE) {}

IRGridPoint::IRGridPoint(const ExpiryConstSP& tenor, const ExpiryConstSP& expiry): 
    CObject(TYPE), tenor(tenor.clone()), expiry(expiry.clone()) {}


/** returns the expiry associated with this result */
ExpiryConstSP IRGridPoint::getExpiry() const{
    return expiry;
}

/** returns the tenor associated with this result */
ExpiryConstSP IRGridPoint::getTenor() const{
    return tenor;
}


class IRGridPointHelper{
public:
    /** Invoked when class is 'loaded' */
    static void load(CClassSP& clazz){
        clazz->setPublic(); // make visible to EAS/spreadsheet
        REGISTER(IRGridPoint, clazz);
        SUPERCLASS(CObject);
        EMPTY_SHELL_METHOD(defaultIRGridPoint);
        FIELD(tenor, "tenor");
        FIELD(expiry, "expiry");
    }
    static IObject* defaultIRGridPoint(){
        return new IRGridPoint();
    }
};

CClassConstSP const IRGridPoint::TYPE = CClass::registerClassLoadMethod(
    "IRGridPoint", typeid(IRGridPoint), IRGridPointHelper::load);

DEFINE_TEMPLATE_TYPE(IRGridPointArray);

DRLIB_END_NAMESPACE
