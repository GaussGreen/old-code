//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : CorrelationSkew.cpp
//
//   Description : Holds implied correlation parameters
//
//   Author      : Oliver Brockhaus
//
//   Date        : 15 Oct 2002
//
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#define QLIB_CORRELATIONSKEW_CPP
#include "edginc/CorrelationSkew.hpp"

DRLIB_BEGIN_NAMESPACE

CorrelationSkew::~CorrelationSkew(){}

/** Validation */
void CorrelationSkew::validatePop2Object(){
    static const string method("CorrelationSkew::validatePop2Object");
    if ( (correlationSkew >= 1.0) || (correlationSkew < 0.0) ){
        throw ModelException( method, "CorrelationSkew "+name+
                              " must be positive and less than 1.0.");
    }
    if ( correlationSkewPower < 0.0 ){
        throw ModelException( method, "CorrelationSkewPower "+name+
                              " must be positive.");
    }
}

double CorrelationSkew::getCorrelationSkew() const
{
    return correlationSkew;
}

double CorrelationSkew::getCorrelationSkewPower() const
{
    return correlationSkewPower;
}

// Sensitivity methods

/** Returns the name of the protected equity - used to determine
    whether to tweak the object */
string CorrelationSkew::sensName(PhiSkew* shift)const
{
    return getName();
}

string CorrelationSkew::sensName(PhiSkewPower* shift)const
{
    return getName();
}

/** Shifts the object using given shift */    
bool CorrelationSkew::sensShift(PhiSkew* shift)
{
    // must set the intitial value
    shift->setInitialValue(correlationSkew);
    
    double shiftSize = shift->getShiftSize();
 
    // then shift
    if (correlationSkew > 0.5)
    {
        correlationSkew -= shiftSize;
    }
    else
    {
        correlationSkew += shiftSize; 
    }
    return false; // none of our components has a PHI_SKEW type sensitivity
}

bool CorrelationSkew::sensShift(PhiSkewPower* shift)
{
    // must set the intitial value
    shift->setInitialValue(correlationSkewPower);
    
    double shiftSize = shift->getShiftSize();
 
    // then shift
    if (correlationSkewPower > 0.5)
    {
        correlationSkewPower -= shiftSize;
    }
    else
    {
        correlationSkewPower += shiftSize; 
    }
    return false; // none of our components has a PHI_SKEW_POWER type sensitivity
}

/** Restores the object to its original form */
void CorrelationSkew::sensRestore(PhiSkew* shift)
{
    correlationSkew = shift->getInitialValue();
}

void CorrelationSkew::sensRestore(PhiSkewPower* shift)
{
    correlationSkewPower = shift->getInitialValue();
}

string CorrelationSkew::getName() const
{
    return name;
}

CorrelationSkew::CorrelationSkew(const string& name,
                                 double        correlationSkew,
                                 double        correlationSkewPower):MarketObject(TYPE), 
    name(name), 
    correlationSkew(correlationSkew),
    correlationSkewPower(correlationSkewPower) {
    validatePop2Object();
}

CorrelationSkew::CorrelationSkew():MarketObject(TYPE)
{
    this->correlationSkew       = 0.0;
    this->correlationSkewPower  = 0.0;
}

class CorrelationSkewHelper{
public:
    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        clazz->setPublic(); // make visible to EAS/spreadsheet
        REGISTER(CorrelationSkew, clazz);
        SUPERCLASS(MarketObject);
        IMPLEMENTS(PhiSkew::RestorableShift);
        IMPLEMENTS(PhiSkewPower::RestorableShift);
        EMPTY_SHELL_METHOD(defaultCorrelationSkew);
        FIELD(name, "name for correlation skew");
        FIELD(correlationSkew, "correlation skew value ");
        FIELD(correlationSkewPower, "correlation skew power ");
    }

    static IObject* defaultCorrelationSkew(){
        return new CorrelationSkew();
    }
};

CClassConstSP const CorrelationSkew::TYPE = CClass::registerClassLoadMethod(
    "CorrelationSkew", typeid(CorrelationSkew), CorrelationSkewHelper::load);

//support for arrays
DEFINE_TEMPLATE_TYPE(CorrelationSkewArray);

///// support for wrapper class /////

// definition of TYPE for MarketWrapper template class
DEFINE_TEMPLATE_TYPE(CorrelationSkewWrapper);

/** Casts array element to an IObject */
IObjectConstSP arrayObjectCast<CorrelationSkewWrapper>::toIObject(
    const CorrelationSkewWrapper& value){
    IObjectConstSP objValue(IObjectConstSP::attachToRef(&value));
    return objValue;
}

/** Casts array element to an IObject */
IObjectSP arrayObjectCast<CorrelationSkewWrapper>::toIObject(CorrelationSkewWrapper& value){
    return IObjectSP::attachToRef(&value);
}

/** Turns the IObjectSP into a DateTime */
CorrelationSkewWrapper arrayObjectCast<CorrelationSkewWrapper>::fromIObject(IObjectSP& value){
    CorrelationSkewWrapper *ptr = dynamic_cast<CorrelationSkewWrapper *>(value.get());
    if (!ptr){
        throw ModelException("arrayObjectCast::fromIObject", "Object is not"
                             " a CorrelationSkewWrapper but a "+
                             value->getClass()->getName());
    }
    return *ptr;
}

DEFINE_TEMPLATE_TYPE(CorrelationSkewWrapperArray);

bool CorrelationSkewLoad() {
    return (CorrelationSkew::TYPE != 0);
}

DRLIB_END_NAMESPACE
