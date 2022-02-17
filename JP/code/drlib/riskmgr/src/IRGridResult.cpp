//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : IRGridResult.cpp
//
//   Description : Result for pointwise type IR vega sensitivities
//
//   Author      : Andrew J Swain
//
//   Date        : 26 February 2002
//
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/IRGridResult.hpp"


DRLIB_BEGIN_NAMESPACE

/** Captures result of a result at a given expiry, ie stores an expiry
    and a double */

IRGridResult::~IRGridResult(){}

/* for simplicity allow array of IRGridResults to be an array of
   structures rather than an array of pointers - dictates public
   default constructor */
IRGridResult::IRGridResult(): CObject(TYPE), result(0){}

IRGridResult::IRGridResult(IRGridPointAbsConstSP point, double result): 
    CObject(TYPE), point(point->toGridPoint()), result(result) {}

IRGridResult::IRGridResult(IRGridPointSP point, double result): 
    CObject(TYPE), point(point), result(result) {}

/** returns the grid point associated with this result */
IRGridPointSP IRGridResult::getGridPoint() const{
    return point;
}

/** returns the double detailing this result */
double IRGridResult::getResult() const{
    return result;
}

/** write object out in 'output' format - ie suitable for comparing
    regression files with */
void IRGridResult::outputWrite(const string& linePrefix,
                               const string& prefix,
                               ostream&      stream) const{
    ios::fmtflags oldFlags     = stream.flags();     // save settings
    long          oldPrecision = stream.precision(); // save settings
    if (result < 1 && result > -1){
        // 8 dp
        stream.flags(ios::fixed);
    } else {
        // 8 sf
        stream.unsetf(ios::fixed);
    }        
    stream.precision(8);

    stream << linePrefix << prefix << "_" << point->getTenor()->toString() 
           << "_" << point->getExpiry()->toString() << ": " << result << endl;
    stream.flags(oldFlags); // restore settings
    stream.precision(oldPrecision);
}

/** scale by factor x */
void IRGridResult::scale(double x) {
    result *= x;
}

/** add IRGridResult object to this result (Implementation of
    CombinableResult) */
void IRGridResult::add(const CombinableResult& x, double scaleFactor) {
    // gcc bug: force to IObject before dynamic cast
    const IRGridResult& ex = dynamic_cast<const IRGridResult&>(static_cast<const IObject&>(x));
    if (!point->getExpiry()->equals(ex.point->getExpiry().get()) ||
        !point->getTenor()->equals(ex.point->getTenor().get())) {
        // shouldn't happen really
        throw ModelException("IRGridResult::add", 
                             "can't add [" + 
                             ex.point->getExpiry()->toString() + "," + 
                             ex.point->getTenor()->toString() + "] to [" +
                             point->getExpiry()->toString() + "," + 
                             point->getTenor()->toString() + "]");
    }
    result += scaleFactor * ex.result;           
}

/** adds the supplied value to this object */
void IRGridResult::addToResult(double value){
    result += value;
}

/** specialisations of arrayObjectCast */
/** Casts array element to an IObject */
IObjectConstSP arrayObjectCast<IRGridResult>::toIObject(
    const IRGridResult& value){
    IObjectConstSP objValue(IObjectConstSP::attachToRef(&value));
    return objValue;
}

/** Casts array element to an IObject */
IObjectSP arrayObjectCast<IRGridResult>::toIObject(IRGridResult& value){
    return IObjectSP::attachToRef(&value);
}

/** Turns the IObjectSP into a IRGridResult */
IRGridResult arrayObjectCast<IRGridResult>::fromIObject(IObjectSP& value){
    IRGridResult *ptr = dynamic_cast<IRGridResult *>(value.get());
    if (!ptr){
        throw ModelException("arrayObjectCast::fromIObject", "Object is not"
                             " an IRGridResult");
    }
    return *ptr;
}

class IRGridResultHelper{
public:
    /** Invoked when class is 'loaded' */
    static void load(CClassSP& clazz){
        clazz->setPublic(); // make visible to EAS/spreadsheet
        REGISTER(IRGridResult, clazz);
        SUPERCLASS(CObject);
        IMPLEMENTS(CombinableResult);
        EMPTY_SHELL_METHOD(defaultIRGridResult);
        FIELD(point, "Indicates which grid point was tweaked ");
        FIELD(result, "The result of the expiry tweak");
    }
    static IObject* defaultIRGridResult(){
        return new IRGridResult();
    }
};

CClassConstSP const IRGridResult::TYPE = CClass::registerClassLoadMethod(
    "IRGridResult", typeid(IRGridResult), IRGridResultHelper::load);

DEFINE_TEMPLATE_TYPE(IRGridResultArray);

DRLIB_END_NAMESPACE
