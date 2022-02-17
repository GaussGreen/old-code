//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : Interval.hpp
//
//   Description : A Date interval
//
//   Author      : André Segger
//
//   Date        : 05 October 2001
//
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#define QLIB_INTERVAL_CPP
#include "edginc/Interval.hpp"

DRLIB_BEGIN_NAMESPACE

/** returns the start date of the interval */
const DateTime& Interval::getStartDate() const
{
    return startDate;
}

/** returns the start date of the interval */
const DateTime& Interval::getEndDate() const
{
    return endDate;
}

void Interval::validatePop2Object()
{
   static const string method = "Interval::validatePop2Object";
   if (endDate < startDate) {
       throw ModelException(method, "Start date (" + startDate.toString() +
                             ") is after end date (" + startDate.toString() + ")");
   }
}

/** returns the number of days between the two dates */
int Interval::daysDiff() const
{
    return endDate.getDate() - startDate.getDate();
}

/** public as DateTimeArray is an array of structures */
Interval::Interval(): CObject(TYPE) {
    // empty
}

//// Added to allow array template to be instantiated
//// Implementation just compares two dates.
bool Interval::operator==(const Interval& rhs) const{
    return (startDate == rhs.startDate && endDate == rhs.endDate);
}

IObjectConstSP arrayObjectCast<Interval>::toIObject(
    const Interval& value){
    IObjectConstSP objValue(IObjectConstSP::attachToRef(&value));
    return objValue;
}

/** Casts array element to an IObject */
IObjectSP arrayObjectCast<Interval>::toIObject(Interval& value){
    return IntervalSP::attachToRef(&value);
}

/** Turns the IObjectSP into a DateTime */
const Interval& arrayObjectCast<Interval>::fromIObject(IObjectSP& value){
    Interval* dtPtr = dynamic_cast<Interval*>(value.get());
    if (!dtPtr){
        throw ModelException("arrayObjectCast::fromIObject", "Object is not"
                             " an Interval");
    }
    return *dtPtr;
}

class IntervalHelper{
public:

    /** Invoked when DateTime class is 'loaded' */
    static void load(CClassSP& clazz){
        clazz->setPublic(); // make visible to EAS/spreadsheet
        REGISTER(Interval, clazz);
        SUPERCLASS(CObject);
        EMPTY_SHELL_METHOD(defaultInterval);
        FIELD(startDate, "start date of the interval");
        FIELD(endDate,   "end date of the interval");
    }

    static IObject* defaultInterval(){
        return new Interval();
    }

};

CClassConstSP const Interval::TYPE = CClass::registerClassLoadMethod(
    "Interval", typeid(Interval), IntervalHelper::load);

//template<> CClassConstSP const IntervalArray::TYPE = CClass::registerClassLoadMethod(
//    "IntervalArray", typeid(IntervalArray), load);
DEFINE_TEMPLATE_TYPE(IntervalArray);

DRLIB_END_NAMESPACE
