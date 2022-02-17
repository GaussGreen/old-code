//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : YCWeighted.cpp
//
//   Description : Scenario shift: Yield Curve benchmarks can be shifted by
//                 different amounts. If the user defined shiftSizes lie between
//                 Yield Curve benchmarks then we interpolate. Extrapolate flat off
//                 the ends.
//
//   Author      : Jose Hilera
//
//   Date        : 15-June-2005
//
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/YCWeightedShift.hpp"
#include "edginc/CashFlow.hpp"

DRLIB_BEGIN_NAMESPACE

YCWeightedShift::IShift::~IShift(){} // empty

YCWeightedShift::YCWeightedShift(const CClassConstSP& clazz):
    MultiExpiryShift(clazz){}


/** returns the interface identifying what an object has to do in order
    to be support the tweak that this object represents */
CClassConstSP YCWeightedShift::shiftInterface() const{
    return IShift::TYPE;
}
 
/** Returns true if the supplied object matches the supplied name
    for this sensitivity. The object must implement this senstivity's
    Shift interface */
bool YCWeightedShift::nameMatches(const OutputName&  name,
                                  IObjectConstSP     obj){
    // cast obj to YCWeightedShift::Shift and then invoke name method
    return name.equals(dynamic_cast<const IShift&>(*obj).sensName(this));
}

/** Appends the name(s) of the supplied object with respect to this
    sensitivity to the supplied list. The object must implement this
    senstivity's Shift interface */
void YCWeightedShift::appendName(OutputNameArray&       namesList,
                                 IObjectConstSP         obj){
    // cast obj to YCWeightedShift::Shift and then invoke name method
    const IShift& volShiftObj = dynamic_cast<const IShift&>(*obj);
    OutputNameSP outputName(new OutputName(volShiftObj.sensName(this)));
    namesList.push_back(outputName);
}

bool YCWeightedShift::shift(IObjectSP obj) {
    // cast obj to YCWeightedShift::Shift and then invoke shift method
    return dynamic_cast<IShift&>(*obj).sensShift(this);
}

void YCWeightedShift::validatePop2Object() {
    try {
        const static string method = "YCWeightedShift::validatePop2Object";
        MultiExpiryShift::validatePop2Object();

        if (shifts.size() <= 1) {
            throw ModelException(method, 
                                 "shifts & expiries arrays must have more than one element");
        }
    }
    catch (exception& e) {
        throw ModelException(e, "YCWeightedShift::validatePop2Object");
    }
}

/**
 * Shifts the DoubleArray passed as 3rd argument in the dates
 * indicated by the 2nd argument. The first argument determines
 * whether the shift is additive or multiplicative.
 * Note this method is not virtual: subclasses are expected to call
 * back to it with the right arguments, rather than overriding it */
void YCWeightedShift::shiftArray(bool isAdditiveShift,
                                 ExpiryArraySP rateExpiries, 
                                 DoubleArray* rates,
                                 DateTime& today)
{
    double shift;

    for (int i=0; i<rateExpiries->getLength(); i++) {
        shift = shiftSize(today, (*rateExpiries)[i]->toDate(today));
        if (isAdditiveShift) { // do additive shift
            (*rates)[i] += shift;
        }
        else { // do multiplicative shift
            (*rates)[i] *= (1.0 + shift);
        }
        (*rates)[i] = Maths::max((*rates)[i], 0.0); // floor to zero
    }
}


/* What's the shift for a given date */
double YCWeightedShift::shiftSize(const DateTime& today,
                                  const DateTime& shiftDate) const
{
    // Simple interpolation - extrapolate flat
    double shift = CashFlow::interpolate(expiries,
                                         shifts,
                                         today,
                                         shiftDate,
                                         true);
    return shift;
}


class YCWeightedShiftHelper{
public:
    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        clazz->setPublic(); // make visible to EAS/spreadsheet
        clazz->setDescription("Absolute shifts to yield curve by benchmark");
        REGISTER(YCWeightedShift, clazz);
        SUPERCLASS(MultiExpiryShift);
    }
};

CClassConstSP const YCWeightedShift::TYPE = CClass::registerClassLoadMethod(
    "YCWeightedShift", typeid(YCWeightedShift), YCWeightedShiftHelper::load);


CClassConstSP const YCWeightedShift::IShift::TYPE =
    CClass::registerInterfaceLoadMethod(
          "YCWeightedShift::IShift", typeid(YCWeightedShift::IShift), 0);

DRLIB_END_NAMESPACE
