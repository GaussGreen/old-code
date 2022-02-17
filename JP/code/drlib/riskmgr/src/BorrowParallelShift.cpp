//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : BorrowParallelShift.cpp
//
//   Description : borrow curve shift scenario - adds a parallel shift
//
//   Author      : Andrew J Swain
//
//   Date        : 10 April 2002
//
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/BorrowParallelShift.hpp"

DRLIB_BEGIN_NAMESPACE

BorrowParallelShift::Shift::~Shift(){} // empty

/** constructor with explicit shift */
BorrowParallelShift::BorrowParallelShift(double shift):
    ScalarPerturbation(TYPE, shift){}

/** for reflection */
BorrowParallelShift::BorrowParallelShift(): ScalarPerturbation(TYPE){}

/** returns the interface identifying what an object has to do in order
    to be support the tweak that this object represents */
CClassConstSP BorrowParallelShift::shiftInterface() const{
    return Shift::TYPE;
}
 
/** Returns true if the supplied object matches the supplied name
    for this sensitivity. The object must implement this senstivity's
    Shift interface */
bool BorrowParallelShift::nameMatches(const OutputName&         name,
                                      IObjectConstSP          obj){
    // cast obj to BorrowParallelShift::Shift and then invoke name method
    const Shift& shiftObj = 
        dynamic_cast<const Shift&>(*obj);
    return name.equals(shiftObj.sensName(this));
}

/** Appends the name(s) of the supplied object with respect to this
    sensitivity to the supplied list. The object must implement this
    senstivity's Shift interface */
void BorrowParallelShift::appendName(OutputNameArray&          namesList,
                                     IObjectConstSP          obj){
    // cast obj to BorrowParallelShift::Shift and then invoke name method
    const Shift& shiftObj = 
        dynamic_cast<const Shift&>(*obj);
    OutputNameSP outputName(new OutputName(shiftObj.sensName(this)));
    namesList.push_back(outputName);
}

bool BorrowParallelShift::shift(IObjectSP obj) {
    // cast obj to BorrowParallelShift::Shift and then invoke shift method
    Shift& shiftObj = 
        dynamic_cast<Shift&>(*obj);
    return shiftObj.sensShift(this);
}

class BorrowParallelShiftHelper{
public:
    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        clazz->setPublic(); // make visible to EAS/spreadsheet
        REGISTER(BorrowParallelShift, clazz);
        SUPERCLASS(ScalarPerturbation);
        EMPTY_SHELL_METHOD(defaultBorrowParallelShift);
    }

    static IObject* defaultBorrowParallelShift(){
        return new BorrowParallelShift();
    }
};

CClassConstSP const BorrowParallelShift::TYPE = CClass::registerClassLoadMethod(
    "BorrowParallelShift", typeid(BorrowParallelShift), BorrowParallelShiftHelper::load);

CClassConstSP const BorrowParallelShift::Shift::TYPE =
CClass::registerInterfaceLoadMethod(
    "BorrowParallelShift::Shift", typeid(BorrowParallelShift::Shift), 0);

DRLIB_END_NAMESPACE
