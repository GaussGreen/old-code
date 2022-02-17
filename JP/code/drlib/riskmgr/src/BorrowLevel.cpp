//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : BorrowLevel.cpp
//
//   Description : borrow level scenario - overwrites the borrow curve with a flat curve
//
//   Author      : André Segger
//
//   Date        : 01 July 2002
//
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/BorrowLevel.hpp"

DRLIB_BEGIN_NAMESPACE

BorrowLevel::Shift::~Shift(){} // empty

/** constructor with explicit shift */
BorrowLevel::BorrowLevel(double shift):
    ScalarPerturbation(TYPE, shift){}

/** for reflection */
BorrowLevel::BorrowLevel(): ScalarPerturbation(TYPE){}

/** returns the interface identifying what an object has to do in order
    to be support the tweak that this object represents */
CClassConstSP BorrowLevel::shiftInterface() const{
    return Shift::TYPE;
}
 
/** Returns true if the supplied object matches the supplied name
    for this sensitivity. The object must implement this senstivity's
    Shift interface */
bool BorrowLevel::nameMatches(const OutputName&         name,
                              IObjectConstSP          obj){
    // cast obj to BorrowLevel::Shift and then invoke name method
    const Shift& shiftObj = 
        dynamic_cast<const Shift&>(*obj);
    return name.equals(shiftObj.sensName(this));
}

/** Appends the name(s) of the supplied object with respect to this
    sensitivity to the supplied list. The object must implement this
    senstivity's Shift interface */
void BorrowLevel::appendName(OutputNameArray&          namesList,
                             IObjectConstSP          obj){
    // cast obj to BorrowLevel::Shift and then invoke name method
    const Shift& shiftObj = 
        dynamic_cast<const Shift&>(*obj);
    OutputNameSP outputName(new OutputName(shiftObj.sensName(this)));
    namesList.push_back(outputName);
}

bool BorrowLevel::shift(IObjectSP obj) {
    // cast obj to BorrowLevel::Shift and then invoke shift method
    Shift& shiftObj = 
        dynamic_cast<Shift&>(*obj);
    return shiftObj.sensShift(this);
}

class BorrowLevelHelper{
public:
    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        clazz->setPublic(); // make visible to EAS/spreadsheet
        REGISTER(BorrowLevel, clazz);
        SUPERCLASS(ScalarPerturbation);
        EMPTY_SHELL_METHOD(defaultBorrowLevel);
    }

    static IObject* defaultBorrowLevel(){
        return new BorrowLevel();
    }
};

CClassConstSP const BorrowLevel::TYPE = CClass::registerClassLoadMethod(
    "BorrowLevel", typeid(BorrowLevel), BorrowLevelHelper::load);

CClassConstSP const BorrowLevel::Shift::TYPE =
CClass::registerInterfaceLoadMethod(
    "BorrowLevel::Shift", typeid(BorrowLevel::Shift), 0);

DRLIB_END_NAMESPACE
