//----------------------------------------------------------------------------
//
//   Group       : QR - Credit Hybrids
//
//   Filename    : CreditMultiBetaLevel.cpp
//
//   Description : Beta substitution scenario -
//                 Supply a list of names and their corresponding new betas.
//
//   Author      : Linus Thand
//
//   Date        : 24 May 2006
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/CreditMultiBetaLevel.hpp"

DRLIB_BEGIN_NAMESPACE

CreditMultiBetaLevel::~CreditMultiBetaLevel() {}

CreditMultiBetaLevel::CreditMultiBetaLevel() : Perturbation (TYPE)  {}     // for reflection

/** returns the interface identifying what an object has to do in order
    to be support the tweak that this object represents */
CClassConstSP CreditMultiBetaLevel::shiftInterface() const{
    return IShift::TYPE;
}
   
StringArrayConstSP CreditMultiBetaLevel::getNames() const {
    return StringArrayConstSP(names);
}

DoubleArrayConstSP CreditMultiBetaLevel::getBetas() const {
    return DoubleArrayConstSP(betas);
}
    
/** Returns true if the supplied object matches the supplied name
    for this sensitivity. The object must implement this senstivity's
    Shift interface */
bool CreditMultiBetaLevel::nameMatches(const OutputName& name,
                                 IObjectConstSP    obj){
    // cast obj to CreditMultiBetaLevel::Shift and then invoke name method
    const IShift& betaShiftObj = dynamic_cast<const IShift&>(*obj);
    return name.equals(betaShiftObj.sensName(this));
}

/** Appends the name(s) of the supplied object with respect to this
    sensitivity to the supplied list. The object must implement this
    senstivity's Shift interface */
void CreditMultiBetaLevel::appendName(OutputNameArray& namesList,
                                IObjectConstSP   obj){
    // cast obj to VolAbsoluteShift::Shift and then invoke name method
    const IShift& betaShiftObj = dynamic_cast<const IShift&>(*obj);
    OutputNameSP outputName(new OutputName(betaShiftObj.sensName(this)));
    namesList.push_back(outputName);
}

bool CreditMultiBetaLevel::shift(IObjectSP obj) {
    // cast obj to CreditMultiBetaLevel::Shift and then invoke shift method
    IShift& betaShiftObj = dynamic_cast<IShift&>(*obj);
    return betaShiftObj.sensShift(this);
    return true;
}

void CreditMultiBetaLevel::validatePop2Object() {  
    static const string method("CreditMultiBetaLevel::validatePop2Object");
    try {
        if (names->empty()) {
            throw ModelException(method, "list of names is empty");
        }
        if (betas->empty()) {
            throw ModelException(method, "list of betas is empty");
        }
        if (names->size() != betas->size()) {
            throw ModelException(method, "names & betas have different lengths");
        }
    }
    catch (exception& e) {
        throw ModelException(e, method);
    }
}

/** Invoked when class is 'loaded' */
void CreditMultiBetaLevel::load(CClassSP& clazz){
    clazz->setPublic(); // make visible to EAS/spreadsheet
    REGISTER(CreditMultiBetaLevel, clazz);
    SUPERCLASS(Perturbation);
    EMPTY_SHELL_METHOD(defaultConstructor);
    FIELD(names, "the portfolio names to change");
    FIELD(betas, "the betas to be replaced");
}

IObject* CreditMultiBetaLevel::defaultConstructor(){
    return new CreditMultiBetaLevel();
}

CClassConstSP const CreditMultiBetaLevel::TYPE = CClass::registerClassLoadMethod(
    "CreditMultiBetaLevel", typeid(CreditMultiBetaLevel), load);


CClassConstSP const CreditMultiBetaLevel::IShift::TYPE = CClass::registerInterfaceLoadMethod(
    "CreditMultiBetaLevel::IShift", typeid(CreditMultiBetaLevel::IShift), 0);


bool CreditMultiBetaLevelLinkIn () {
    return CreditMultiBetaLevel::TYPE != NULL;
}

DRLIB_END_NAMESPACE

