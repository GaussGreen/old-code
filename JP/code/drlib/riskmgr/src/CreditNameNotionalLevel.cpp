//----------------------------------------------------------------------------
//
//   Group       : QR - Credit Hybrids
//
//   Filename    : NotionalLevel.hpp
//
//   Description : Change in name notional scenario 
//
//   Author      : Linus Thand
//
//   Date        : 25 May 2006
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/CreditNameNotionalLevel.hpp"

DRLIB_BEGIN_NAMESPACE

/** returns the interface identifying what an object has to do in order
    to be support the tweak that this object represents */
CClassConstSP CreditNameNotionalLevel::shiftInterface() const 
{ 
    return Shift::TYPE; 
}

CreditNameNotionalLevel::CreditNameNotionalLevel() : ScalarPerturbation(TYPE) {} /** for reflection */

CreditNameNotionalLevel::~CreditNameNotionalLevel() {} 

/** Returns true if the supplied object matches the supplied name
    for this sensitivity. The object must implement this senstivity's
    Shift interface */
bool CreditNameNotionalLevel::nameMatches(const OutputName& name, 
                                IObjectConstSP obj) {
    const Shift& shiftObj = dynamic_cast<const Shift&>(*obj);
    return name.equals(shiftObj.sensName(this));
}

void CreditNameNotionalLevel::appendName(OutputNameArray& namesList,
                               IObjectConstSP   obj){
    const Shift& CreditNameNotionalLevelObj = dynamic_cast<const Shift&>(*obj);
    OutputNameSP outputName(new OutputName(CreditNameNotionalLevelObj.sensName(this)));
    namesList.push_back(outputName);
}

bool CreditNameNotionalLevel::shift(IObjectSP obj) {
    Shift& CreditNameNotionalLevelObj = dynamic_cast<Shift&>(*obj);
    return CreditNameNotionalLevelObj.sensShift(this);
}

void CreditNameNotionalLevel::load(CClassSP& clazz){
   clazz->setPublic(); // make visible to EAS/spreadsheet
   REGISTER(CreditNameNotionalLevel, clazz);
   SUPERCLASS(ScalarPerturbation);
   EMPTY_SHELL_METHOD(defaultConstructor);
}

IObject* CreditNameNotionalLevel::defaultConstructor(){
   return new CreditNameNotionalLevel();
}

CClassConstSP const CreditNameNotionalLevel::TYPE = 
    CClass::registerClassLoadMethod("CreditNameNotionalLevel", 
                                    typeid(CreditNameNotionalLevel), 
                                    load);

CClassConstSP const CreditNameNotionalLevel::Shift::TYPE =
    CClass::registerInterfaceLoadMethod("CreditNameNotionalLevel::Shift", 
                                        typeid(CreditNameNotionalLevel::Shift), 
                                        0);

bool CreditNameNotionalLevelLinkIn ()  {
    return CreditNameNotionalLevel::TYPE != NULL;
}

DRLIB_END_NAMESPACE
