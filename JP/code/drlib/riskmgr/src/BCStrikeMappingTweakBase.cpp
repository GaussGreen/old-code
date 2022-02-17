//----------------------------------------------------------------------------
//
//   Group       : QR
//
//   Filename    : BCStrikeMappingTweakBase.cpp
//
//   Description : Tweak to shift/set the Base Correlation Strike Mapping
//
//   Author      : Jose Hilera
//
//   Date        : 11 October 2005
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/BCStrikeMappingTweakBase.hpp"

DRLIB_BEGIN_NAMESPACE


BCStrikeMappingTweakBase::IShift::~IShift()
{} // empty

/** constructor with no shift size */
BCStrikeMappingTweakBase::BCStrikeMappingTweakBase(CClassConstSP clazz) :
    ScalarPerturbation(clazz)
{}


/** returns the interface identifying what an object has to do in order
    to support the tweak that this object represents */
CClassConstSP BCStrikeMappingTweakBase::shiftInterface() const {
    return IShift::TYPE;
}


/** Returns true if the supplied object matches the supplied name
    for this sensitivity. The object must implement this sensitivity's
    Shift interface */
bool BCStrikeMappingTweakBase::nameMatches(const OutputName&  name,
                                           IObjectConstSP     obj)
{
    // cast obj to BCStrikeMappingTweakBase::Shift and then invoke name method
    const IShift& bcStrikeMappingTweakObj = dynamic_cast<const IShift&>(*obj);
    return name.equals(bcStrikeMappingTweakObj.sensName(this));
}


/** Appends the name(s) of the supplied object with respect to this
    sensitivity to the supplied list. The object must implement this
    sensitivity's Shift interface */
void BCStrikeMappingTweakBase::appendName(OutputNameArray&   namesList,
                                          IObjectConstSP     obj)
{
    // cast obj to BCStrikeMappingTweakBase::Shift and then invoke name method
    const IShift& bcStrikeMappingTweakObj = dynamic_cast<const IShift&>(*obj);
    OutputNameSP outputName(new OutputName(bcStrikeMappingTweakObj.sensName(this)));
    namesList.push_back(outputName);
}

bool BCStrikeMappingTweakBase::shift(IObjectSP obj) {
    // cast obj to BCStrikeMappingTweakBase::Shift and then invoke shift method
    IShift& bcStrikeMappingTweakObj = dynamic_cast<IShift&>(*obj);
    return bcStrikeMappingTweakObj.sensShift(this);
}


/** Invoked when class is 'loaded' */
void BCStrikeMappingTweakBase::load(CClassSP& clazz) {
    clazz->setPublic(); // make visible to EAS/spreadsheet
    REGISTER(BCStrikeMappingTweakBase, clazz);
    SUPERCLASS(ScalarPerturbation);
};


CClassConstSP const BCStrikeMappingTweakBase::TYPE = 
    CClass::registerClassLoadMethod("BCStrikeMappingTweakBase", 
                                    typeid(BCStrikeMappingTweakBase), 
                                    load);

CClassConstSP const BCStrikeMappingTweakBase::IShift::TYPE =
    CClass::registerInterfaceLoadMethod("BCStrikeMappingTweakBase::Shift", 
                                        typeid(BCStrikeMappingTweakBase::IShift), 
                                        0);

/**
 * Included in RiskMgrLib::linkInClasses() to force linkage into the Windows exe
 */
bool BCStrikeMappingTweakBaseLinkIn() {
    return BCStrikeMappingTweakBase::TYPE != NULL;
}

DRLIB_END_NAMESPACE
