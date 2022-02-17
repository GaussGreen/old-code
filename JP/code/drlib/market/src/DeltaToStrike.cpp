//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : DeltaToStrike.cpp
//
//   Description : Abstract base for IR vols
//
//   Author      : Andrew J Swain
//
//   Date        : 11 January 2002
//
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/DeltaToStrike.hpp"

DRLIB_BEGIN_NAMESPACE

IDeltaToStrikeMaker::IDeltaToStrike::~IDeltaToStrike(){}

void IDeltaToStrikeMaker::load(CClassSP& clazz){
    clazz->setPublic(); // make visible to EAS/spreadsheet
    REGISTER_INTERFACE(IDeltaToStrikeMaker, clazz);
    EXTENDS(IObject);
}

CClassConstSP const IDeltaToStrikeMaker::TYPE = 
CClass::registerInterfaceLoadMethod("IDeltaToStrikeMaker", typeid(IDeltaToStrikeMaker), load);

bool DeltaToStrikeLoad() {
    return IDeltaToStrikeMaker::TYPE != NULL;
}

DRLIB_END_NAMESPACE
