//----------------------------------------------------------------------------
//
//   Group       : Quantitative Research
//
//   Filename    : IAuxSpreadProcess.hpp
//
//   Description : Interface for auxiliary spread porcess needed by spreadLossTree
//
//   Author      : Matthias Arnsdorf
//
//   Date        : June 2006
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/IAuxSpreadProcess.hpp"

DRLIB_BEGIN_NAMESPACE


IAuxSpreadProcess::~IAuxSpreadProcess()
{}


void IAuxSpreadProcess::load (CClassSP& clazz) {
	REGISTER_INTERFACE(IAuxSpreadProcess, clazz);
	EXTENDS(IObject);
}


CClassConstSP const IAuxSpreadProcess::TYPE = 
CClass::registerInterfaceLoadMethod("IAuxSpreadProcess", 
									typeid(IAuxSpreadProcess), 
									load);

DRLIB_END_NAMESPACE
