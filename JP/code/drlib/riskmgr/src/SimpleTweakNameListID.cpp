/**
 * @file SimpleTweakNameListID.cpp
 */

#include "edginc/config.hpp"
#include "edginc/MarketObject.hpp"
#include "edginc/SimpleTweakNameListID.hpp"

DRLIB_BEGIN_NAMESPACE

SimpleTweakNameListID::SimpleTweakNameListID(CClassConstSP subjectType):
    CObject(TYPE),
    subjectType(subjectType)
{
    if (!MarketObject::TYPE->isAssignableFrom(subjectType)) {
        throw ModelException(
            "SimpleTweakNameListID::SimpleTweakNameListID",
            "'subjectType' should be a subclass of MarketObject but "
            "is " + (subjectType ? string("null") : subjectType->getName()));
    }
}

SimpleTweakNameListID::~SimpleTweakNameListID() {}

CClassConstSP SimpleTweakNameListID::shiftInterface() const {
    return subjectType;
}

void SimpleTweakNameListID::appendName(OutputNameArray& namesList,
                                       IObjectConstSP obj) {
    try {
        namesList.push_back(OutputNameSP(
            new OutputName(MarketObjectConstSP::dynamicCast(obj)->getName())));
    }
    catch (exception& e) {
        throw ModelException(e, "SimpleTweakNameListID::appendName");
    }
}

void SimpleTweakNameListID::load(CClassSP& clazz) {
    REGISTER(SimpleTweakNameListID, clazz);
    SUPERCLASS(CObject);
    //IMPLEMENTS(ITweakNameListID);
}

CClassConstSP const SimpleTweakNameListID::TYPE = CClass::registerClassLoadMethod(
    "SimpleTweakNameListID", typeid(SimpleTweakNameListID), load);

DRLIB_END_NAMESPACE

// Local Variables: ***
// c-basic-offset:4 ***
// End: ***
