//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : MappingFunction.cpp
//
//   Description : Mapping function implementations
//
//   Author      : Antoine Gregoire
//
//   Date        : March 2005
//
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/MappingFunction.hpp"
#include "edginc/Format.hpp"
#include "edginc/TweakableWith.hpp"
#include "edginc/SqueezeParallelTweak.hpp"

DRLIB_BEGIN_NAMESPACE

/** Loads IMappingFunction interface */
static void loadIMappingFunction(CClassSP& clazz){
    REGISTER_INTERFACE(IMappingFunction, clazz);
    EXTENDS(IObject);
}

/** TYPE for IMappingFunction */
CClassConstSP const IMappingFunction::TYPE =
    CClass::registerInterfaceLoadMethod(
        "IMappingFunction", typeid(IMappingFunction), loadIMappingFunction);

/** Only build instances of that class using reflection */
//IMappingFunction::IMappingFunction() : CObject(TYPE) {}

DRLIB_END_NAMESPACE
