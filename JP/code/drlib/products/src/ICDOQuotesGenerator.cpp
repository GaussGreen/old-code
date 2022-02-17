//----------------------------------------------------------------------------
//
//   Group       : Quantitative Research
//
//   Filename    : ICDOQuotesGenerator.hpp
//
//   Description : ICDOQuotesGenerator is an interface for
//                  generators of pseudo CDO quotes to be used
//                  afterwards in a calibration.
//
//   Author      : Sebastien Gay    
//
//   Date        : October 2006
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/ICDOQuotesGenerator.hpp"
#ifndef CDO_QUOTES_HPP
#include "edginc/CDOQuotes.hpp"
#endif

DRLIB_BEGIN_NAMESPACE

DEFINE_TEMPLATE_TYPE(ICDOQuotesGeneratorWrapper);

/** Destructor */
ICDOQuotesGenerator::~ICDOQuotesGenerator() {}

/** Invoked when Class is 'loaded' */
void ICDOQuotesGenerator::load(CClassSP& clazz) {
    REGISTER_INTERFACE(ICDOQuotesGenerator, clazz);
    EXTENDS(IObject);
}

/** TYPE for CDOQuotesGenerator */
CClassConstSP const ICDOQuotesGenerator::TYPE = CClass::registerInterfaceLoadMethod(
    "ICDOQuotesGenerator", typeid(ICDOQuotesGenerator), ICDOQuotesGenerator::load);
    

/* external symbol to allow class to be forced to be linked in */
bool ICDOQuotesGeneratorLoad(){
    return (ICDOQuotesGenerator::TYPE != 0);
}

DRLIB_END_NAMESPACE
