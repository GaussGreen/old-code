//----------------------------------------------------------------------------
//
//   Group       : Credit Hybrids QR
//
//   Filename    : PiecewiseFlatMappingFunction.cpp
//
//   Description : Piecewise Flat Mapping function implementations
//
//   Author      : Sebastien Gay
//
//   Date        : March 2006
//
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/PiecewiseFlatMappingFunction.hpp"
#include "edginc/Format.hpp"
//#include "edginc/LessEqualGreaterEps.hpp" // in utils

DRLIB_BEGIN_NAMESPACE 

/** The mapping function */
double PiecewiseFlatMappingFunction::map(double x) const
{
        // result
        double y;
        
        int nbVal = xPoints->size();
        
        if (x <= (*xPoints)[0]) {
                y = (*yPoints)[0];
        } else if (x >= (*xPoints)[nbVal-1]) {
                y = (*yPoints)[nbVal-1];
        } else {
            int k=0;
            while ((x > (*xPoints)[k]) && (k < nbVal - 1)) k++;
            return (*yPoints)[k-1];
        }
        return y;
    }

/** Called immediately after object constructed */
void PiecewiseFlatMappingFunction::validatePop2Object() {
    try{
        PiecewiseMappingFunction::validatePop2Object();
        if (yPoints->size() - xPoints->size() > 1 || yPoints->size() - xPoints->size() < 0) {
            throw ModelException(
                "xPoints array (size = " +
                Format::toString(xPoints->size()) +
                ") and yPoints array (size = " + 
                Format::toString(yPoints->size()) +
                ") are not compatible");
        }
    } catch (exception& e){
        throw ModelException(e, "PiecewiseFlatMappingFunction::validatePop2Object");
    }
}

/** Destructor */
PiecewiseFlatMappingFunction::~PiecewiseFlatMappingFunction(){}

/** Invoked when Class is 'loaded' */
void PiecewiseFlatMappingFunction::load(CClassSP& clazz) {
    clazz->setPublic(); // make visible to EAS/spreadsheet
    REGISTER(PiecewiseFlatMappingFunction, clazz);
    SUPERCLASS(PiecewiseMappingFunction);
    EMPTY_SHELL_METHOD(defaultConstructor);
}
    
/** Only build instances of that class using reflection */
PiecewiseFlatMappingFunction::PiecewiseFlatMappingFunction() : PiecewiseMappingFunction(TYPE) {}
    
/** Constructor from the fields */
PiecewiseFlatMappingFunction::PiecewiseFlatMappingFunction(CDoubleArraySP x, CDoubleArraySP y): PiecewiseMappingFunction(TYPE) 
{
    xPoints = x;
    yPoints = y;
}

/** Default constructor */
IObject* PiecewiseFlatMappingFunction::defaultConstructor() {
    return new PiecewiseFlatMappingFunction();
}

CClassConstSP const PiecewiseFlatMappingFunction::TYPE =
    CClass::registerClassLoadMethod(
        "PiecewiseFlatMappingFunction",
      typeid(PiecewiseFlatMappingFunction),
        PiecewiseFlatMappingFunction::load);


DRLIB_END_NAMESPACE
