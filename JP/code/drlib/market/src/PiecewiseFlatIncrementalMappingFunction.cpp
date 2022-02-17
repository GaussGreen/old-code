//----------------------------------------------------------------------------
//
//   Group       : Credit Hybrids QR
//
//   Filename    : PiecewiseFlatIncrementalMappingFunction.cpp
//
//   Description : Piecewise Flat Incremental Mapping function implementations
//
//   Author      : Sebastien Gay
//
//   Date        : April 2006
//
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/PiecewiseFlatIncrementalMappingFunction.hpp"
#include "edginc/Format.hpp"

DRLIB_BEGIN_NAMESPACE 

/** The mapping function */
double PiecewiseFlatIncrementalMappingFunction::map(double x) const
{
    return mappingFunction->map(x);
}

/** Called immediately after object constructed */
void PiecewiseFlatIncrementalMappingFunction::validatePop2Object() {
    try{
        PiecewiseIncrementalMappingFunction::validatePop2Object();

        //update mapping function
        CDoubleArraySP x = CDoubleArraySP(new CDoubleArray);
        CDoubleArraySP y = CDoubleArraySP(new CDoubleArray);
        this->computeXandY(x, y);
    
        // create the mapping function corresponding to the increments
        mappingFunction = PiecewiseFlatMappingFunctionSP(new PiecewiseFlatMappingFunction(x,y));
        mappingFunction->validatePop2Object();
    } catch (exception& e){
        throw ModelException(e, "PiecewiseFlatIncrementalMappingFunction::validatePop2Object");
    }
}

/* Utility to get the input arrays and fields */
CDoubleArrayConstSP PiecewiseFlatIncrementalMappingFunction::getX() const {
    return mappingFunction->getX();
}

/* Utility to get the input arrays and fields */
CDoubleArrayConstSP PiecewiseFlatIncrementalMappingFunction::getY() const {
    return mappingFunction->getY();
}

/** Called immediately after object modified */
void PiecewiseFlatIncrementalMappingFunction::fieldsUpdated(const CFieldArray& fields)
{
    //update mapping function
    CDoubleArraySP x = CDoubleArraySP(new CDoubleArray);
    CDoubleArraySP y = CDoubleArraySP(new CDoubleArray);
    this->computeXandY(x, y);

    // create the mapping function corresponding to the increments
    mappingFunction = PiecewiseFlatMappingFunctionSP(new PiecewiseFlatMappingFunction(x,y));
}

/** Destructor */
PiecewiseFlatIncrementalMappingFunction::~PiecewiseFlatIncrementalMappingFunction(){}

/** Invoked when Class is 'loaded' */
void PiecewiseFlatIncrementalMappingFunction::load(CClassSP& clazz) {
    clazz->setPublic(); // make visible to EAS/spreadsheet
    REGISTER(PiecewiseFlatIncrementalMappingFunction, clazz);
    SUPERCLASS(PiecewiseIncrementalMappingFunction);
    EMPTY_SHELL_METHOD(defaultConstructor);
    FIELD(mappingFunction,"actual mapping function computed from the incremental mapping function");
    FIELD_MAKE_TRANSIENT(mappingFunction);
}
    
/** Only build instances of that class using reflection */
PiecewiseFlatIncrementalMappingFunction::PiecewiseFlatIncrementalMappingFunction() : PiecewiseIncrementalMappingFunction(TYPE) {}
    
/** Default constructor */
IObject* PiecewiseFlatIncrementalMappingFunction::defaultConstructor() {
    return new PiecewiseFlatIncrementalMappingFunction();
}

CClassConstSP const PiecewiseFlatIncrementalMappingFunction::TYPE =
    CClass::registerClassLoadMethod(
        "PiecewiseFlatIncrementalMappingFunction",
      typeid(PiecewiseFlatIncrementalMappingFunction),
        PiecewiseFlatIncrementalMappingFunction::load);


DRLIB_END_NAMESPACE
