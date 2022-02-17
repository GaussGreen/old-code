//----------------------------------------------------------------------------
//
//   Group       : Credit Hybrids QR
//
//   Filename    : PiecewiseLinearIncrementalMappingFunction.cpp
//
//   Description : Piecewise Linear Incremental Mapping function implementations
//
//   Author      : Sebastien Gay
//
//   Date        : April 2006
//
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/PiecewiseLinearIncrementalMappingFunction.hpp"
#include "edginc/Format.hpp"

DRLIB_BEGIN_NAMESPACE 

/** The mapping function */
double PiecewiseLinearIncrementalMappingFunction::map(double x) const
{
    return mappingFunction->map(x);
}

/** Called immediately after object constructed */
void PiecewiseLinearIncrementalMappingFunction::validatePop2Object() {
    try{
        PiecewiseIncrementalMappingFunction::validatePop2Object();

        CDoubleArraySP x = CDoubleArraySP(new CDoubleArray);
        CDoubleArraySP y = CDoubleArraySP(new CDoubleArray);
        this->computeXandY(x, y);
    
        // create the mapping function corresponding to the increments
        mappingFunction = PiecewiseLinearMappingFunctionSP(new PiecewiseLinearMappingFunction(x,y));
        mappingFunction->validatePop2Object();
    } catch (exception& e){
        throw ModelException(e, "PiecewiseLinearIncrementalMappingFunction::validatePop2Object");
    }
}

/* Utility to get the input arrays and fields */
CDoubleArrayConstSP PiecewiseLinearIncrementalMappingFunction::getX() const {
    return mappingFunction->getX();
}

/* Utility to get the input arrays and fields */
CDoubleArrayConstSP PiecewiseLinearIncrementalMappingFunction::getY() const {
    return mappingFunction->getY();
}

/** Called immediately after object modified */
void PiecewiseLinearIncrementalMappingFunction::fieldsUpdated(const CFieldArray& fields)
{
    //update mapping function
    CDoubleArraySP x = CDoubleArraySP(new CDoubleArray);
    CDoubleArraySP y = CDoubleArraySP(new CDoubleArray);
    this->computeXandY(x, y);

    // create the mapping function corresponding to the increments
    mappingFunction = PiecewiseLinearMappingFunctionSP(new PiecewiseLinearMappingFunction(x,y));
}


/** Destructor */
PiecewiseLinearIncrementalMappingFunction::~PiecewiseLinearIncrementalMappingFunction(){}

/** Invoked when Class is 'loaded' */
void PiecewiseLinearIncrementalMappingFunction::load(CClassSP& clazz) {
    clazz->setPublic(); // make visible to EAS/spreadsheet
    REGISTER(PiecewiseLinearIncrementalMappingFunction, clazz);
    SUPERCLASS(PiecewiseIncrementalMappingFunction);
    EMPTY_SHELL_METHOD(defaultConstructor);
    FIELD(mappingFunction,"actual mapping function computed from the incremental mapping function");
    FIELD_MAKE_TRANSIENT(mappingFunction);
}
    
/** Only build instances of that class using reflection */
PiecewiseLinearIncrementalMappingFunction::PiecewiseLinearIncrementalMappingFunction() : PiecewiseIncrementalMappingFunction(TYPE) {}
    
/** Default constructor */
IObject* PiecewiseLinearIncrementalMappingFunction::defaultConstructor() {
    return new PiecewiseLinearIncrementalMappingFunction();
}

CClassConstSP const PiecewiseLinearIncrementalMappingFunction::TYPE =
    CClass::registerClassLoadMethod(
        "PiecewiseLinearIncrementalMappingFunction",
      typeid(PiecewiseLinearIncrementalMappingFunction),
        PiecewiseLinearIncrementalMappingFunction::load);

DRLIB_END_NAMESPACE
