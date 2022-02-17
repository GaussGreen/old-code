//----------------------------------------------------------------------------
//
//   Group       : Credit Hybrids QR
//
//   Filename    : PiecewiseLinearMappingFunction.cpp
//
//   Description : Piecewise Linear Mapping function implementations
//
//   Author      : Sebastien Gay
//
//   Date        : March 2006
//
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/PiecewiseLinearMappingFunction.hpp"
#include "edginc/Format.hpp"

DRLIB_BEGIN_NAMESPACE

/** Linear interpolation between 2 points */
static void linearInterpolation(
    double x, 
    double *y, 
    double x0, 
    double x1, 
    double y0, 
    double y1)
{
    double a,b,l0,l1;

    /* basic check */
    if (x0 ==x1) {
        throw ModelException("linterp",
                "x inputs supplied for interpolation are the same.");
    }

    if ((x == x0) || (y0 == y1)) {
        *y = y0;
    } else if (x == x1) {
        *y = y1;
    } else {
        a=x-x1;
        b=x0-x1;
        l0=a/b;

        a=x-x0;
        b=x1-x0;
        l1=a/b;

        *y=l0*y0+l1*y1;
    }
}

/** Destructor */
PiecewiseLinearMappingFunction::~PiecewiseLinearMappingFunction(){}

/** The mapping function */
double PiecewiseLinearMappingFunction::map(double x) const
{
        // result
        double y;
        
        int nbVal = xPoints->size();
        
        if (nbVal == 1) {
            /* if only one point, we extrapolate flat */
            y = (*yPoints)[0];
        } else if (x <= (*xPoints)[0]) {
                y = (*yPoints)[0];
        } else if (x >= (*xPoints)[nbVal-1]) {
                y = (*yPoints)[nbVal-1];
        } else {
            int k=0;
            while ((x > (*xPoints)[k]) && (k < nbVal - 1)) k++;
            linearInterpolation(x,&y,(*xPoints)[k-1],(*xPoints)[k],(*yPoints)[k-1],(*yPoints)[k]);
        }
        return y;
    }

/** Called immediately after object constructed */
void PiecewiseLinearMappingFunction::validatePop2Object() {
    try{
        PiecewiseMappingFunction::validatePop2Object();
        if (xPoints->size() != yPoints->size()) {
            throw ModelException(
                "xPoints array (size = " +
                Format::toString(xPoints->size()) +
                ") and yPoints array (size = " + 
                Format::toString(yPoints->size()) +
                ") have not the same size");
        }
    } catch (exception& e){
        throw ModelException(e, "PiecewiseLinearMappingFunction::validatePop2Object");
    }
}

/** Invoked when Class is 'loaded' */
void PiecewiseLinearMappingFunction::load(CClassSP& clazz) {
    clazz->setPublic(); // make visible to EAS/spreadsheet
    REGISTER(PiecewiseLinearMappingFunction, clazz);
    SUPERCLASS(PiecewiseMappingFunction);
    EMPTY_SHELL_METHOD(defaultConstructor);
}
    
/** Only build instances of that class using reflection */
PiecewiseLinearMappingFunction::PiecewiseLinearMappingFunction() : PiecewiseMappingFunction(TYPE) {}

/** Constructor from the fields */
PiecewiseLinearMappingFunction::PiecewiseLinearMappingFunction(CDoubleArraySP x, CDoubleArraySP y) : PiecewiseMappingFunction(TYPE) 
{
    xPoints = x;
    yPoints = y;
}

/** Default constructor */
IObject* PiecewiseLinearMappingFunction::defaultConstructor() {
    return new PiecewiseLinearMappingFunction();
}

CClassConstSP const PiecewiseLinearMappingFunction::TYPE =
    CClass::registerClassLoadMethod(
        "PiecewiseLinearMappingFunction",
      typeid(PiecewiseLinearMappingFunction),
        PiecewiseLinearMappingFunction::load);


DRLIB_END_NAMESPACE
