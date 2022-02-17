//----------------------------------------------------------------------------
//
//   Group       : Credit Hybrids QR
//
//   Filename    : PiecewiseIncrementalMappingFunction.cpp
//
//   Description : Piecewise Incremental Mapping function implementations
//
//   Author      : Sebastien Gay
//
//   Date        : April 2006
//
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/PiecewiseIncrementalMappingFunction.hpp"
#include "edginc/Format.hpp"
#include "edginc/TweakableWith.hpp"
#include "edginc/SqueezeParallelTweak.hpp"
#include "edginc/Calibrator.hpp"

DRLIB_BEGIN_NAMESPACE

/** Destructor */
PiecewiseIncrementalMappingFunction::~PiecewiseIncrementalMappingFunction(){}

/** utility function to compute the xPoints and yPoints of a PiecewiseMappingFunction from
  * the base points, increments and values of a PiecewiseIncrementalMappingFunction
  */
void PiecewiseIncrementalMappingFunction::computeXandY(
                    CDoubleArraySP x,
                    CDoubleArraySP y) const
{    
    int size_x = xPoints->size();
    int size_y = yPoints->size();

    x->resize(size_x + 1);
    y->resize(size_y);

    (*x)[0] = getBasePoint();
    for (int i=1;i<size_x+1;++i)
    {
        // increment the x points
        (*x)[i] = (*x)[i-1] + (*xPoints)[i-1];
    }    

    for (int i=0;i<size_y;++i)
    {
        (*y)[i] = (*yPoints)[i];
    }
}

/** Called immediately after object constructed */
void PiecewiseIncrementalMappingFunction::validatePop2Object() {
    try{
        if (xPoints->size() == 0) {
            throw ModelException(
                "xPoints and yPoints arrays are empty !");
        }
        for (int i=0; i<xPoints->size(); ++i)
        {
            if ((*xPoints)[i]<=0)
            {
                throw ModelException(
                    "Increments in the incremental mapping function are not all > 0");
            }
        }
    } catch (exception& e){
        throw ModelException(e, "PiecewiseIncrementalMappingFunction::validatePop2Object");
    }
}

// Squeeze parallel tweak (sensitivity) support
string PiecewiseIncrementalMappingFunction::sensName(SqueezeParallelTwk* shift) const {
    return name;
}

// utility to get y points
double PiecewiseIncrementalMappingFunction::getBasePoint() const
{
    return basePoint;
}

/** Returns the name of this object. This is the name with which
    it is stored in the market data cache and is the name with
    which results (eg tweaks) should be reported against */
string PiecewiseIncrementalMappingFunction::getName() const{
    return name;
}

/** Only build instances of that class using reflection */
PiecewiseIncrementalMappingFunction::PiecewiseIncrementalMappingFunction(CClassConstSP clazz) : PiecewiseMappingFunction(clazz) {}

// Squeeze parallel tweak (sensitivity) support
bool PiecewiseIncrementalMappingFunction::sensShift(SqueezeParallelTwk* shift) {
    double multiplier = shift->getShiftSize();
    for (int i = 0; i < yPoints->size(); i++) {
        (*yPoints)[i] *= multiplier;			
	}
    return false;
}


void PiecewiseIncrementalMappingFunction::load(CClassSP& clazz) {
    clazz->setPublic(); // make visible to EAS/spreadsheet
    REGISTER(PiecewiseIncrementalMappingFunction, clazz);
    SUPERCLASS(PiecewiseMappingFunction);
    IMPLEMENTS(TweakableWith<SqueezeParallelTwk>);
    IMPLEMENTS(Calibrator::IAdjustable);
    FIELD(basePoint, "initial x point");

    // Register "xPoints" as a field that can be calibrated
    Calibrator::IAdjustable::registerField(
        clazz,
        "basePoint",
        new Range(OpenBoundary(-7.5),  OpenBoundary(7.5)));
}

CClassConstSP const PiecewiseIncrementalMappingFunction::TYPE =
    CClass::registerClassLoadMethod(
        "PiecewiseIncrementalMappingFunction",
      typeid(PiecewiseIncrementalMappingFunction),
        PiecewiseIncrementalMappingFunction::load);

DRLIB_END_NAMESPACE
