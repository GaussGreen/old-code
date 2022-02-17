//----------------------------------------------------------------------------
//
//   Group       : Credit Hybrids QR
//
//   Filename    : PiecewiseMappingFunction.cpp
//
//   Description : Piecewise Mapping function implementations
//
//   Author      : Sebastien Gay
//
//   Date        : March 2006
//
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/PiecewiseMappingFunction.hpp"
#include "edginc/Format.hpp"
#include "edginc/TweakableWith.hpp"
#include "edginc/SqueezeParallelTweak.hpp"
#include "edginc/Calibrator.hpp"

DRLIB_BEGIN_NAMESPACE

/** Destructor */
PiecewiseMappingFunction::~PiecewiseMappingFunction(){}

/** Set a name for this function (useful for tweaking) */
void PiecewiseMappingFunction::setName(string name) {
    this->name = name;
}

/** Called immediately after object constructed */
void PiecewiseMappingFunction::validatePop2Object() {
    try{
        if (xPoints->size() == 0) {
            throw ModelException(
                "xPoints and yPoints arrays are empty !");
        }
        for (int i=0; i<xPoints->size()-1; ++i)
        {
            if((*xPoints)[i] >= (*xPoints)[i+1])
            {
                throw ModelException(
                    "xPoints should be increasing ! xPoint[" + Format::toString(i) + "] is greater or equal to xPoint[" + Format::toString(i+1) + "]");                
            }
        }
    } catch (exception& e){
        throw ModelException(e, "PiecewiseMappingFunction::validatePop2Object");
    }
}

// Squeeze parallel tweak (sensitivity) support
string PiecewiseMappingFunction::sensName(SqueezeParallelTwk* shift) const {
    return name;
}

// utility to get x points
CDoubleArrayConstSP PiecewiseMappingFunction::getX() const
{
    return xPoints;
}

// utility to get y points
CDoubleArrayConstSP PiecewiseMappingFunction::getY() const
{
    return yPoints;
}

/** Returns the name of this object. This is the name with which
    it is stored in the market data cache and is the name with
    which results (eg tweaks) should be reported against */
string PiecewiseMappingFunction::getName() const{
    return name;
}

/** Only build instances of that class using reflection */
PiecewiseMappingFunction::PiecewiseMappingFunction(CClassConstSP clazz) : MarketObject(clazz) {}

// Squeeze parallel tweak (sensitivity) support
bool PiecewiseMappingFunction::sensShift(SqueezeParallelTwk* shift) {
    double multiplier = shift->getShiftSize();
    for (int i = 0; i < yPoints->size(); i++) {
        (*yPoints)[i] *= multiplier;			
	}
    return false;
}

// Support to flatten the yPoints
bool PiecewiseMappingFunction::sensShift(QuasiContractualBaseCorrelation* shift) {
    double yPointsLevel = shift->getPiecewiseMappingFunctionLevel();
    for (int i=0; i < yPoints->size(); ++i) {
        (*yPoints)[i] = yPointsLevel;
	}
    return false;
}


void PiecewiseMappingFunction::load(CClassSP& clazz) {
    clazz->setPublic(); // make visible to EAS/spreadsheet
    REGISTER(PiecewiseMappingFunction, clazz);
    SUPERCLASS(MarketObject);
    IMPLEMENTS(IMappingFunction);
    IMPLEMENTS(TweakableWith<SqueezeParallelTwk>);
    IMPLEMENTS(Calibrator::IAdjustable);
    IMPLEMENTS(QuasiContractualBaseCorrelation::IShift);
    FIELD(xPoints, "Array of 'x' points");
    FIELD(yPoints, "Array of 'y'=map(x) points");
    FIELD(name, "name");
    FIELD_MAKE_OPTIONAL(name);

    // Register "skew" as a field that can be calibrated
    Calibrator::IAdjustable::registerField(
        clazz,
        "xPoints",
        new Range(OpenBoundary(-7.5),  OpenBoundary(7.5)));

    // Register "skew" as a field that can be calibrated
    Calibrator::IAdjustable::registerField(
        clazz,
        "yPoints",
        new Range(OpenBoundary(-0.5),  OpenBoundary(1.0)));
}

CClassConstSP const PiecewiseMappingFunction::TYPE =
    CClass::registerClassLoadMethod(
        "PiecewiseMappingFunction",
      typeid(PiecewiseMappingFunction),
        PiecewiseMappingFunction::load);

DRLIB_END_NAMESPACE
