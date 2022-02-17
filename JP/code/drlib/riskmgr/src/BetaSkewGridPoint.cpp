//----------------------------------------------------------------------------
//
//   Group       : Quantitative Research & Development
//
//   Filename    : BetaSkewGridPoint.cpp
//
//   Description : Identifies a point on a beta skew surface
//
//   Author      : Antoine Gregoire
//
//   Date        : 17-Jun-2005
//
//
//   $Log: $
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/BetaSkewGridPoint.hpp"

DRLIB_BEGIN_NAMESPACE

BetaSkewGridPoint::~BetaSkewGridPoint(){}

BetaSkewGridPoint::BetaSkewGridPoint(): CObject(TYPE) {}

BetaSkewGridPoint::BetaSkewGridPoint(const DateTime& maturity, double strike): 
    CObject(TYPE), maturity(maturity), strike(strike) {}

/** returns the maturity associated with this result */
DateTime BetaSkewGridPoint::getMaturity() const {
    return maturity;
}

/** returns the strike associated with this result */
double BetaSkewGridPoint::getStrike() const {
    return strike;
}

/** converts a BetaSkewGridPointSet into an array */
BetaSkewGridPointArrayConstSP BetaSkewGridPoint::toArray(const BetaSkewGridPointSet& points) {
    BetaSkewGridPointArraySP result(new BetaSkewGridPointArray());
    BetaSkewGridPointSet::const_iterator iter = points.begin();
    for(;iter != points.end(); iter++) {
        result->push_back(BetaSkewGridPointSP(
            new BetaSkewGridPoint((*iter)->getMaturity(), (*iter)->getStrike())));
    }
    return result;
}

/** Invoked when class is 'loaded' */
void BetaSkewGridPoint::load(CClassSP& clazz){
    clazz->setPublic(); // make visible to EAS/spreadsheet
    REGISTER(BetaSkewGridPoint, clazz);
    SUPERCLASS(CObject);
    EMPTY_SHELL_METHOD(defaultBetaSkewGridPoint);
    FIELD(maturity, "maturity");
    FIELD(strike, "strike");
}

IObject* BetaSkewGridPoint::defaultBetaSkewGridPoint(){
    return new BetaSkewGridPoint();
}

CClassConstSP const BetaSkewGridPoint::TYPE =
    CClass::registerClassLoadMethod(
        "BetaSkewGridPoint",
        typeid(BetaSkewGridPoint),
        BetaSkewGridPoint::load);

DEFINE_TEMPLATE_TYPE(BetaSkewGridPointArray);

DRLIB_END_NAMESPACE

