//----------------------------------------------------------------------------
//
//   Group       : Quantitative Research & Development
//
//   Filename    : BetaSkewGridPoint.hpp
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

#ifndef QLIB_BETA_SKEW_GRID_POINT_HPP
#define QLIB_BETA_SKEW_GRID_POINT_HPP

#include "edginc/DateTime.hpp"
#include <set>

DRLIB_BEGIN_NAMESPACE

// predeclarations
class BetaSkewGridPoint;
typedef smartConstPtr<BetaSkewGridPoint> BetaSkewGridPointConstSP;
typedef smartPtr<BetaSkewGridPoint> BetaSkewGridPointSP;
typedef array<BetaSkewGridPointSP, BetaSkewGridPoint> BetaSkewGridPointArray;
typedef smartPtr<BetaSkewGridPointArray> BetaSkewGridPointArraySP;
typedef smartConstPtr<BetaSkewGridPointArray> BetaSkewGridPointArrayConstSP;

struct BetaSkewGridPointCompare;
typedef set<BetaSkewGridPointConstSP, BetaSkewGridPointCompare> BetaSkewGridPointSet;

/** Identifies a point on a beta skew surface */
class RISKMGR_DLL BetaSkewGridPoint : public CObject {
public:
    static CClassConstSP const TYPE;

    virtual ~BetaSkewGridPoint();

    BetaSkewGridPoint(const DateTime& maturity, double strike);

    /** returns the maturity associated with this result */
    DateTime getMaturity() const;

    /** returns the strike associated with this result */
    double getStrike() const;
    
    /** converts a BetaSkewGridPointSet into an array */
    static BetaSkewGridPointArrayConstSP toArray(const BetaSkewGridPointSet& points);

private:

    // Reflection
    BetaSkewGridPoint();
    static IObject* defaultBetaSkewGridPoint();
    static void load(CClassSP& clazz);

    // Fields
    DateTime maturity;
    double strike;
};

/**
 * Compares 2 BetaSkewGridPoint, useful to build sets !
 * We arbitrarily choose to sort BetaSkewGridPoint by increasing maturity and then
 * increasing strike.
 * */
struct BetaSkewGridPointCompare {
    bool operator()(const BetaSkewGridPointConstSP p1, const BetaSkewGridPointConstSP p2) const {
        if (p1->getMaturity() == p2->getMaturity()) {
            return (p1->getStrike() < p2->getStrike());
        } else {
            return (p1->getMaturity() < p2->getMaturity());
        }
    }
};

DRLIB_END_NAMESPACE

#endif //QLIB_BETA_SKEW_GRID_POINT_HPP
