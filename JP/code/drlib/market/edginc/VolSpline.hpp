//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : VolSpline.hpp
//
//   Description : Vol Spline
//
//   Author      : Regis Guichard
//
//   Date        : 04 April 02
//
//
//----------------------------------------------------------------------------

#ifndef VOLSPLINE_HPP
#define VOLSPLINE_HPP

#include "edginc/config.hpp"
#include "edginc/VolBaseParamSurface.hpp"
#include "edginc/VolatilityBS.hpp"
#include "edginc/VolatilityDVF.hpp"
#include "edginc/VolProcessedBS.hpp"
#include "edginc/VolSurface.hpp"
#include "edginc/Asset.hpp"
#include "edginc/Format.hpp"
#include "edginc/Maths.hpp"
#include "edginc/mathlib.hpp"
#include "edginc/VolParam.hpp"
#include "edginc/VolParamTweak.hpp"
#include "edginc/Spline.hpp"
#include "edginc/DeltaSurface.hpp"

DRLIB_BEGIN_NAMESPACE

/** A parameterised CVolBase that supports BS and DVF views of the world */
class MARKET_DLL VolSpline: public CVolBaseParamSurface,
                 virtual public IVolatilityBS,
                 virtual public IVolatilityDVF,
                 virtual public DeltaSurface::IShift {
public:
    friend class VSVolParam;
    friend class VolSplineHelper;
    static CClassConstSP const TYPE;

    void validatePop2Object();

    /** method that builds a CVolParam. */
    virtual CVolParam* createVolParam() const;
    
    VolSpline(const VolSurface& volSurface,
              double strikeRef = 0.0,
              double upperSigmaUnscaled = default_upperSigmaUnscaled);

    virtual string sensName(DeltaSurface* shift) const;
    virtual bool sensShift(DeltaSurface* shift);

protected:
    /*** Build the parameterised vol and cache any values **/
    void buildCache();
    void buildCache(const VolSurface* backbone);

private:
    void ComputeImpVol(const CLatticeDouble&      strikes,
                       const DateTimeArray&       maturities,
                       CLatticeDouble&            impV) const;

    CDoubleMatrixSP spotVolMatrixFromStrikes(
            const double*         strikes,
            int                   strikes_size,
            int                   iBackboneDateStart,
            int                   iBackboneDateEnd) const;

    VolSurface* spotVolSurfaceFromStrikes(const CDoubleArray& strikes) const;

    
    double scalessStrike(double strike, 
                         double strikeRefVols,
                         double sqrTimeFracs) const;
    double rightExtrapolate(double strike, 
                            double upperStrike,
                            double upperVol,
                            double upperSlope,
                            double upperCurvature,
                            double strikeRefVol,
                            double sqrTimeFrac) const;
    double leftExtrapolate(double strike, int iMat) const;

    // do some caching to avoid memory reallocation
    CDoubleMatrixSP createVolMatrix(int nbStrikes,
                                    int nbDates) const;

    // registered fields
    double              strikeRef;
    double              upperSigmaUnscaled;

    // transient fields (won't appear in dd interface)
    DateTime            baseDate;
    TimeMetricSP        timeMetric;
    DateTimeArraySP     backboneDates;
    int                 nbBackboneDates;
    DoubleArraySP       backboneStrikes;
    int                 nbBackboneStrikes;
    CDoubleMatrixSP     backboneVols;
    DoubleArraySP       lowerSlopes;
    DoubleArraySP       upperStrikes;
    DoubleArraySP       upperVols;
    DoubleArraySP       upperSlopes;
    DoubleArraySP       upperCurvatures;
    DoubleArraySP       strikeRefVols;
    DoubleArraySP       timeFracs;
    DoubleArraySP       sqrTimeFracs;
    NRSpline::InterpolantArray splines;
    mutable IntArray    splineStrikeIndices;
    mutable CDoubleMatrixSP volMatrix;

    static const double default_upperSigmaUnscaled;
    static const double exphalf;

    VolSpline();

    static IObject* defaultCtor();
};

DRLIB_END_NAMESPACE

#endif