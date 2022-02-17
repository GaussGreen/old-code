//----------------------------------------------------------------------------
//
//   Group       : Quantitative Research
//
//   Filename    : SkewSurfaceBootstrapper.hpp
//
//   Description : Class to bootstrap a SkewSurface.
//
//   Author      : Matthias Arnsdorf
//
//   Date        : December 2006
//
//----------------------------------------------------------------------------

#ifndef QLIB_SKEW_SURFACE_BOOTSTRAPPER_HPP
#define QLIB_SKEW_SURFACE_BOOTSTRAPPER_HPP

#include "edginc/SkewSurface.hpp"
#include "edginc/TrancheIndexLeastSquareFit.hpp"
#include "edginc/DECLARE.hpp"


DRLIB_BEGIN_NAMESPACE


class PRODUCTS_DLL SkewSurfaceBootstrapper:
    public CObject,
    public virtual IBootstrapper
{
public:

    // Possible bootstrap types (used in TrancheIndexLeastSquarFit)
    /** bootstrap in time and strike dimension */
    static const string BOOTSTRAP_TIME_STRIKE;

    /** bootstrap in time dimension only (globcal calibration in strike dim) */
    static const string BOOTSTRAP_TIME;


    /** TYPE (for reflection) */
    static CClassConstSP const TYPE;

    /** Virtual destructor */
    virtual ~SkewSurfaceBootstrapper();

    /** Constructor (external) */
    SkewSurfaceBootstrapper(
        TrancheIndexLeastSquareFitSP  objFunc,
        const Calibrator::InstanceIDArray & ids);
  

    /** Main method that runs bootstrap calibration and returns results */
    void bootstrap(
        OptimizerNDSP optimizer,
        int & nbVars,                                 // (O) number of variables calibrated
        Calibrator::InstanceIDArray & aggregIds,    // (0)
        DoubleArray & aggregVals,                   // (O) calibrated values
        DoubleArray & obFuncValues                  // (O) objective function values
        );


private:
    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz);

    /** Constructor (internal - used by reflection) */
    SkewSurfaceBootstrapper();

    /** Default constructor */
    static IObject* defaultConstructor();


    // ---------
    // Methods
    // --------
    /**
    * Method called before first step of the loop
    * */
    void init();

    /**
    * Method called after each step of the loop
    * */
    void next(); 

    /**
    * Method called to test the end of the loop
    * */
    bool end() const;

    /**
    * Returns "state" corresponding to current step of the loop
    * */
    IObjectSP getCurrentState() const;

    /** get the instance ids to calibrate at current point in bootstrapping */
    Calibrator::InstanceIDArraySP getCurrentInstanceIDs() const;

    

    // ------
    // Fields
    // ------

    /** type of bootstrapping */
    string bootstrapType;

    /** maturities for bootstrapping */
    DateTimeArraySP bootstrapMaturities;

    /** Pointer to array of (maturity, strike1, strike2) points */
    // needed currently for validation
    TimePoint2DArraySP timePoints;

    /** Current index for bootstrapping*/
    int currBootstrapIdx;

    /** objective function */
    TrancheIndexLeastSquareFitSP  objFunc;

    /**
    * Array containing the index of the points in SkewSurface for
    * each bootstrap maturity
    * */
    IntArrayArray pointsToCalibrate;

    /** iterator for CDO quotes */
    CDOQuotesBootstrapperSP cdoQuotesIterator;

    /** Name of corresponding SkewSurface */
    string name;

    /** name of field to calibrate */
    string fieldToCalibrate; 

    /** override for range to calibrate */
    RangeSP rangeOverride;
};

DECLARE(SkewSurfaceBootstrapper);

DRLIB_END_NAMESPACE

#endif
