//----------------------------------------------------------------------------
//
//   Group       : Quantitative Research
//
//   Filename    : BSLPBootstrapper.hpp
//
//   Description : Class to bootstrap a SkewSurface for BSLP model.
//
//   Author      : Matthias Arnsdorf
//
//   Date        : December 2006
//
//----------------------------------------------------------------------------

#ifndef QLIB_BSLP_BOOTSTRAPPER_HPP
#define QLIB_BSLP_BOOTSTRAPPER_HPP

#include "edginc/SkewSurface.hpp"
#include "edginc/TrancheIndexLeastSquareFit.hpp"
#include "edginc/DECLARE.hpp"


DRLIB_BEGIN_NAMESPACE


class PRODUCTS_DLL BSLPBootstrapper:
    public CObject,
    public virtual IBootstrapper
{
public:

    /** botstrap type corresponding to this algorithm */
    static const string BOOTSTRAP_BSLP;

    /** TYPE (for reflection) */
    static CClassConstSP const TYPE;

    /** Virtual destructor */
    virtual ~BSLPBootstrapper();

    /** Constructor (external) */
    BSLPBootstrapper(
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
    BSLPBootstrapper();

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

DECLARE(BSLPBootstrapper);

DRLIB_END_NAMESPACE

#endif
