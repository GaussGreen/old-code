//----------------------------------------------------------------------------
//
//   Group       : Quantitative Research
//
//   Filename    : RFLMixtureModelCalibrator.hpp
//
//   Description : Class to calibrate the RFLMixtureModel
//
//   Author      : Jakob Sidenius
//
//   Date        : December 2006
//
//----------------------------------------------------------------------------

#ifndef QLIB_RFLMIXTURE_MODEL_CALIBRATOR_HPP
#define QLIB_RFLMIXTURE_MODEL_CALIBRATOR_HPP

#include "edginc/RFLMixtureDefaultsModel.hpp"
#include "edginc/TrancheIndexLeastSquareFit.hpp"
#include "edginc/DECLARE.hpp"


DRLIB_BEGIN_NAMESPACE


class PRODUCTS_DLL RFLMixtureModelCalibrator:
    public CObject,
    public virtual IBootstrapper
{
public:

    // Possible bootstrap types (used in TrancheIndexLeastSquarFit)
    /** bootstrap in time and strike dimension */
//    static const string BOOTSTRAP_TIME_STRIKE;

    /** bootstrap in time dimension only (globcal calibration in strike dim) */
//    static const string BOOTSTRAP_TIME;


    /** TYPE (for reflection) */
    static CClassConstSP const TYPE;

    /** Constructor (external) */
    RFLMixtureModelCalibrator(
        TrancheIndexLeastSquareFitSP  objFunc,
        const Calibrator::InstanceIDArray & ids);
  

    /** Main method that runs bootstrap calibration and returns results */
    void bootstrap(
        OptimizerNDSP optimizer,
        int & nbVars,                                 // (O) number of variables calibrated
        Calibrator::InstanceIDArray & aggregIds,    // (0)
        DoubleArray & aggregVals,                   // (O) calibrated values
        DoubleArray & obFuncValues);                  // (O) objective function values

private:
    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz);

    /** Constructor (internal - used by reflection) */
    RFLMixtureModelCalibrator();

    /** Default constructor */
    static IObject* defaultConstructor();

    /**
    * Function called before first step of the loop
    * */
    void init();

    /**
    * Function called after each step of the loop
    * */
    void next(); 

    /**
    * Function called to test the end of the loop
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

    Calibrator::InstanceIDArraySP theIds;

    /** objective function */
    TrancheIndexLeastSquareFitSP  objFunc;

};

DECLARE(RFLMixtureModelCalibrator);

DRLIB_END_NAMESPACE

#endif
