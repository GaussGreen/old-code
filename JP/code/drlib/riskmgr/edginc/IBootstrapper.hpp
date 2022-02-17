//----------------------------------------------------------------------------
//
//   Group       : Quantitative Research
//
//   Filename    : IBootstrapper.hpp
//
//   Description : Generic interface for objects capable to produce a "state"
//                 for each step of a loop
//
//   Author      : Antoine Gregoire
//
//   Date        : January 2006
//
//----------------------------------------------------------------------------


#ifndef I_BOOTSTRAPPER_HPP
#define I_BOOTSTRAPPER_HPP

#include "edginc/Array.hpp"
#include "edginc/Calibrator.hpp"

DRLIB_BEGIN_NAMESPACE

// Forward declaration of class
class IBootstrapper;


// Support for smart pointers
typedef smartPtr<IBootstrapper> IBootstrapperSP;
typedef smartConstPtr<IBootstrapper> IBootstrapperConstSP;

// Support for arrays
typedef array<IBootstrapperSP, IBootstrapper> IBootstrapperArray;
typedef smartPtr<IBootstrapperArray> IBootstrapperArraySP;


/**
 * Generic interface for objects capable of 'bootstrapping'
 * Bootstrapping can any type of bespoke calibration algortihm
 * */
class RISKMGR_DLL IBootstrapper: public virtual IObject {
public:    
    /** Virtual destructor */
    virtual ~IBootstrapper();

    /** Main method that runs bootstrap calibration and returns results */
    virtual void bootstrap(
        OptimizerNDSP optimizer,               // (I) calibrator to be used for calibration
        int & nbVars,                                 // (O) number of variables calibrated
        Calibrator::InstanceIDArray & aggregIds,    // (0)
        DoubleArray & aggregVals,                   // (O) calibrated values
        DoubleArray & obFuncValues                  // (O) objective function values
        ) = 0;            

    
    /** TYPE */
    static CClassConstSP const TYPE;
};

DRLIB_END_NAMESPACE

#endif /*I_BOOTSTRAPPER_HPP*/

