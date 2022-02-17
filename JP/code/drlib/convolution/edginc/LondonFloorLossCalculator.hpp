//----------------------------------------------------------------------------
//
//   Group       : QR
//
//   Author      : Mark A Robson
//
//   Date        : 23 Dec 2005
//
//----------------------------------------------------------------------------

#ifndef QR_LONDONFLOORLOSSCALCULATOR_HPP
#define QR_LONDONFLOORLOSSCALCULATOR_HPP

#include "edginc/TrancheLossCalculator.hpp"
#include "edginc/DateTime.hpp"
#include "edginc/FORWARD_DECLARE.hpp"
#include "edginc/Control_forward.hpp"
#include "edginc/Results_forward.hpp"

DRLIB_BEGIN_NAMESPACE

FORWARD_DECLARE(ICDSParSpreads);

/** Implemenetation of ITrancheLossCalculator using 'London Floor' adjustment
    on top of an existing ITrancheLossCalculator */
class CONVOLUTION_DLL LondonFloorLossCalculator: public ITrancheLossCalculator{
public:
    virtual ~LondonFloorLossCalculator();

    /** Constructor - takes in the loss calculator which does the unadjusted  
        loss. */
    LondonFloorLossCalculator(
        ITrancheLossCalculatorConstSP originalLossCalculator,     /* (I) */
        const DateTimeArray&          timeline,                   /* (I) */
        const DateTime&               today,                      /* (I) */
        bool                          computeCondCurve,           /* (I) */
        ICDSParSpreadsConstSP         londonFloor,                /* (I) */
        double                        portfolioNotional,          /* (I) */
        double                        equityStrike,               /* (I) */
        double                        seniorStrike);              /* (I) */

    /** Calculate the expected loss for specified timepoint and strikes. */
    virtual void loss(
        int     timePoint,        // (I) do the calculation for this timepoint
        double  k1,               /* (I) lower strike      */
        double  k2,               /* (I) upper strike      */
        double& loss,             /* (O) tranche loss amt  */
        double& lossCond) const;  /* (O) tranche loss amt cond on cpty
                                     surviving */

    /** store calculator result */
    void storeResults(Results* result, Control* control) const {};

    /** Returns a key used to optimise repeated calculations of losses
        at the same timepoint. */
    virtual IKey* lossKey(
        int timePoint) const; // (I) do the calculation for this timepoint
private:
    class Key;
    friend class Key;
    void loss(
        IKeySP  originalKey,          // (I) access to original calculator
        int     timePoint,            // (I) for accessing londonFloorProb
        double  k1,                   /* (I) lower strike      */
        double  k2,                   /* (I) upper strike      */
        double& loss,                 /* (O) tranche loss amt  */
        double& lossCond) const;      /* (O) tranche loss amt cond on cpty
                                         surviving */
    /// fields ///////
    ITrancheLossCalculatorConstSP originalLossCalculator;
    bool                          computeCondCurve;
    DoubleArray                   londonFloorProb;
    double                        totalNotional;
    double                        equityStrike;   
    double                        seniorStrike;
};
typedef smartConstPtr<
    LondonFloorLossCalculator> LondonFloorLossCalculatorConstSP;

DRLIB_END_NAMESPACE
#endif
