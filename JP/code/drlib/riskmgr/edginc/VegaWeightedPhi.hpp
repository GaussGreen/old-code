//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : VegaWeightedPhi.hpp
//
//   Description : Vega weighted phi sensitivity
//
//   Author      : Andrew McCleery
//
//   Date        : 14 Oct 2004
//
//
//----------------------------------------------------------------------------


#ifndef VEGA_WEIGHTED_PHI_HPP
#define VEGA_WEIGHTED_PHI_HPP
#include "edginc/Sensitivity.hpp"
#include "edginc/smartPtr.hpp"
#include "edginc/Additive.hpp"

DRLIB_BEGIN_NAMESPACE

class RISKMGR_DLL VegaWeightedPhi: public Sensitivity, 
                       public virtual Additive {
public:
    static CClassConstSP const TYPE;
    const static string NAME;
    const static double DEFAULT_SHIFT;

    /** constructor with explicit shift size */
    VegaWeightedPhi(double shiftSize);

    /** identifies the name used for storing associated results in the output*/
    const string& getSensOutputName() const;

    /** Is this sensitivity made using a discrete shift (ie a jump) or a
        an approximately continuous one (return value: false) */
    bool discreteShift() const;

    // Calculate vega-weighted phi (parallel phi weighted by vega parallel)
    // (questionable theoretical justification)
    void calculate(TweakGroup*      tweakGroup,
                   Results*         results);
private:
    double shiftSize;      // not used

    friend class VegaWeightedPhiHelper;
    /** for reflection */
    VegaWeightedPhi();
    VegaWeightedPhi(const VegaWeightedPhi &rhs);
    VegaWeightedPhi& operator=(const VegaWeightedPhi& rhs);
};

typedef smartConstPtr<VegaWeightedPhi> VegaWeightedPhiConstSP;
typedef smartPtr<VegaWeightedPhi> VegaWeightedPhiSP;

DRLIB_END_NAMESPACE

#endif
