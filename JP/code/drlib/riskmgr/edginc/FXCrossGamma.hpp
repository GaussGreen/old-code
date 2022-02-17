//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : CrossGamma.hpp
//
//   Description : CrossGamma sensitivity
//
//   Author      : Mark A Robson
//
//   Date        : 16 Mar 2001
//
//
//----------------------------------------------------------------------------

#ifndef EDG_FX_CROSS_GAMMA_HPP
#define EDG_FX_CROSS_GAMMA_HPP
#include "edginc/ScalarShift.hpp"
#include "edginc/Additive.hpp"
#include "edginc/TwoSidedDeriv.hpp"
#include "edginc/AssetSpotGreek.hpp"

DRLIB_BEGIN_NAMESPACE

/** Sens Control for CROSS_GAMMA - a second order cross
    derivative. The Delta::Shift interface and the
    FXDelta::Shiftdefines what needs to be implemented for CROSS_GAMMA
    to be calculated */
class RISKMGR_DLL FXCrossGamma: public Sensitivity,
                    virtual public Additive,
                    virtual public IAssetSpotGreek,
                    virtual public ITwoSidedDeriv {
public:
    friend class FXCrossGammaHelper;
    static CClassConstSP const TYPE;
    const static string NAME;
    const static double DEFAULT_SHIFT;

    /** constructor with explicit shift size */
    FXCrossGamma(double     shiftSize);

    /** Is this sensitivity made using a discrete shift (ie a jump) or a
        an approximately continuous one. Returns false */
    virtual bool discreteShift() const;

    /** calculates FXCrossGamma - override delta calculation  */
    virtual void calculate(TweakGroup*      tweakGroup,
                           Results*         results);

    /** identifies the name used for storing associated results in the output*/
    virtual const string& getSensOutputName() const;

    /** When calculating cross gamma, several pricings have to be done. These
        either involve a single shift or a double shift. This routine, called
        upon the Sensitivity returned from Control::getSens(), returns the
        shifts which have been made for the current pricing call */
    virtual ScalarShiftArray getComponentShifts() const;

private:
    /** for reflection */
    FXCrossGamma();
    FXCrossGamma(const FXCrossGamma &rhs);
    FXCrossGamma& operator=(const FXCrossGamma& rhs);

    // **** fields ****
    double shiftSize;
    // not registered
    ScalarShiftArray   shifts; // $unregistered

};

typedef smartConstPtr<FXCrossGamma> FXCrossGammaConstSP;
typedef smartPtr<FXCrossGamma> FXCrossGammaSP;

DRLIB_END_NAMESPACE

#endif
