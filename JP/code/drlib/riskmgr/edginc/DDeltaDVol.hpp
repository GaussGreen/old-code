//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : DDeltaDVol.Hpp
//
//   Description : DDeltaDVol sensitivity
//
//   Author      : Jay R Blumenstein
//
//   Date        : 15 Apr 2003
//
//
//----------------------------------------------------------------------------

#ifndef EDG_DDELTADVOL_HPP
#define EDG_DDELTADVOL_HPP

#include "edginc/ScalarShift.hpp"
#include "edginc/Additive.hpp"
#include "edginc/TwoSidedDeriv.hpp"

DRLIB_BEGIN_NAMESPACE

/** Sens Control for CROSS_GAMMA - a second order cross
    derivative. The Delta::Shift interface and the
    FXDelta::Shiftdefines what needs to be implemented for CROSS_GAMMA
    to be calculated */
class RISKMGR_DLL DDeltaDVol: public Sensitivity,
                  public virtual Additive,
                  public virtual ITwoSidedDeriv {
public:
    friend class DDeltaDVolHelper;
    static CClassConstSP const TYPE;
    const static string NAME;
    const static double DEFAULT_SHIFT;

    /** constructor with explicit shift size */
    DDeltaDVol(double     shiftSize);

    /** Is this sensitivity made using a discrete shift (ie a jump) or a
        an approximately continuous one. Returns false */
    virtual bool discreteShift() const;

    /** calculates DDeltaDVol - override delta calculation  */
    virtual void calculate(TweakGroup*      tweakGroup,
                           Results*         results);

    /** identifies the name used for storing associated results in the output*/
    virtual const string& getSensOutputName() const;

    /** When calculating cross gamma, several pricings have to be done. These
        either involve a single shift or a double shift. This routine, called
        upon the Sensitivity returned from Control::getSens(), returns the
        shifts which have been made for the current pricing call */
    virtual ScalarShiftArray getComponentShifts() const;

    /** class just to get the name pairs */

private:
    /** for reflection */
    DDeltaDVol();
    DDeltaDVol(const DDeltaDVol &rhs);
    DDeltaDVol& operator=(const DDeltaDVol& rhs);

    // **** fields ****
    double shiftSize;
    // not registered
    ScalarShiftArray   shifts; // $unregistered

};

typedef smartConstPtr<DDeltaDVol> DDeltaDVolConstSP;
typedef smartPtr<DDeltaDVol> DDeltaDVolSP;

DRLIB_END_NAMESPACE

#endif
