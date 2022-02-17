//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : ThetaNIE.hpp
//
//   Description : Theta Net Interest Earned
//
//   Author      : Stephen Hope
//
//   Date        : 8 Nov 2001
//
//
//----------------------------------------------------------------------------


#ifndef EDR_THETA_NIE_HPP
#define EDR_THETA_NIE_HPP
#include "edginc/Theta.hpp"

DRLIB_BEGIN_NAMESPACE

/** Theta shift */
class RISKMGR_DLL ThetaNIE: public Theta{
public:
    static CClassConstSP const TYPE;
    static const string NAME;

    /** If timeOffsetNIE = END_OF_DAY_TIME return same date at EOD,
        If timeOffsetNIE = PRE_START_OF_DAY_TIME return date + offset @ PSOD */
    virtual DateTime rollDate(const DateTime &date) const;

    /** calculates given sensitivity - invoked by calculateSens */
    virtual void calculate(TweakGroup*  tweakGroup,
                           CResults*    results);

    /** Public constructor */
    ThetaNIE(int offset, HolidaySP hols);

private:
    friend class ThetaHelper;

    int            timeOffsetNIE;     // only used for ThetaNIE

    /** for reflection */
    ThetaNIE();
    ThetaNIE(const ThetaNIE &rhs);
    ThetaNIE& operator=(const ThetaNIE& rhs);

    static void load(CClassSP& clazz);
    static IObject* defaultThetaNIE();
};

typedef smartConstPtr<ThetaNIE> ThetaNIEConstSP;
typedef smartPtr<ThetaNIE> ThetaNIESP;


DRLIB_END_NAMESPACE

#endif
