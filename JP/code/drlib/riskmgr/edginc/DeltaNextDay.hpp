//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : DeltaNextDay.hpp
//
//   Description : T+1 delta shift
//
//   Author      : Andrew J Swain
//
//   Date        : 16 May 2001
//
//
//----------------------------------------------------------------------------

#ifndef DELTANEXTDAY_HPP
#define DELTANEXTDAY_HPP

#include "edginc/T1Delta.hpp"

DRLIB_BEGIN_NAMESPACE

/** T+1 delta shift */
class RISKMGR_DLL DeltaNextDay: public T1Delta {
public:
    static CClassConstSP const TYPE;
    const static string NAME;
    const static double DEFAULT_SHIFT;
    const static int    DEFAULT_OFFSET;

    /** constructor */
    DeltaNextDay(double shift, int offset, HolidaySP hols);

    /** identifies the name used storing associated results in the output 
        (Implements pure virtual function in Sensitivity class) */
    const string& getSensOutputName() const;

    /** From INextDaySensitivity interface */
    virtual SensitivitySP requiredSensitivity(TweakGroup* tweakGroup) const;

private:
    /** for reflection */
    DeltaNextDay();
    DeltaNextDay(const DeltaNextDay &rhs);
    DeltaNextDay& operator=(const DeltaNextDay& rhs);
    friend class DeltaNextDayHelper;
};

typedef smartConstPtr<DeltaNextDay> DeltaNextDayConstSP;
typedef smartPtr<DeltaNextDay> DeltaNextDaySP;

DRLIB_END_NAMESPACE

#endif
