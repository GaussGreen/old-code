//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : DeltaToCreditNextDay.hpp
//
//   Description : T+1 delta to credit shift
//
//   Author      : André Segger
//
//   Date        : 25 October 2002
//
//
//----------------------------------------------------------------------------

#ifndef DELTA_TO_CREDIT_NEXTDAY_HPP
#define DELTA_TO_CREDIT_NEXTDAY_HPP

#include "edginc/T1Delta.hpp"

DRLIB_BEGIN_NAMESPACE

/** T+1 delta shift */
class RISKMGR_DLL DeltaToCreditNextDay: public T1Delta {
public:
    static CClassConstSP const TYPE;
    const static string NAME;
    const static double DEFAULT_SHIFT;
    const static int    DEFAULT_OFFSET;

    /** constructor */
    DeltaToCreditNextDay(double shift, int offset, HolidaySP hols);

    /** identifies the name used storing associated results in the output 
        (Implements pure virtual function in Sensitivity class) */
    const string& getSensOutputName() const;

    /** From INextDaySensitivity interface */
    virtual SensitivitySP requiredSensitivity(TweakGroup* tweakGroup) const;
      
    /** From INextDaySensitivity interface */
    virtual void writeResults(const SensitivitySP&          reqdSens,
                              Results*                      dest,
                              const UntweakableSP&          untweakable,
                              const Results*                src) const;
private:
    /** for reflection */
    DeltaToCreditNextDay();
    DeltaToCreditNextDay(const DeltaToCreditNextDay &rhs);
    DeltaToCreditNextDay& operator=(const DeltaToCreditNextDay& rhs);
    friend class DeltaToCreditNextDayHelper;
};

typedef smartConstPtr<DeltaToCreditNextDay> DeltaToCreditNextDayConstSP;
typedef smartPtr<DeltaToCreditNextDay> DeltaToCreditNextDaySP;

DRLIB_END_NAMESPACE

#endif
