//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : DeltaProxyNextDay.hpp
//
//   Description : Fund proxy T+1 delta sensitivity
//
//   Author      : Andrew J Swain
//
//   Date        : 8 January 2001
//
//
//----------------------------------------------------------------------------

#ifndef _DELTAPROXYNEXTDAY_HPP
#define _DELTAPROXYNEXTDAY_HPP

#include "edginc/T1Delta.hpp"

DRLIB_BEGIN_NAMESPACE

/** T+1 delta shift for proxies */
class RISKMGR_DLL DeltaProxyNextDay: public T1Delta {
public:
    static CClassConstSP const TYPE;
    const static string NAME;
    const static double DEFAULT_SHIFT;
    const static int    DEFAULT_OFFSET;

    /** constructor */
    DeltaProxyNextDay(double shift, int offset, HolidaySP hols);

    /** identifies the name used storing associated results in the output 
        (Implements pure virtual function in Sensitivity class) */
    const string& getSensOutputName() const;

    /** From INextDaySensitivity interface */
    virtual SensitivitySP requiredSensitivity(TweakGroup* tweakGroup) const;
  
    /** From INextDaySensitivity interface */
    virtual void writeResults(const SensitivitySP& reqdSens,
                              Results*             dest,
                              const UntweakableSP& untweakable,
                              const Results*       src) const;
private:
    /** for reflection */
    DeltaProxyNextDay();
    DeltaProxyNextDay(const DeltaProxyNextDay &rhs);
    DeltaProxyNextDay& operator=(const DeltaProxyNextDay& rhs);
    friend class DeltaProxyNextDayHelper;
};

typedef smartConstPtr<DeltaProxyNextDay> DeltaProxyNextDayConstSP;
typedef smartPtr<DeltaProxyNextDay> DeltaProxyNextDaySP;

DRLIB_END_NAMESPACE

#endif
