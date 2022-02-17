//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : VolParamTweak.hpp
//
//   Description : Interface for Tweaked Parametrized Vol Surface
//
//   Author      : Jean-Noël Juston
//
//   Date        : 25 Oct 2001
//
//
//----------------------------------------------------------------------------

#ifndef EDG_VOL_PARAM_TWEAK_HPP
#define EDG_VOL_PARAM_TWEAK_HPP

#include "edginc/VolParam.hpp"

#include "edginc/VolParallel.hpp"
#include "edginc/VegaMatrix.hpp"

DRLIB_BEGIN_NAMESPACE

/** The instances of CVolParamTweak are implementing each of the Vol tweaks
    on top of VolParam.
    They basically add a functional form (the tweak) to the VolParam they wrap.
    Still, you could do other things than just adding. */
class MARKET_DLL CVolParamTweak: public CVolParam{
public:
    static CClassConstSP const TYPE;

    /** constructor: takes clone of shift (casts to IObjcect), but not
        of paramVol. Can be removed once all the derived classes are
        written */
    CVolParamTweak(const CVolParamSP&         paramVol,
                   const DateTime&            valueDate,
                   const TimeMetricConstSP&   timeMetric,
                   const ITweakID*            shift);

    /** Implementation of method defined in CVolParam */
    virtual VolSurface* spotVolSurfaceFromStrikes(
        const CVolBase*       vol,
        const CDoubleArray&   strikes) const;

    /** Implementation of method defined in CVolParam.  */
    virtual void ComputeImpVol(const CVolBase*          vol,
                               const CLatticeDouble&    strikes,
                               const DateTimeArray&     maturities,
                               CLatticeDouble&          impV) const;

protected:
    CVolParamTweak(const CClassConstSP& clazz);

private:
    static void load(CClassSP& clazz);
    friend class TweakLattice;

    // fields //
    CVolParamConstSP      paramVol; // $unregistered
    DateTime              valueDate; // $unregistered
    TimeMetricConstSP     timeMetric; // $unregistered
    ITweakID*             theShift; // $unregistered
    IObjectSP             shiftSP; // ensure we delete it $unregistered
};


DRLIB_END_NAMESPACE
#endif //EDG_VOL_PARAM_TWEAK_HPP
