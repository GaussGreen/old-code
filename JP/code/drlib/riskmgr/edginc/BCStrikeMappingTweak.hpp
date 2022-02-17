//----------------------------------------------------------------------------
//
//   Group       : QR
//
//   Filename    : BCStrikeMappingTweak.hpp
//
//   Description : Tweak to shift/set the Base Correlation Strike Mapping
//
//   Author      : Jose Hilera
//
//   Date        : 11 October 2005
//
//----------------------------------------------------------------------------


#ifndef QLIB_BCSTRIKEMAPPINGTWEAK_HPP
#define QLIB_BCSTRIKEMAPPINGTWEAK_HPP

#include "edginc/BCStrikeMappingTweakBase.hpp"

DRLIB_BEGIN_NAMESPACE

/** Perturbation for strike mapping scenarios */
class RISKMGR_DLL BCStrikeMappingTweak : public BCStrikeMappingTweakBase {
public:
    static CClassConstSP const TYPE;

    /** Shifts the original StrikeMapping according to the tweakType and
        returns the new one */
    virtual double applyShift(double& unadjStrikeMapping);
  
private:
    string tweakType; // optional - default to "S"et (as opposed to 
                      // "A"bsolute or "R"elative shift)

    /** for reflection */
    static void load(CClassSP& clazz);
    static IObject* defaultConstructor();
    BCStrikeMappingTweak();
    BCStrikeMappingTweak(const BCStrikeMappingTweak&);
    BCStrikeMappingTweak& operator=(const BCStrikeMappingTweak&);
};


DRLIB_END_NAMESPACE

#endif
