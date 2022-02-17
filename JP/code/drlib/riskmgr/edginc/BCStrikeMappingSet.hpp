//----------------------------------------------------------------------------
//
//   Group       : QR
//
//   Filename    : BCStrikeMappingSet.hpp
//
//   Description : Tweak to shift/set the Base Correlation Strike Mapping
//
//   Author      : Jose Hilera
//
//   Date        : 11 October 2005
//
//----------------------------------------------------------------------------


#ifndef QLIB_BCSTRIKEMAPPINGSET_HPP
#define QLIB_BCSTRIKEMAPPINGSET_HPP

#include "edginc/BCStrikeMappingTweakBase.hpp"

DRLIB_BEGIN_NAMESPACE

/** Perturbation for strike mapping scenarios */
class RISKMGR_DLL BCStrikeMappingSet: public BCStrikeMappingTweakBase {
public:
    static CClassConstSP const TYPE;

    /** Shifts the original StrikeMapping according to the tweakType and
        returns the new one */
    virtual double applyShift(double& unadjStrikeMapping);
  
private:
    /** for reflection */
    static void load(CClassSP& clazz);
    static IObject* defaultConstructor();
    BCStrikeMappingSet();
    BCStrikeMappingSet(const BCStrikeMappingSet&);
    BCStrikeMappingSet& operator=(const BCStrikeMappingSet&);
};


DRLIB_END_NAMESPACE

#endif
