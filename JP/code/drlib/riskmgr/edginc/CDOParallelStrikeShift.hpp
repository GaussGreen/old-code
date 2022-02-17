//----------------------------------------------------------------------------
//
//   Group       : QR - Credit Hybrids
//
//   Filename    : CDOParallelStrikeShift.hpp
//
//   Description : Sensitivity to a parallel shift in CDO strikes
//
//   Author      : Linus Thand
//
//   Date        : 13 July 2006
//
//----------------------------------------------------------------------------

/* A sensitivity to a parallel shift in CDO strikes, used by managed trades
 * through ImpliedScalarShiftMulti.
 * Inheriting from AllNamesRiskPropertySensitivity, so it's a declarative greek.
 * Implementing IScalarPerNameShift, so it can be used for solving.
 */

#ifndef QLIB_CDO_PARALLEL_STRIKE_SHIFT_H
#define QLIB_CDO_PARALLEL_STRIKE_SHIFT_H

#include "edginc/CDOParallelStrike.hpp"
#include "edginc/AllNamesRiskPropertySensitivity.hpp"
#include "edginc/ScalarRiskPropertySensitivity.hpp"

DRLIB_BEGIN_NAMESPACE

class RISKMGR_DLL CDOParallelStrikeShift : public AllNamesRiskPropertySensitivity, 
                                           virtual public IScalarPerNameShift {
 public:
    static CClassConstSP const TYPE;
    static const double DEFAULT_SHIFT;
    CDOParallelStrikeShift(double shiftSize = DEFAULT_SHIFT);
        
    //// To implement IScalarPerNameShift
    IHypothesis::AlternateWorldSP appliedTo(OutputNameConstSP name,
                                            double shiftSize,
                                            IObjectSP world);
    OutputNameArrayConstSP allNames(const IObject* object) const;
    
 private:
    double getShiftSize() const;
    static const double SENSITIVITY_UNIT;
    static const string NAME;
    ~CDOParallelStrikeShift();
    CDOParallelStrikeShift(const CDOParallelStrikeShift& rhs); 
    CDOParallelStrikeShift& operator=(const CDOParallelStrikeShift& rhs);
    static void load(CClassSP& clazz);
    Deriv deriv() const;
};

FORWARD_DECLARE(CDOParallelStrikeShift)

DRLIB_END_NAMESPACE

#endif
