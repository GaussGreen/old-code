//----------------------------------------------------------------------------
//
//   Group       : QR
//
//   Filename    : RecoveryPerturb.hpp
//
//   Author      : Mark A Robson
//
//   Date        : 26 July 2005
//
//
//----------------------------------------------------------------------------


#ifndef EDG_RECOVERYSHIFT_HPP
#define EDG_RECOVERYSHIFT_HPP
#include "edginc/ScalarShift.hpp"
#include "edginc/Additive.hpp"

DRLIB_BEGIN_NAMESPACE

/** What the shift/restore methods for changing the [CDS] recovery take as a
    parameter */
class RISKMGR_DLL IRecoveryPerturb {
public:
    virtual ~IRecoveryPerturb(){}

    /** Returns the new recovery level given the original one */
    virtual double applyShift(double unadjRecovery) = 0;

    /** Returns the original recovery level given the adjusted one */
    virtual double undoShift(double adjRecovery) = 0;
};


DRLIB_END_NAMESPACE

#endif
