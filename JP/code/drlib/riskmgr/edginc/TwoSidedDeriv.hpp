//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : TwoSidedDeriv.hpp
//
//   Description : Interface for sensitivities that are calculated via a two
//                 sided tweaking algorithm
//
//   Author      : Mark A Robson
//
//   Date        : 11 July 2002
//
//
//----------------------------------------------------------------------------

#ifndef EDR_TWOSIDEDDERIV_HPP
#define EDR_TWOSIDEDDERIV_HPP

#include "edginc/ScalarShift.hpp"

DRLIB_BEGIN_NAMESPACE

/** Interface for sensitivities that are calculated via a two sided
    tweaking algorithm (eg Delta, FXCrossGamma, CrossGamma) */
class RISKMGR_DLL ITwoSidedDeriv {
public:
    /** Returns the shift(s) which have been made for the current pricing
        call */
    virtual ScalarShiftArray getComponentShifts() const = 0;

    virtual ~ITwoSidedDeriv(); // empty

protected:
    ITwoSidedDeriv(); // empty
};


DRLIB_END_NAMESPACE
#endif
