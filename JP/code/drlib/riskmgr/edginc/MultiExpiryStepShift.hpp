//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : MultiExpiryStepShift.hpp
//
//   Description : Specialised Perturbation where shifts for a range of 
//                 expiries are defined together in buckets. 
//                 Really meant as a scenario not as a greek.
//                 Shifts are defined for dates <= corresponding expiry
//                 in a step-wise manner.
//
//   Author      : Andrew J Swain
//
//   Date        : 24 April 2003
//
//
//----------------------------------------------------------------------------

#ifndef _MULTIEXPIRYSTEPSHIFT__HPP
#define _MULTIEXPIRYSTEPSHIFT__HPP
#include "edginc/MultiExpiryShift.hpp"

DRLIB_BEGIN_NAMESPACE

/** Specialised Perturbation where shifts for a range of expiries are defined 
together in buckets.  Really meant as a scenario not as a greek.
Shifts are defined for dates <= corresponding expiry in a step-wise manner.
*/

class RISKMGR_DLL MultiExpiryStepShift: public MultiExpiryShift {
public:    
    static CClassConstSP const TYPE;

    virtual ~MultiExpiryStepShift();

    // what's the shift for a given date ? */
    virtual double shiftSize(const DateTime& today,     // to anchor expiries
                             const DateTime& shiftDate) const;

protected:
    /** Note MultiExpiryStepShift is abstract. */
    MultiExpiryStepShift(const CClassConstSP& clazz);

private:
    friend class MultiExpiryStepShiftHelper;
    double upperShift;   // used after last expiry
};


typedef smartConstPtr<MultiExpiryStepShift> MultiExpiryStepShiftConstSP;
typedef smartPtr<MultiExpiryStepShift> MultiExpiryStepShiftSP;

DRLIB_END_NAMESPACE

#endif
