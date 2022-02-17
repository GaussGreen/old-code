//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : MultiExpiryShift.hpp
//
//   Description : Specialised Perturbation where shifts for a range of 
//                 expiries are defined together in buckets. The interpretation
//                 of the buckets is not defined at this level - up to any
//                 derived classes. Really meant as a scenario not as a greek.
//
//   Author      : Andrew J Swain
//
//   Date        : 24 April 2003
//
//
//----------------------------------------------------------------------------

#ifndef _MULTIEXPIRYSHIFT__HPP
#define _MULTIEXPIRYSHIFT__HPP
#include "edginc/TweakNameListID.hpp"
#include "edginc/IRiskProperty.hpp"
#include "edginc/Perturbation.hpp"
#include "edginc/DateTime.hpp"
#include "edginc/Expiry.hpp"

DRLIB_BEGIN_NAMESPACE

/** Specialised Perturbation where shifts for a range of expiries are defined 
together in buckets. The interpretation of the buckets is not defined at this
level - up to any derived classes. Really meant as a scenario not as a greek.
*/
class RISKMGR_DLL MultiExpiryShift: public Perturbation,
                                    public virtual ITweakNameListID {
public:    
    static CClassConstSP const TYPE;

    virtual ~MultiExpiryShift();

    // what's the shift for a given date ? */
    virtual double shiftSize(const DateTime& today,     // to anchor expiries
                             const DateTime& shiftDate) const = 0;

    virtual OutputNameArrayConstSP names(IObjectConstSP tweakGroup) const;

    virtual IScalarRiskPropertyConstSP asRiskProperty() const;

    virtual void validatePop2Object();

protected:
    /** Note MultiExpiryShift is abstract. */
    MultiExpiryShift(const CClassConstSP& clazz);

    // fields
    ExpiryArray expiries;
    DoubleArray shifts;

private:
    OutputNameArraySP toTweak; /* to be removed - ignored. Legacy implementation
                                  reasons */
    friend class MultiExpiryShiftHelper;
};


typedef smartConstPtr<MultiExpiryShift> MultiExpiryShiftConstSP;
typedef smartPtr<MultiExpiryShift> MultiExpiryShiftSP;

DRLIB_END_NAMESPACE

#endif
