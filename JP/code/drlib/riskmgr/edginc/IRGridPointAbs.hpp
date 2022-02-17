//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : IRGridPoint.hpp
//
//   Description : Identifies a point on an IR vol surface
//
//   Author      : Andrew J Swain
//
//   Date        : 22 February 2002
//
//
//----------------------------------------------------------------------------

#ifndef _IRGRIDPOINT_ABS_HPP
#define _IRGRIDPOINT_ABS_HPP

#include "edginc/IRGridPoint.hpp"

DRLIB_BEGIN_NAMESPACE
class IRGridPointAbs;
typedef smartConstPtr<IRGridPointAbs> IRGridPointAbsConstSP;
typedef smartPtr<IRGridPointAbs> IRGridPointAbsSP;

typedef vector<IRGridPointAbsSP> IRGridPointAbsArray;
typedef refCountPtr<IRGridPointAbsArray> IRGridPointAbsArraySP;
typedef refCountPtr<const IRGridPointAbsArray> IRGridPointAbsArrayConstSP;

/** Identifies a point on an IR vol surface absolutely ie includes DateTime
    of expiry and swap maturity - this is useful when computing 'endDates'
    used to determine when to stop tweaking since mapping the expiries and
    tenors to dates can be non-trivial */
class RISKMGR_DLL IRGridPointAbs : public virtual VirtualDestructorBase {
public:

    ~IRGridPointAbs();

    IRGridPointAbs(ExpiryConstSP expiry, const DateTime& swaptionEnd,
                   ExpiryConstSP tenor,  const DateTime& swapMaturity);

    /** Returns an IRGridPoint from this object (IRGridPoint is the type
        returned to the client) */
    IRGridPointSP toGridPoint() const;

    /** returns the expiry associated with this result */
    ExpiryConstSP getExpiry() const;

    /** returns the tenor associated with this result */
    ExpiryConstSP getTenor() const;

    /** Modify gridPts so that all points whose expiries are after endDate,
        except for the very first, are removed. This is useful when determining
        when to stop tweaking */
    static void trim(IRGridPointAbsArray& gridPts, /* (M) */
                     const DateTime&      endDate);

    /** Finds the greatest swapMaturity date in the supplied list of points. The
        supplied array must not be empty */
    static DateTime maxSwapMaturity(const IRGridPointAbsArray& gridPts);

private:
    ExpiryConstSP expiry;
    DateTime swaptionEnd;
    ExpiryConstSP tenor;
    DateTime swapMaturity;
};


DRLIB_END_NAMESPACE
#endif
