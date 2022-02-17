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

#ifndef _IRGRIDPOINT_HPP
#define _IRGRIDPOINT_HPP

#include "edginc/Expiry.hpp"

DRLIB_BEGIN_NAMESPACE

/** Identifies a point on an IR vol surface */
class RISKMGR_DLL IRGridPoint: public CObject {
public:
    static CClassConstSP const TYPE;

    virtual ~IRGridPoint();

    IRGridPoint(const ExpiryConstSP& tenor, const ExpiryConstSP& expiry);

    /** returns the expiry associated with this result */
    ExpiryConstSP getExpiry() const;

    /** returns the tenor associated with this result */
    ExpiryConstSP getTenor() const;

private:
    friend class IRGridPointHelper;
    IRGridPoint();

    ExpirySP tenor;
    ExpirySP expiry;
};

typedef smartConstPtr<IRGridPoint> IRGridPointConstSP;
typedef smartPtr<IRGridPoint> IRGridPointSP;

typedef array<IRGridPointSP, IRGridPoint> IRGridPointArray;
typedef smartPtr<IRGridPointArray> IRGridPointArraySP;
typedef smartConstPtr<IRGridPointArray> IRGridPointArrayConstSP;

DRLIB_END_NAMESPACE
#endif
