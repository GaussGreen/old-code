//----------------------------------------------------------------------------
//
//   Group       : Credit QR
//
//   Filename    : IABCDSDecretion.hpp
//
//   Description : Interface to the ABCDS decretion information
//
//   Author      : Keith Jia
//
//   Date        : 03 January 2006
//
//
//----------------------------------------------------------------------------

#ifndef IABCDSDECRETION_HPP
#define IABCDSDECRETION_HPP

#include "edginc/Object.hpp"
#include "edginc/FORWARD_DECLARE_WRAPPER.hpp"

DRLIB_BEGIN_NAMESPACE

FORWARD_DECLARE_WRAPPER(IDecretionCurve);
FORWARD_DECLARE_WRAPPER(DecretionCurve);


/** Interface to expose Asset Backed CDS DecretionCurve (IDecretionCurve).
 */

class MARKET_DLL IABCDSDecretion : public virtual IObject {

public:
    virtual ~IABCDSDecretion();
    static CClassConstSP const TYPE;

    /** Currently we don't have a good model for loss. */
    //virtual IDecretionCurveSP getLossCurve() const = 0;

    /** Return prepayment curve. */
    virtual IDecretionCurveConstSP getPrepayCurve() const = 0;
    virtual DecretionCurveConstSP getPrepayCurveObj() const = 0;
};


DRLIB_END_NAMESPACE

#endif
