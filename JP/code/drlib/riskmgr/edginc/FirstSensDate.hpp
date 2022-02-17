//----------------------------------------------------------------------------
//
//   Group       : GCCG Derivatives Research
//
//   Filename    : FirstSensDate.hpp
//
//   Description : Interface for instruments that know when to begin
//                 tweaking pointwise greeks. Modeled after LastSensDate.hpp
//
//   Author      : Sean Chen
//
//   Date        : 5 Aug 2005
//
//
//----------------------------------------------------------------------------

#ifndef FIRSTSENSDATE_HPP
#define FIRSTSENSDATE_HPP

#include "edginc/Object.hpp"
#include "edginc/DateTime.hpp"
#include "edginc/SensControl.hpp"
#include "edginc/Model.hpp"
#include "edginc/Expiry.hpp"

DRLIB_BEGIN_NAMESPACE

/** Interface for instruments that know when to begin
    tweaking pointwise greeks
*/
class RISKMGR_DLL FirstSensDate {
public:
    static CClassConstSP const TYPE;

    virtual ~FirstSensDate();

    /** when to stop tweaking */
    virtual DateTime beginDate(const SensControl* sensControl) const = 0;

protected:
    FirstSensDate();
};


DRLIB_END_NAMESPACE
#endif
