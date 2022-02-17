//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : FwdLevel.hpp
//
//   Description : spot level scenario - set spot to supplied value
//                 at a future date
//
//   Author      : Andrew J Swain
//
//   Date        : 9 June 2003
//
//
//----------------------------------------------------------------------------

#ifndef FWDLEVEL__HPP
#define FWDLEVEL__HPP
#include "edginc/SpotLevel.hpp"
#include "edginc/DateTime.hpp"

DRLIB_BEGIN_NAMESPACE

/** spot level scenario - set spot to supplied value at a future date */
class RISKMGR_DLL FwdLevel: public SpotLevel {
public:
    static CClassConstSP const TYPE;

    /** constructor with explicit spot level & date */
    FwdLevel(double spot, const DateTime& fwdDate);
 
    // given today's date, when is spot defined
    virtual DateTime spotDate(const DateTime& today) const;

private:
    friend class FwdLevelHelper;
    /** for reflection */
    FwdLevel();
    FwdLevel(const FwdLevel &rhs);
    FwdLevel& operator=(const FwdLevel& rhs);

    DateTime fwdDate;
};


typedef smartConstPtr<FwdLevel> FwdLevelConstSP;
typedef smartPtr<FwdLevel> FwdLevelSP;

DRLIB_END_NAMESPACE

#endif
