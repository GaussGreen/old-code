//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : BadDayFollowing.hpp
//
//   Description : "Following" bad day convention 
//
//   Author      : Andrew J Swain
//
//   Date        : 29 January 2001
//
//
//----------------------------------------------------------------------------

#ifndef BADDAYFOLLOWING_HPP
#define BADDAYFOLLOWING_HPP

#include "edginc/BadDayConvention.hpp"
#include <string>

using namespace std;   // string

DRLIB_BEGIN_NAMESPACE
    
/** "Following" bad day convention */

class MARKET_DLL BadDayFollowing:public BadDayConvention {
public:
    static CClassConstSP const TYPE;
    friend class BadDayFollowingHelper;

    BadDayFollowing();
    virtual ~BadDayFollowing();

    /** adjust date for bad days */
    virtual DateTime adjust(const DateTime &date, const Holiday *hols) const;

    /** returns a string description */
    virtual string toString() const;
};

typedef smartConstPtr<BadDayFollowing> BadDayFollowingConstSP;
typedef smartPtr<BadDayFollowing> BadDayFollowingP;

DRLIB_END_NAMESPACE
#endif
