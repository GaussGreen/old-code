//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : BadDayPrevious.hpp
//
//   Description : "Previous" bad day convention 
//
//   Author      : Andrew J Swain
//
//   Date        : 29 January 2001
//
//
//----------------------------------------------------------------------------

#ifndef BADDAYPREVIOUS_HPP
#define BADDAYPREVIOUS_HPP

#include "edginc/BadDayConvention.hpp"
#include <string>

using namespace std;   // string

DRLIB_BEGIN_NAMESPACE
 
/**  "Previous" bad day convention */
  
class MARKET_DLL BadDayPrevious:public BadDayConvention {
public:
    static CClassConstSP const TYPE;
    friend class BadDayPreviousHelper;

    BadDayPrevious();
    virtual ~BadDayPrevious();

    /** adjust date for bad days */
    virtual DateTime adjust(const DateTime &date, const Holiday *hols) const;

    /** returns a string description */
    virtual string toString() const;
};

typedef smartConstPtr<BadDayPrevious> BadDayPreviousConstSP;
typedef smartPtr<BadDayPrevious> BadDayPreviousP;

DRLIB_END_NAMESPACE
#endif
