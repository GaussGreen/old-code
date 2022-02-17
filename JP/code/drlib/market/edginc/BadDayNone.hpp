//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : BadDayNone.hpp
//
//   Description : "None" bad day convention 
//
//   Author      : Andrew J Swain
//
//   Date        : 29 January 2001
//
//
//----------------------------------------------------------------------------

#ifndef BADDAYNONE_HPP
#define BADDAYNONE_HPP

#include "edginc/BadDayConvention.hpp"
#include <string>

using namespace std;   // string

DRLIB_BEGIN_NAMESPACE
    
/** "None" bad day convention */

class MARKET_DLL BadDayNone:public BadDayConvention {
public:
    static CClassConstSP const TYPE;
    friend class BadDayNoneHelper;

    BadDayNone();
    virtual ~BadDayNone();

    /** adjust date for bad days */
    virtual DateTime adjust(const DateTime &date, const Holiday *hols) const;

    /** returns a string description */
    virtual string toString() const;
};

typedef smartConstPtr<BadDayNone> BadDayNoneConstSP;
typedef smartPtr<BadDayNone> BadDayNoneP;

DRLIB_END_NAMESPACE
#endif
