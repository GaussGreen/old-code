//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : BadDayModified.hpp
//
//   Description : "Modified Following" bad day convention 
//
//   Author      : Andrew J Swain
//
//   Date        : 29 January 2001
//
//
//----------------------------------------------------------------------------

#ifndef BADDAYMODIFIED_HPP
#define BADDAYMODIFIED_HPP

#include "edginc/BadDayConvention.hpp"
#include <string>

using namespace std;   // string

DRLIB_BEGIN_NAMESPACE

/** "Modified Following" bad day convention */   

class MARKET_DLL BadDayModified:public BadDayConvention {
public:
    static CClassConstSP const TYPE;
    friend class BadDayModifiedHelper;

    BadDayModified();
    virtual ~BadDayModified();

    /** adjust date for bad days */
    virtual DateTime adjust(const DateTime &date, const Holiday *hols) const;

    /** returns a string description */
    virtual string toString() const;
};

typedef smartConstPtr<BadDayModified> BadDayModifiedConstSP;
typedef smartPtr<BadDayModified> BadDayModifiedP;

DRLIB_END_NAMESPACE
#endif
