//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : BadDayConvention.hpp
//
//   Description : Bad day convention interface
//
//   Author      : Andrew J Swain
//
//   Date        : 29 January 2001
//
//
//----------------------------------------------------------------------------

#ifndef BADDAYCONVENTION_HPP
#define BADDAYCONVENTION_HPP

#include "edginc/DateTime.hpp"
#include "edginc/Holiday.hpp"
#include "edginc/config.hpp"
#include <string>

using namespace std;   // string

DRLIB_BEGIN_NAMESPACE

class MaturityPeriod;


/** defines an interface to be implemented by concrete classes */
/** A BadDayConvention determines how to adjust a date for holidays */     
    
class MARKET_DLL BadDayConvention:public CObject {
public:
    static CClassConstSP const TYPE;
    friend class BadDayConventionHelper;

    virtual ~BadDayConvention();

    /** adjust date for bad days */
    virtual DateTime adjust(const DateTime &date,
                            const Holiday *hols) const = 0;

    /** Bad day adjust an array of dates in place. Default implementation just
        does a simple loop */
    virtual void adjust(vector<DateTime>::iterator start, 
                        vector<DateTime>::iterator end,
                        const Holiday*             hols) const;

    /** returns a string description e.g. Following */
    virtual string toString() const = 0;

protected:
    BadDayConvention(CClassConstSP clazz);

};

typedef smartConstPtr<BadDayConvention> BadDayConventionConstSP;
typedef smartPtr<BadDayConvention> BadDayConventionSP;
#ifndef QLIB_BADDAYCONVENTION_CPP
EXTERN_TEMPLATE(class MARKET_DLL_SP smartConstPtr<BadDayConvention>);
EXTERN_TEMPLATE(class MARKET_DLL_SP smartPtr<BadDayConvention>);
#else
INSTANTIATE_TEMPLATE(class MARKET_DLL smartConstPtr<BadDayConvention>);
INSTANTIATE_TEMPLATE(class MARKET_DLL smartPtr<BadDayConvention>);
#endif


/**
 * Utility function for date adjustment.
 */
DateTime DateFwdThenAdjust(
    const DateTime&         date, 
    const string&           tenor, 
    int                     count, 
    const BadDayConvention& bdc, 
    const Holiday&          holidays);


/**
 * Utility function for date adjustment.
 */
DateTime DateFwdThenAdjust(
    const DateTime&         date, 
    const MaturityPeriod&   tenor, 
    int                     count, 
    const BadDayConvention& bdc, 
    const Holiday&          holidays);


DRLIB_END_NAMESPACE
#endif
